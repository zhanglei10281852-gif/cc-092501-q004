from __future__ import annotations

import hashlib
import json
import math
import random
import sqlite3
from datetime import datetime, timezone
from typing import Any

from app.database import get_connection, transaction


SCHEMA = """
CREATE TABLE IF NOT EXISTS hydro_wells (
 id INTEGER PRIMARY KEY AUTOINCREMENT, code TEXT NOT NULL UNIQUE, name TEXT NOT NULL,
 latitude REAL NOT NULL, longitude REAL NOT NULL, aquifer TEXT NOT NULL, screen_depth_m REAL NOT NULL,
 status TEXT NOT NULL DEFAULT 'active', created_at TEXT NOT NULL, updated_at TEXT NOT NULL
);
CREATE TABLE IF NOT EXISTS hydro_endmembers (
 id INTEGER PRIMARY KEY AUTOINCREMENT, name TEXT NOT NULL, isotope_d18o REAL NOT NULL,
 isotope_d2h REAL NOT NULL, solute_mg_l REAL NOT NULL, uncertainty REAL NOT NULL,
 version TEXT NOT NULL, active INTEGER NOT NULL DEFAULT 1 CHECK(active IN (0,1)), created_at TEXT NOT NULL,
 UNIQUE(name,version)
);
CREATE TABLE IF NOT EXISTS hydro_samples (
 id INTEGER PRIMARY KEY AUTOINCREMENT, well_id INTEGER NOT NULL REFERENCES hydro_wells(id) ON DELETE RESTRICT,
 sample_code TEXT NOT NULL UNIQUE, sampled_at TEXT NOT NULL, isotope_d18o REAL, isotope_d2h REAL,
 solute_mg_l REAL, detection_limit REAL NOT NULL, measurement_error REAL NOT NULL,
 quality_status TEXT NOT NULL DEFAULT 'pending', created_at TEXT NOT NULL
);
CREATE TABLE IF NOT EXISTS hydro_inversions (
 id INTEGER PRIMARY KEY AUTOINCREMENT, sample_id INTEGER NOT NULL REFERENCES hydro_samples(id) ON DELETE RESTRICT,
 task_key TEXT NOT NULL UNIQUE, model_version TEXT NOT NULL, method TEXT NOT NULL,
 input_json TEXT NOT NULL, status TEXT NOT NULL DEFAULT 'queued', attempts INTEGER NOT NULL DEFAULT 0,
 worker_id TEXT NOT NULL DEFAULT '', result_json TEXT NOT NULL DEFAULT '{}', error TEXT NOT NULL DEFAULT '',
 created_at TEXT NOT NULL, updated_at TEXT NOT NULL
);
CREATE TABLE IF NOT EXISTS hydro_transport_runs (
 id INTEGER PRIMARY KEY AUTOINCREMENT, well_id INTEGER NOT NULL REFERENCES hydro_wells(id) ON DELETE RESTRICT,
 task_key TEXT NOT NULL UNIQUE, model_version TEXT NOT NULL, input_json TEXT NOT NULL,
 status TEXT NOT NULL DEFAULT 'queued', result_json TEXT NOT NULL DEFAULT '{}', created_at TEXT NOT NULL, updated_at TEXT NOT NULL
);
CREATE TABLE IF NOT EXISTS hydro_interval_runs (
 id INTEGER PRIMARY KEY AUTOINCREMENT,
 inversion_id INTEGER NOT NULL REFERENCES hydro_inversions(id) ON DELETE RESTRICT,
 task_key TEXT NOT NULL UNIQUE,
 method TEXT NOT NULL CHECK(method IN ('parametric-bootstrap','deterministic-profile')),
 confidence_level REAL NOT NULL, replicates INTEGER NOT NULL, seed INTEGER,
 model_version TEXT NOT NULL, input_json TEXT NOT NULL,
 status TEXT NOT NULL DEFAULT 'queued' CHECK(status IN ('queued','running','done','failed')),
 attempts INTEGER NOT NULL DEFAULT 0, worker_id TEXT NOT NULL DEFAULT '',
 progress INTEGER NOT NULL DEFAULT 0, result_json TEXT NOT NULL DEFAULT '{}',
 error TEXT NOT NULL DEFAULT '', created_at TEXT NOT NULL, updated_at TEXT NOT NULL
);
CREATE TABLE IF NOT EXISTS hydro_interval_replicates (
 id INTEGER PRIMARY KEY AUTOINCREMENT,
 run_id INTEGER NOT NULL REFERENCES hydro_interval_runs(id) ON DELETE CASCADE,
 replicate_index INTEGER NOT NULL, fractions_json TEXT NOT NULL,
 objective REAL NOT NULL, converged INTEGER NOT NULL CHECK(converged IN (0,1)),
 created_at TEXT NOT NULL, UNIQUE(run_id,replicate_index)
);
CREATE TABLE IF NOT EXISTS hydro_audit (
 id INTEGER PRIMARY KEY AUTOINCREMENT, resource_type TEXT NOT NULL, resource_id INTEGER,
 action TEXT NOT NULL, actor TEXT NOT NULL, payload_json TEXT NOT NULL DEFAULT '{}', created_at TEXT NOT NULL
);
CREATE INDEX IF NOT EXISTS idx_hydro_samples_well ON hydro_samples(well_id,sampled_at);
CREATE INDEX IF NOT EXISTS idx_hydro_inversions_status ON hydro_inversions(status,created_at);
CREATE INDEX IF NOT EXISTS idx_hydro_interval_runs_inversion ON hydro_interval_runs(inversion_id,confidence_level,model_version);
CREATE INDEX IF NOT EXISTS idx_hydro_interval_replicates_run ON hydro_interval_replicates(run_id,replicate_index);
"""


def _now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def ensure_schema() -> None:
    get_connection().executescript(SCHEMA)


def _digest(value: Any) -> str:
    return hashlib.sha256(json.dumps(value, sort_keys=True, ensure_ascii=False).encode()).hexdigest()


def _dict(row: sqlite3.Row | None) -> dict[str, Any] | None:
    return dict(row) if row else None


def _tracer_scale(observed: list[float | None]) -> list[float]:
    return [20.0, 100.0, max(1.0, float(observed[2]) if observed[2] is not None else 1.0)]


def _seeded_rng(seed: int, index: int) -> random.Random:
    """每个重复样本使用独立确定的随机流，与执行顺序和分块方式无关。"""
    digest = hashlib.sha256(f"hydro-interval:{seed}:{index}".encode()).digest()
    return random.Random(int.from_bytes(digest[:8], "big"))


def _quantile(sorted_values: list[float], probability: float) -> float | None:
    if not sorted_values: return None
    count = len(sorted_values)
    if count == 1: return sorted_values[0]
    position = probability * (count - 1)
    lower = int(math.floor(position))
    upper = min(lower + 1, count - 1)
    fraction = position - lower
    return sorted_values[lower] * (1 - fraction) + sorted_values[upper] * fraction


def _pearson(xs: list[float], ys: list[float]) -> float | None:
    count = len(xs)
    if count < 3: return None
    mean_x = sum(xs) / count
    mean_y = sum(ys) / count
    sxx = sum((x - mean_x) ** 2 for x in xs)
    syy = sum((y - mean_y) ** 2 for y in ys)
    if sxx <= 1e-18 or syy <= 1e-18: return None
    return sum((x - mean_x) * (y - mean_y) for x, y in zip(xs, ys)) / math.sqrt(sxx * syy)


def _normal_ppf(probability: float) -> float:
    """Acklam 有理近似的标准正态分位数，最大误差约 1.15e-9，完全确定。"""
    a = [-3.969683028665376e+01, 2.209460984245205e+02, -2.759285104469687e+02, 1.383577518672690e+02, -3.066479806614716e+01, 2.506628277459239e+00]
    b = [-5.447609879822406e+01, 1.615858368580409e+02, -1.556989798598866e+02, 6.680131188771972e+01, -1.328068155288572e+01]
    c = [-7.784894002430293e-03, -3.223964580411365e-01, -2.400758277161838e+00, -2.549732539343734e+00, 4.374664141464968e+00, 2.938163982698783e+00]
    d = [7.784695709041462e-03, 3.224671290700398e-01, 2.445134137142996e+00, 3.754408661907416e+00]
    plow = 0.02425
    if probability < plow:
        q = math.sqrt(-2 * math.log(probability))
        return (((((c[0]*q+c[1])*q+c[2])*q+c[3])*q+c[4])*q+c[5]) / ((((d[0]*q+d[1])*q+d[2])*q+d[3])*q+1)
    if probability > 1 - plow:
        q = math.sqrt(-2 * math.log(1 - probability))
        return -(((((c[0]*q+c[1])*q+c[2])*q+c[3])*q+c[4])*q+c[5]) / ((((d[0]*q+d[1])*q+d[2])*q+d[3])*q+1)
    q = probability - 0.5
    r = q * q
    return (((((a[0]*r+a[1])*r+a[2])*r+a[3])*r+a[4])*r+a[5])*q / (((((b[0]*r+b[1])*r+b[2])*r+b[3])*r+b[4])*r+1)


def _chi2_df1(confidence_level: float) -> float:
    """单参数剖面似然阈值：chi2(1) 分位数 = 正态分位数的平方。"""
    z = _normal_ppf((1.0 + confidence_level) / 2.0)
    return z * z


class HydroService:
    def __init__(self, connection: sqlite3.Connection | None = None):
        self.connection = connection or get_connection()
        ensure_schema()

    def create_well(self, payload: dict[str, Any], actor: str = "researcher") -> dict[str, Any]:
        now = _now()
        with transaction(immediate=True) as connection:
            cursor = connection.execute("INSERT INTO hydro_wells(code,name,latitude,longitude,aquifer,screen_depth_m,created_at,updated_at) VALUES(?,?,?,?,?,?,?,?)", (payload["code"],payload["name"],payload["latitude"],payload["longitude"],payload["aquifer"],payload["screen_depth_m"],now,now))
            well_id = cursor.lastrowid
            connection.execute("INSERT INTO hydro_audit(resource_type,resource_id,action,actor,payload_json,created_at) VALUES('well',?,?,?,?,?)", (well_id,"create",actor,json.dumps(payload,ensure_ascii=False),now))
            return dict(connection.execute("SELECT * FROM hydro_wells WHERE id=?",(well_id,)).fetchone())

    def get_well(self, well_id: int) -> dict[str, Any] | None:
        well = self.connection.execute("SELECT * FROM hydro_wells WHERE id=?",(well_id,)).fetchone()
        if well is None: return None
        result = dict(well)
        result["samples"] = [dict(r) for r in self.connection.execute("SELECT * FROM hydro_samples WHERE well_id=? ORDER BY sampled_at,id",(well_id,)).fetchall()]
        return result

    def delete_well(self, well_id: int) -> bool:
        with transaction(immediate=True) as connection:
            cursor = connection.execute("DELETE FROM hydro_wells WHERE id=?",(well_id,))
            if cursor.rowcount == 0: raise KeyError("well_not_found")
            return True

    def create_endmember(self, payload: dict[str, Any]) -> dict[str, Any]:
        now=_now()
        with transaction(immediate=True) as connection:
            cursor=connection.execute("INSERT INTO hydro_endmembers(name,isotope_d18o,isotope_d2h,solute_mg_l,uncertainty,version,created_at) VALUES(?,?,?,?,?,?,?)",(payload["name"],payload["isotope_d18o"],payload["isotope_d2h"],payload["solute_mg_l"],payload["uncertainty"],payload["version"],now))
            return dict(connection.execute("SELECT * FROM hydro_endmembers WHERE id=?",(cursor.lastrowid,)).fetchone())

    def add_sample(self, well_id: int, payload: dict[str, Any]) -> dict[str, Any]:
        if self.connection.execute("SELECT id FROM hydro_wells WHERE id=?",(well_id,)).fetchone() is None: raise KeyError("well_not_found")
        values=[payload.get("isotope_d18o"),payload.get("isotope_d2h"),payload.get("solute_mg_l")]
        quality="usable" if sum(v is not None for v in values)>=2 else "incomplete"
        now=_now()
        with transaction(immediate=True) as connection:
            cursor=connection.execute("INSERT INTO hydro_samples(well_id,sample_code,sampled_at,isotope_d18o,isotope_d2h,solute_mg_l,detection_limit,measurement_error,quality_status,created_at) VALUES(?,?,?,?,?,?,?,?,?,?)",(well_id,payload["sample_code"],payload["sampled_at"],payload.get("isotope_d18o"),payload.get("isotope_d2h"),payload.get("solute_mg_l"),payload["detection_limit"],payload["measurement_error"],quality,now))
            return dict(connection.execute("SELECT * FROM hydro_samples WHERE id=?",(cursor.lastrowid,)).fetchone())

    def _project_simplex(self, values: list[float], fixed: dict[int, float] | None = None) -> list[float]:
        if not fixed:
            clipped=[max(0.0,v) for v in values]
            total=sum(clipped)
            return [1/len(values)]*len(values) if total<=1e-15 else [v/total for v in clipped]
        result=list(values)
        for index,value in fixed.items(): result[index]=min(1.0,max(0.0,value))
        pinned=sum(result[index] for index in fixed)
        free=[i for i in range(len(values)) if i not in fixed]
        rest=max(0.0,1.0-pinned)
        clipped=[max(0.0,result[i]) for i in free]
        total=sum(clipped)
        if total<=1e-15:
            share=rest/len(free) if free else 0.0
            for i in free: result[i]=share
        else:
            for i,value in zip(free,clipped): result[i]=value/total*rest
        return result

    def _objective(self, fractions: list[float], observed: list[float | None], vectors: list[list[float]], active: list[int], scale: list[float]) -> float:
        predicted=[sum(fractions[j]*vectors[j][k] for j in range(len(vectors))) for k in range(3)]
        return sum(((predicted[k]-float(observed[k]))/scale[k])**2 for k in active)+((sum(fractions)-1.0)*10)**2

    def _solve_fractions(self, observed: list[float | None], vectors: list[list[float]], max_iterations: int, tolerance: float, fixed: dict[int, float] | None = None) -> dict[str, Any]:
        """投影梯度求解端元比例；fixed 可固定部分端元（剖面法使用）。纯函数，无随机性。"""
        active=[i for i,v in enumerate(observed) if v is not None]
        if len(active)<2: raise ValueError("insufficient_measurements")
        count=len(vectors)
        fractions=[1/count]*count
        if fixed:
            rest=max(0.0,1.0-sum(min(1.0,max(0.0,v)) for v in fixed.values()))
            share=rest/(count-len(fixed)) if count>len(fixed) else 0.0
            for i in range(count): fractions[i]=min(1.0,max(0.0,fixed[i])) if i in fixed else share
        scale=_tracer_scale(observed)
        rate=0.08
        last=float("inf")
        objective=float("inf")
        iteration=-1
        for iteration in range(max_iterations):
            predicted=[sum(fractions[j]*vectors[j][k] for j in range(count)) for k in range(3)]
            residual=[(predicted[k]-float(observed[k]))/scale[k] if k in active else 0.0 for k in range(3)]
            objective=sum(r*r for r in residual)+((sum(fractions)-1.0)*10)**2
            if abs(last-objective)<tolerance: break
            last=objective
            gradient=[]
            for j in range(count):
                gradient.append(0.0 if fixed and j in fixed else 2*sum(residual[k]*vectors[j][k]/scale[k] for k in active))
            fractions=self._project_simplex([f-rate*g for f,g in zip(fractions,gradient)],fixed)
        predicted=[sum(fractions[j]*vectors[j][k] for j in range(count)) for k in range(3)]
        rmse=math.sqrt(sum(((predicted[k]-float(observed[k]))/scale[k])**2 for k in active)/len(active))
        return {"fractions":fractions,"predicted":predicted,"rmse":rmse,"iterations":iteration+1,"converged":abs(last-objective)<tolerance,"objective":self._objective(fractions,observed,vectors,active,scale)}

    def solve_mixture(self, sample: sqlite3.Row, endmembers: list[sqlite3.Row], max_iterations: int, tolerance: float) -> dict[str, Any]:
        observed=[sample["isotope_d18o"],sample["isotope_d2h"],sample["solute_mg_l"]]
        vectors=[[e["isotope_d18o"],e["isotope_d2h"],e["solute_mg_l"]] for e in endmembers]
        core=self._solve_fractions(observed,vectors,max_iterations,tolerance)
        fractions=core["fractions"]
        return {"fractions":[{"endmember_id":e["id"],"name":e["name"],"fraction":round(f,8)} for e,f in zip(endmembers,fractions)],"mass_balance":round(sum(fractions),10),"predicted":core["predicted"],"rmse":core["rmse"],"iterations":core["iterations"],"converged":core["converged"]}

    def enqueue_inversion(self, sample_id: int, payload: dict[str, Any]) -> dict[str, Any]:
        sample=self.connection.execute("SELECT * FROM hydro_samples WHERE id=?",(sample_id,)).fetchone()
        if sample is None: raise KeyError("sample_not_found")
        ids=sorted(set(payload["endmember_ids"]))
        endmembers=self.connection.execute(f"SELECT * FROM hydro_endmembers WHERE active=1 AND id IN ({','.join('?' for _ in ids)}) ORDER BY id",ids).fetchall()
        if len(endmembers)!=len(ids): raise ValueError("endmember_not_found")
        input_data={**payload,"endmember_ids":ids,"sample":dict(sample),"endmembers":[dict(e) for e in endmembers]}
        key=_digest(input_data); now=_now()
        with transaction(immediate=True) as connection:
            old=connection.execute("SELECT * FROM hydro_inversions WHERE task_key=?",(key,)).fetchone()
            if old: return dict(old)
            cursor=connection.execute("INSERT INTO hydro_inversions(sample_id,task_key,model_version,method,input_json,created_at,updated_at) VALUES(?,?,?,?,?,?,?)",(sample_id,key,payload["model_version"],payload["method"],json.dumps(input_data,ensure_ascii=False),now,now))
            return dict(connection.execute("SELECT * FROM hydro_inversions WHERE id=?",(cursor.lastrowid,)).fetchone())

    def run_inversion(self, task_id: int, worker_id: str) -> dict[str, Any]:
        with transaction(immediate=True) as connection:
            task=connection.execute("SELECT * FROM hydro_inversions WHERE id=?",(task_id,)).fetchone()
            if task is None: raise KeyError("task_not_found")
            if task["status"]=="done": return dict(task)
            connection.execute("UPDATE hydro_inversions SET status='running',attempts=attempts+1,worker_id=?,updated_at=? WHERE id=?",(worker_id,_now(),task_id))
        data=json.loads(task["input_json"])
        sample=self.connection.execute("SELECT * FROM hydro_samples WHERE id=?",(task["sample_id"],)).fetchone()
        ids=data["endmember_ids"]
        endmembers=self.connection.execute(f"SELECT * FROM hydro_endmembers WHERE id IN ({','.join('?' for _ in ids)}) ORDER BY id",ids).fetchall()
        try: result=self.solve_mixture(sample,endmembers,data["max_iterations"],data["tolerance"])
        except Exception as exc:
            with transaction(immediate=True) as connection: connection.execute("UPDATE hydro_inversions SET status='failed',error=?,updated_at=? WHERE id=?",(str(exc),_now(),task_id))
            raise
        with transaction(immediate=True) as connection:
            connection.execute("UPDATE hydro_inversions SET status='done',result_json=?,error='',updated_at=? WHERE id=?",(json.dumps(result,ensure_ascii=False),_now(),task_id))
            return dict(connection.execute("SELECT * FROM hydro_inversions WHERE id=?",(task_id,)).fetchone())

    def enqueue_interval(self, inversion_id: int, payload: dict[str, Any]) -> dict[str, Any]:
        inversion=self.connection.execute("SELECT * FROM hydro_inversions WHERE id=?",(inversion_id,)).fetchone()
        if inversion is None: raise KeyError("inversion_not_found")
        if inversion["status"]!="done": raise ValueError("inversion_not_done")
        method=payload["method"]
        seed=payload.get("seed") if method=="parametric-bootstrap" else None
        input_data={"inversion_id":inversion_id,"inversion_task_key":inversion["task_key"],"method":method,"confidence_level":payload["confidence_level"],"replicates":payload["replicates"],"seed":seed,"model_version":payload["model_version"]}
        key=_digest(input_data); now=_now()
        with transaction(immediate=True) as connection:
            old=connection.execute("SELECT * FROM hydro_interval_runs WHERE task_key=?",(key,)).fetchone()
            if old: return dict(old)
            cursor=connection.execute("INSERT INTO hydro_interval_runs(inversion_id,task_key,method,confidence_level,replicates,seed,model_version,input_json,created_at,updated_at) VALUES(?,?,?,?,?,?,?,?,?,?)",(inversion_id,key,method,payload["confidence_level"],payload["replicates"],seed,payload["model_version"],json.dumps(input_data,ensure_ascii=False),now,now))
            return dict(connection.execute("SELECT * FROM hydro_interval_runs WHERE id=?",(cursor.lastrowid,)).fetchone())

    def run_interval(self, task_id: int, worker_id: str, chunk_size: int | None = None) -> dict[str, Any]:
        with transaction(immediate=True) as connection:
            task=connection.execute("SELECT * FROM hydro_interval_runs WHERE id=?",(task_id,)).fetchone()
            if task is None: raise KeyError("task_not_found")
            if task["status"]=="done": return dict(task)
            connection.execute("UPDATE hydro_interval_runs SET status='running',attempts=attempts+1,worker_id=?,error='',updated_at=? WHERE id=?",(worker_id,_now(),task_id))
        try: self._execute_interval(task,chunk_size)
        except Exception as exc:
            with transaction(immediate=True) as connection: connection.execute("UPDATE hydro_interval_runs SET status='failed',error=?,updated_at=? WHERE id=?",(str(exc),_now(),task_id))
            raise
        return dict(self.connection.execute("SELECT * FROM hydro_interval_runs WHERE id=?",(task_id,)).fetchone())

    def _execute_interval(self, task: sqlite3.Row, chunk_size: int | None) -> None:
        task_id=task["id"]
        inversion=self.connection.execute("SELECT * FROM hydro_inversions WHERE id=?",(task["inversion_id"],)).fetchone()
        if inversion is None or inversion["status"]!="done": raise ValueError("inversion_not_done")
        inversion_input=json.loads(inversion["input_json"])
        sample=inversion_input["sample"]; endmembers=inversion_input["endmembers"]
        observed=[sample.get("isotope_d18o"),sample.get("isotope_d2h"),sample.get("solute_mg_l")]
        vectors=[[e["isotope_d18o"],e["isotope_d2h"],e["solute_mg_l"]] for e in endmembers]
        uncertainties=[float(e.get("uncertainty") or 0.0) for e in endmembers]
        measurement_error=float(sample.get("measurement_error") or 0.0)
        max_iterations=int(inversion_input["max_iterations"]); tolerance=float(inversion_input["tolerance"])
        total=task["replicates"] if task["method"]=="parametric-bootstrap" else task["replicates"]*len(endmembers)
        existing={row[0] for row in self.connection.execute("SELECT replicate_index FROM hydro_interval_replicates WHERE run_id=?",(task_id,)).fetchall()}
        pending=[index for index in range(total) if index not in existing]
        if chunk_size is not None: pending=pending[:chunk_size]
        computed=[]
        for index in pending:
            fractions,objective,converged=self._interval_replicate(task,index,observed,vectors,uncertainties,measurement_error,max_iterations,tolerance)
            computed.append((index,json.dumps([round(f,10) for f in fractions]),round(objective,12),1 if converged else 0))
        now=_now()
        with transaction(immediate=True) as connection:
            for index,fractions_json,objective,converged in computed:
                connection.execute("INSERT OR IGNORE INTO hydro_interval_replicates(run_id,replicate_index,fractions_json,objective,converged,created_at) VALUES(?,?,?,?,?,?)",(task_id,index,fractions_json,objective,converged,now))
            stored=connection.execute("SELECT COUNT(*) FROM hydro_interval_replicates WHERE run_id=?",(task_id,)).fetchone()[0]
            if stored<total:
                connection.execute("UPDATE hydro_interval_runs SET status='queued',progress=?,updated_at=? WHERE id=?",(stored,_now(),task_id))
                return
            rows=connection.execute("SELECT * FROM hydro_interval_replicates WHERE run_id=? ORDER BY replicate_index",(task_id,)).fetchall()
            result=self._aggregate_interval(task,rows,endmembers,inversion,observed,vectors,measurement_error)
            connection.execute("UPDATE hydro_interval_runs SET status='done',progress=?,result_json=?,error='',updated_at=? WHERE id=?",(stored,json.dumps(result,ensure_ascii=False),_now(),task_id))

    def _interval_replicate(self, task: sqlite3.Row, index: int, observed: list[float | None], vectors: list[list[float]], uncertainties: list[float], measurement_error: float, max_iterations: int, tolerance: float) -> tuple[list[float], float, bool]:
        if task["method"]=="parametric-bootstrap":
            rng=_seeded_rng(task["seed"],index)
            scale=_tracer_scale(observed)
            active=[i for i,v in enumerate(observed) if v is not None]
            perturbed=list(observed)
            for k in active: perturbed[k]=float(observed[k])+rng.gauss(0.0,measurement_error*scale[k])
            perturbed_vectors=[[value+rng.gauss(0.0,uncertainty*max(1.0,abs(value))) for value in vector] for vector,uncertainty in zip(vectors,uncertainties)]
            core=self._solve_fractions(perturbed,perturbed_vectors,max_iterations,tolerance)
        else:
            replicates=task["replicates"]
            endmember_index=index//replicates
            grid_value=(index%replicates)/(replicates-1)
            core=self._solve_fractions(observed,vectors,max_iterations,tolerance,fixed={endmember_index:grid_value})
        return core["fractions"],core["objective"],core["converged"]

    def _interval_entry(self, endmember: dict[str, Any], point: dict[int, float], lower: float | None, upper: float | None) -> dict[str, Any]:
        return {"endmember_id":endmember["id"],"name":endmember["name"],"point":point.get(endmember["id"]),"lower":round(lower,8) if lower is not None else None,"upper":round(upper,8) if upper is not None else None,"width":round(upper-lower,8) if lower is not None and upper is not None else None}

    def _aggregate_interval(self, task: sqlite3.Row, rows: list[sqlite3.Row], endmembers: list[dict[str, Any]], inversion: sqlite3.Row, observed: list[float | None], vectors: list[list[float]], measurement_error: float) -> dict[str, Any]:
        level=task["confidence_level"]; method=task["method"]; replicates=task["replicates"]
        successes=[(row["replicate_index"],json.loads(row["fractions_json"]),row["objective"]) for row in rows if row["converged"]==1]
        failures=len(rows)-len(successes)
        point={item["endmember_id"]:item["fraction"] for item in json.loads(inversion["result_json"]).get("fractions",[])}
        alpha=(1.0-level)/2.0
        intervals=[]; notes=[]; profile=None
        if method=="parametric-bootstrap":
            for j,endmember in enumerate(endmembers):
                values=sorted(fractions[j] for _,fractions,_ in successes)
                intervals.append(self._interval_entry(endmember,point,_quantile(values,alpha),_quantile(values,1.0-alpha)))
        else:
            active=[i for i,v in enumerate(observed) if v is not None]
            scale=_tracer_scale(observed)
            base_point=[point.get(e["id"],1.0/len(endmembers)) for e in endmembers]
            base_objective=self._objective(base_point,observed,vectors,active,scale)
            chi2=_chi2_df1(level)
            noise=max(measurement_error,1e-6)**2
            threshold=base_objective+chi2*noise
            profile={"base_objective":round(base_objective,12),"chi2_1":round(chi2,8),"noise_variance":noise,"threshold":round(threshold,12)}
            for j,endmember in enumerate(endmembers):
                entries=sorted(((index%replicates)/(replicates-1),objective) for index,_,objective in successes if index//replicates==j)
                accepted=[grid for grid,objective in entries if objective<=threshold]
                if not accepted and entries:
                    accepted=[min(entries,key=lambda item:item[1])[0]]
                    notes.append(f"端元 {endmember['name']} 的剖面网格均未低于阈值，区间退化为网格最优点")
                intervals.append(self._interval_entry(endmember,point,min(accepted) if accepted else None,max(accepted) if accepted else None))
        correlations=[]
        for a in range(len(endmembers)):
            for b in range(a+1,len(endmembers)):
                r=_pearson([fractions[a] for _,fractions,_ in successes],[fractions[b] for _,fractions,_ in successes])
                entry={"endmember_a":endmembers[a]["id"],"name_a":endmembers[a]["name"],"endmember_b":endmembers[b]["id"],"name_b":endmembers[b]["name"],"pearson":round(r,8) if r is not None else None}
                if r is not None and abs(r)>=0.8:
                    entry["warning"]="端元比例强相关，反演结果此消彼长，点估计可能掩盖不确定性"
                    notes.append(f"端元 {endmembers[a]['name']} 与 {endmembers[b]['name']} 相关系数 {round(r,4)}，接近简并")
                correlations.append(entry)
        degenerate=[item["name"] for item in intervals if item["lower"] is not None and item["upper"] is not None and abs(item["upper"]-item["lower"])<=1e-8]
        if degenerate: notes.append("以下端元区间退化为点："+"、".join(degenerate)+"，请检查测量误差与端元不确定性设置")
        return {"inversion_id":task["inversion_id"],"inversion_task_key":inversion["task_key"],"method":method,"confidence_level":level,"model_version":task["model_version"],"seed":task["seed"],"evaluations":len(rows),"samples":len(successes),"failures":failures,"intervals":intervals,"correlations":correlations,"diagnostics":{"correlation_source":"bootstrap-replicates" if method=="parametric-bootstrap" else "profile-ensemble","degenerate_endmembers":degenerate,"notes":notes,"profile":profile}}

    def _interval_view(self, row: sqlite3.Row) -> dict[str, Any]:
        data=dict(row)
        data["result"]=json.loads(data["result_json"] or "{}")
        return data

    def get_interval(self, task_id: int) -> dict[str, Any] | None:
        row=self.connection.execute("SELECT * FROM hydro_interval_runs WHERE id=?",(task_id,)).fetchone()
        return self._interval_view(row) if row else None

    def list_intervals(self, inversion_id: int, confidence_level: float | None = None, model_version: str | None = None) -> list[dict[str, Any]]:
        query="SELECT * FROM hydro_interval_runs WHERE inversion_id=?"
        params: list[Any]=[inversion_id]
        if confidence_level is not None: query+=" AND confidence_level=?"; params.append(confidence_level)
        if model_version is not None: query+=" AND model_version=?"; params.append(model_version)
        query+=" ORDER BY confidence_level,model_version,id"
        return [self._interval_view(row) for row in self.connection.execute(query,params).fetchall()]

    def run_transport(self, well_id: int, payload: dict[str, Any]) -> dict[str, Any]:
        if self.connection.execute("SELECT id FROM hydro_wells WHERE id=?",(well_id,)).fetchone() is None: raise KeyError("well_not_found")
        key=_digest({"well_id":well_id,**payload}); now=_now()
        old=self.connection.execute("SELECT * FROM hydro_transport_runs WHERE task_key=?",(key,)).fetchone()
        if old: return dict(old)
        points=[]; t=payload["step_days"]
        while t<=payload["duration_days"]+1e-12:
            d=payload["dispersion_m2_day"]; x=payload["distance_m"]; v=payload["velocity_m_day"]
            c=payload["source_concentration"]*math.exp(-((x-v*t)**2)/(4*d*t))*math.exp(-payload["decay_per_day"]*t)/math.sqrt(4*math.pi*d*t)
            points.append({"time_days":round(t,8),"concentration":c}); t+=payload["step_days"]
        peak=max(points,key=lambda p:p["concentration"])
        result={"points":points,"peak":peak,"arrival_time_days":payload["distance_m"]/payload["velocity_m_day"],"model_version":payload["model_version"]}
        with transaction(immediate=True) as connection:
            cursor=connection.execute("INSERT INTO hydro_transport_runs(well_id,task_key,model_version,input_json,status,result_json,created_at,updated_at) VALUES(?,?,?,?,?,?,?,?)",(well_id,key,payload["model_version"],json.dumps(payload,ensure_ascii=False),"done",json.dumps(result,ensure_ascii=False),now,now))
            return dict(connection.execute("SELECT * FROM hydro_transport_runs WHERE id=?",(cursor.lastrowid,)).fetchone())
