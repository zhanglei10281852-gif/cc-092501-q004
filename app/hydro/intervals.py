"""混合反演端元比例的区间估计纯计算模块。

本模块不访问数据库、不读取时钟，所有结果只由输入决定：

* ``bootstrap_intervals`` —— 固定随机种子的参数自助法（parametric bootstrap），
  给出经验分位数置信区间与端元比例的 Pearson 相关矩阵；
* ``profile_intervals`` —— 确定性剖面似然法（profile likelihood），
  给出单参数剖面置信区间，并在最优点用数值 Hessian 给出局部相关性诊断。

随机数只使用 ``random.Random`` 且以整数播种（Python 文档保证 MT 序列跨版本
稳定），子样本流由整数运算派生，避免依赖受 ``PYTHONHASHSEED`` 影响的对象哈希。
"""

from __future__ import annotations

import math
import random
from typing import Any, Mapping, Sequence

TRACERS = ("isotope_d18o", "isotope_d2h", "solute_mg_l")
BOOTSTRAP = "parametric-bootstrap"
PROFILE = "deterministic-profile"


# ---------------------------------------------------------------------------
# 基础数值工具
# ---------------------------------------------------------------------------

def _observed(sample: Mapping[str, Any]) -> list[float | None]:
    return [sample["isotope_d18o"], sample["isotope_d2h"], sample["solute_mg_l"]]


def _vectors(endmembers: Sequence[Mapping[str, Any]]) -> list[list[float]]:
    return [[float(e[key]) for key in TRACERS] for e in endmembers]


def _scales(observed: Sequence[float | None]) -> list[float]:
    return [20.0, 100.0, max(1.0, float(observed[2] or 1.0))]


def _active(observed: Sequence[float | None]) -> list[int]:
    return [k for k, value in enumerate(observed) if value is not None]


def _predict(fractions: Sequence[float], vectors: Sequence[Sequence[float]]) -> list[float]:
    return [sum(f * v[k] for f, v in zip(fractions, vectors)) for k in range(3)]


def _rss(fractions: Sequence[float], observed: Sequence[float], active: Sequence[int],
         vectors: Sequence[Sequence[float]], scales: Sequence[float]) -> float:
    predicted = _predict(fractions, vectors)
    return sum(((predicted[k] - observed[k]) / scales[k]) ** 2 for k in active)


def _project_sum(values: Sequence[float], total: float) -> list[float]:
    """投影到 ``x>=0 且 sum(x)=total`` 的单纯形（截距平移投影）。"""
    if total <= 1e-15:
        return [0.0] * len(values)
    clipped = [max(0.0, v) for v in values]
    if sum(clipped) <= 1e-15:
        return [total / len(values)] * len(values)
    # 二分求截断水线，使裁剪后总和恰好为 total
    low, high = 0.0, max(clipped)
    for _ in range(60):
        mid = (low + high) / 2
        if sum(max(0.0, v - mid) for v in clipped) <= total:
            high = mid
        else:
            low = mid
    scale = sum(max(0.0, v - high) for v in clipped)
    if scale <= 1e-18:
        return [total / len(values)] * len(values)
    return [total * max(0.0, v - high) / scale for v in clipped]


def solve_fractions(observed: Sequence[float | None], vectors: Sequence[Sequence[float]],
                    active: Sequence[int], scales: Sequence[float],
                    max_iterations: int, tolerance: float,
                    fixed: Mapping[int, float] | None = None) -> dict[str, Any]:
    """投影梯度求解混合比例。

    ``fixed`` 给定时把指定端元比例锁在给定值，其余端元分配剩余质量
    （质量守恒由构造满足）；不给定时与点估计 ``solve_mixture`` 同算法，
    含质量守恒罚项，保证自助分布围绕原始点估计。
    """
    n = len(vectors)
    rate = 0.08
    if not fixed:
        fractions = [1.0 / n] * n
        last = float("inf")
        objective = float("inf")
        iteration = 0
        for iteration in range(max_iterations):
            predicted = _predict(fractions, vectors)
            residual = [(predicted[k] - float(observed[k])) / scales[k] if k in active else 0.0
                        for k in range(3)]
            objective = sum(r * r for r in residual) + ((sum(fractions) - 1.0) * 10) ** 2
            if abs(last - objective) < tolerance:
                break
            last = objective
            gradient = [2 * sum(residual[k] * vectors[j][k] / scales[k] for k in active)
                        for j in range(n)]
            fractions = _project_sum([f - rate * g for f, g in zip(fractions, gradient)], 1.0)
        return {"fractions": fractions, "converged": abs(last - objective) < tolerance,
                "iterations": iteration + 1, "rss": _rss(fractions, observed, active, vectors, scales)}

    fixed = dict(fixed)
    fixed_mass = sum(fixed.values())
    if fixed_mass > 1.0 + 1e-12 or any(not (0.0 <= p <= 1.0 + 1e-12) for p in fixed.values()):
        raise ValueError("fixed_fraction_infeasible")
    free_idx = [i for i in range(n) if i not in fixed]
    free_mass = max(0.0, 1.0 - fixed_mass)
    fractions = [0.0] * n
    for j, p in fixed.items():
        fractions[j] = min(1.0, max(0.0, p))
    current = [free_mass / len(free_idx)] * len(free_idx) if free_idx else []
    if free_mass <= 1e-15 or not free_idx:
        for slot, j in enumerate(free_idx):
            fractions[j] = 0.0
        return {"fractions": fractions, "converged": True, "iterations": 0,
                "rss": _rss(fractions, observed, active, vectors, scales)}

    def assemble(parts: Sequence[float]) -> list[float]:
        out = [0.0] * n
        for j, p in fixed.items():
            out[j] = p
        for slot, j in enumerate(free_idx):
            out[j] = parts[slot]
        return out

    last = float("inf")
    objective = float("inf")
    iteration = 0
    for iteration in range(max_iterations):
        fractions = assemble(current)
        predicted = _predict(fractions, vectors)
        residual = [(predicted[k] - float(observed[k])) / scales[k] if k in active else 0.0
                    for k in range(3)]
        objective = sum(r * r for r in residual)
        if abs(last - objective) < tolerance:
            break
        last = objective
        gradient = [2 * sum(residual[k] * vectors[j][k] / scales[k] for k in active)
                    for j in free_idx]
        current = _project_sum([u - rate * g for u, g in zip(current, gradient)], free_mass)
    fractions = assemble(current)
    return {"fractions": fractions, "converged": abs(last - objective) < tolerance,
            "iterations": iteration + 1,
            "rss": _rss(fractions, observed, active, vectors, scales)}


# ---------------------------------------------------------------------------
# 分位数与相关性
# ---------------------------------------------------------------------------

def quantile_type7(sorted_values: Sequence[float], probability: float) -> float:
    """线性插值经验分位数（与 numpy 默认 quantile 的线性插值定义一致）。"""
    n = len(sorted_values)
    if n == 0:
        raise ValueError("empty_sample")
    if n == 1:
        return sorted_values[0]
    position = (n - 1) * probability
    low = int(math.floor(position))
    fraction = position - low
    upper = min(low + 1, n - 1)
    return sorted_values[low] * (1.0 - fraction) + sorted_values[upper] * fraction


def pearson_matrix(columns: Sequence[Sequence[float]]) -> tuple[list[list[float]], list[int]]:
    """返回 Pearson 相关矩阵与常数列下标（常数列无法定义相关，记为 0）。"""
    n_cols = len(columns)
    means = [sum(col) / len(col) for col in columns]
    centered = [[v - means[j] for v in columns[j]] for j in range(n_cols)]
    stds = [math.sqrt(sum(v * v for v in col) / len(col)) for col in centered]
    constant = [j for j, std in enumerate(stds) if std <= 1e-14]
    matrix = [[0.0] * n_cols for _ in range(n_cols)]
    for j in range(n_cols):
        matrix[j][j] = 1.0
    for i in range(n_cols):
        for j in range(i + 1, n_cols):
            if i in constant or j in constant:
                value = 0.0
            else:
                cov = sum(a * b for a, b in zip(centered[i], centered[j])) / len(columns[i])
                value = cov / (stds[i] * stds[j])
            matrix[i][j] = matrix[j][i] = max(-1.0, min(1.0, value))
    return matrix, constant


def normal_ppf(probability: float) -> float:
    """标准正态分位数（Acklam 有理逼近，|误差| < ~1.2e-7）。

    剖面法只需要 χ²(1) 分位数，而 χ²(1,p)=z((1+p)/2)²，因此无需通用 χ² 库。
    """
    p = min(1.0 - 1e-12, max(1e-12, probability))
    a = (-3.969683028665376e01, 2.209460984245205e02, -2.759285104469687e02,
         1.383577518672690e02, -3.066479806614716e01, 2.506628277459239e00)
    b = (-5.447609879822406e01, 1.615858368580409e02, -1.556989798598866e02,
         6.680131188771972e01, -1.328068155288572e01)
    c = (-7.784894002430293e-03, -3.223964580411365e-01, -2.400758277161838e00,
         -2.549732539343734e00, 4.374664141464968e00, 2.938163982698783e00)
    d = (7.784695709041462e-03, 3.224671290700398e-01, 2.445134137142996e00,
         3.754408661907416e00)
    p_low, p_high = 0.02425, 0.97575
    if p < p_low:
        q = math.sqrt(-2.0 * math.log(p))
        return (((((c[0] * q + c[1]) * q + c[2]) * q + c[3]) * q + c[4]) * q + c[5]) \
            / ((((d[0] * q + d[1]) * q + d[2]) * q + d[3]) * q + 1.0)
    if p <= p_high:
        q = p - 0.5
        r = q * q
        return (((((a[0] * r + a[1]) * r + a[2]) * r + a[3]) * r + a[4]) * r + a[5]) * q \
            / (((((b[0] * r + b[1]) * r + b[2]) * r + b[3]) * r + b[4]) * r + 1.0)
    q = math.sqrt(-2.0 * math.log(1.0 - p))
    return -(((((c[0] * q + c[1]) * q + c[2]) * q + c[3]) * q + c[4]) * q + c[5]) \
        / ((((d[0] * q + d[1]) * q + d[2]) * q + d[3]) * q + 1.0)


# ---------------------------------------------------------------------------
# 端元标签与扰动模型
# ---------------------------------------------------------------------------

def _labels(endmembers: Sequence[Mapping[str, Any]]) -> list[dict[str, Any]]:
    return [{"endmember_id": e["id"], "name": e["name"]} for e in endmembers]


def _point_map(point_result: Mapping[str, Any]) -> dict[int, float]:
    return {int(item["endmember_id"]): float(item["fraction"])
            for item in point_result["fractions"]}


def _perturb_observed(observed: Sequence[float | None], active: Sequence[int],
                      scales: Sequence[float], sample: Mapping[str, Any],
                      rng: random.Random) -> list[float | None]:
    rel_error = float(sample.get("measurement_error") or 0.0)
    floor = float(sample.get("detection_limit") or 0.0)
    result: list[float | None] = [None] * 3
    for k in active:
        sigma = max(abs(float(observed[k])) * rel_error, floor)
        draw = float(observed[k]) + rng.gauss(0.0, sigma)
        result[k] = max(0.0, draw) if k == 2 else draw  # 溶质浓度非负
    return result


def _perturb_endmembers(vectors: Sequence[Sequence[float]], scales: Sequence[float],
                        endmembers: Sequence[Mapping[str, Any]],
                        rng: random.Random) -> list[list[float]]:
    drawn = []
    for e, vector in zip(endmembers, vectors):
        uncertainty = float(e.get("uncertainty") or 0.0)
        row = []
        for k, value in enumerate(vector):
            sigma = max(abs(value), scales[k] * 0.02) * uncertainty
            draw = value + rng.gauss(0.0, sigma)
            row.append(max(0.0, draw) if k == 2 else draw)
        drawn.append(row)
    return drawn


# ---------------------------------------------------------------------------
# 参数自助法
# ---------------------------------------------------------------------------

def bootstrap_intervals(*, sample: Mapping[str, Any],
                        endmembers: Sequence[Mapping[str, Any]],
                        point_result: Mapping[str, Any],
                        confidence_level: float, random_seed: int,
                        n_bootstrap: int, max_iterations: int,
                        tolerance: float) -> dict[str, Any]:
    observed = _observed(sample)
    active = _active(observed)
    vectors = _vectors(endmembers)
    scales = _scales(observed)
    n_members = len(endmembers)
    seed = int(random_seed) % (2 ** 63)

    draws: list[list[float]] = []
    failures = 0
    for replicate in range(n_bootstrap):
        # 整数派生独立子流：同一 (seed, 编号) 永远产生同一次扰动
        rng = random.Random((seed * 1_000_003 + replicate) % (2 ** 63))
        syn_observed = _perturb_observed(observed, active, scales, sample, rng)
        syn_vectors = _perturb_endmembers(vectors, scales, endmembers, rng)
        try:
            solution = solve_fractions(syn_observed, syn_vectors, active, scales,
                                       max_iterations, tolerance)
        except Exception:
            failures += 1
            continue
        if not solution["converged"]:
            failures += 1
            continue
        draws.append(solution["fractions"])

    n_success = len(draws)
    warnings: list[str] = []
    if n_bootstrap and failures / n_bootstrap > 0.1:
        warnings.append(f"bootstrap_failure_rate_high:{failures}/{n_bootstrap}")
    if n_success < 30:
        warnings.append(f"low_effective_sample:{n_success}")
    if n_success * (1.0 - confidence_level) < 10:
        warnings.append("low_tail_resolution_for_confidence_level")
    if n_success == 0:
        raise RuntimeError("bootstrap_all_replicates_failed")

    point = _point_map(point_result)
    labels = _labels(endmembers)
    tail = (1.0 - confidence_level) / 2.0
    intervals = []
    for j, label in enumerate(labels):
        column = sorted(row[j] for row in draws)
        lower = quantile_type7(column, tail)
        upper = quantile_type7(column, 1.0 - tail)
        point_value = point.get(int(label["endmember_id"]), sum(row[j] for row in draws) / n_success)
        boundary = []
        if lower <= 1e-6:
            boundary.append("lower=0")
        if upper >= 1.0 - 1e-6:
            boundary.append("upper=1")
        if boundary:
            warnings.append(f"interval_touches_boundary:endmember_id={label['endmember_id']}")
        intervals.append({
            "endmember_id": label["endmember_id"],
            "name": label["name"],
            "fraction_point": round(point_value, 8),
            "lower": round(max(0.0, lower), 8),
            "upper": round(min(1.0, upper), 8),
            "width": round(upper - lower, 8),
            "contains_point": bool(lower - 1e-9 <= point_value <= upper + 1e-9),
            "boundary": boundary,
        })

    columns = [[row[j] for row in draws] for j in range(n_members)]
    matrix, constant = pearson_matrix(columns)
    for j in constant:
        warnings.append(f"constant_fraction:endmember_id={labels[j]['endmember_id']}")
    sigmas = _known_sigmas(observed, active, vectors, scales,
                           sample, endmembers,
                           [sum(col) / len(col) for col in columns])
    separation, separation_warnings = _separation_matrix(vectors, active, sigmas, labels)
    warnings.extend(separation_warnings)
    correlations = {
        "method": "pearson",
        "labels": labels,
        "matrix": [[round(value, 8) for value in row] for row in matrix],
        "endmember_std_distance": separation,
    }

    return {
        "method": BOOTSTRAP,
        "confidence_level": confidence_level,
        "random_seed": int(random_seed),
        "summary": {
            "n_requested": n_bootstrap,
            "n_samples": n_bootstrap,
            "n_success": n_success,
            "n_failures": failures,
            "success_rate": round(n_success / n_bootstrap, 8),
        },
        "intervals": intervals,
        "correlations": correlations,
        "diagnostics": {
            "warnings": warnings,
            "point_estimate_covered": all(item["contains_point"] for item in intervals),
            "quantile_method": "type7_linear",
            "perturbation": "observed_normal_measurement_error;endmember_normal_uncertainty",
        },
    }


# ---------------------------------------------------------------------------
# 确定性剖面法
# ---------------------------------------------------------------------------

def _profile_curve(j: int, grid_points: int, observed, active, vectors, scales,
                   max_iterations, tolerance, rss0, threshold, hat_value):
    """粗网格扫描 + 穿越点二分细化，返回 (曲线, 可行区间, 失败点数, 评估计数)。

    scales 即已知噪声标准差，故标准化残差平方和本身就是 χ² 尺度，
    delta = RSS(f_j) - RSS0 直接与 χ²(1) 阈值比较。
    """
    coarse = sorted({0.0, 1.0,
                     *(step / (grid_points - 1) for step in range(grid_points)),
                     min(1.0, max(0.0, hat_value))})
    curve = []
    deltas = []
    failures = 0
    evaluations = 0
    for value in coarse:
        solution = solve_fractions(observed, vectors, active, scales,
                                   max_iterations, tolerance, fixed={j: value})
        evaluations += 1
        if not solution["converged"]:
            failures += 1
            curve.append({"fraction": round(value, 8), "delta": None})
            deltas.append(None)
            continue
        delta = solution["rss"] - rss0
        curve.append({"fraction": round(value, 8), "delta": round(delta, 8)})
        deltas.append(delta)

    def delta_at(value: float) -> float:
        solution = solve_fractions(observed, vectors, active, scales,
                                   max_iterations, tolerance, fixed={j: value})
        return solution["rss"] - rss0

    def feasible(delta: float | None) -> bool:
        return delta is not None and delta <= threshold

    pairs = sorted(zip(coarse, deltas), key=lambda pair: pair[0])
    feasible_values = [v for v, d in pairs if feasible(d)]
    roots: list[float] = []
    for index in range(len(pairs) - 1):
        x0, left = pairs[index]
        x1, right = pairs[index + 1]
        if left is None or right is None or (left <= threshold) == (right <= threshold):
            continue
        lo, hi = (x0, x1) if left <= threshold else (x1, x0)
        for _ in range(14):
            mid = (lo + hi) / 2.0
            if delta_at(mid) <= threshold:
                lo = mid
            else:
                hi = mid
        evaluations += 14
        roots.append((lo + hi) / 2.0)

    if feasible_values or roots:
        boundary_points = feasible_values + roots
        lower = min(boundary_points)
        upper = max(boundary_points)
    else:
        # 粗网格上没有任何可行点：从点估计向外扩张，直到越过阈值或单纯形边界
        step = 1.0 / grid_points
        hat = min(1.0, max(0.0, hat_value))
        try:
            hat_delta = delta_at(hat)
            evaluations += 1
        except Exception:
            hat_delta = float("inf")
        if not feasible(hat_delta):
            lower = upper = 0.0
        else:
            lo = hat
            moved = step
            while lo - moved > 0.0 and delta_at(lo - moved) <= threshold:
                lo -= moved
                moved = min(2 * moved, 1.0)
                evaluations += 1
            lo_edge = max(0.0, lo - moved)
            if delta_at(lo_edge) <= threshold:
                lower = 0.0
            else:
                a, b = lo_edge, lo
                for _ in range(14):
                    mid = (a + b) / 2.0
                    if delta_at(mid) <= threshold:
                        b = mid
                    else:
                        a = mid
                evaluations += 14
                lower = (a + b) / 2.0
            hi = hat
            moved = step
            while hi + moved < 1.0 and delta_at(hi + moved) <= threshold:
                hi += moved
                moved = min(2 * moved, 1.0)
                evaluations += 1
            hi_edge = min(1.0, hi + moved)
            if delta_at(hi_edge) <= threshold:
                upper = 1.0
            else:
                a, b = hi, hi_edge
                for _ in range(14):
                    mid = (a + b) / 2.0
                    if delta_at(mid) <= threshold:
                        a = mid
                    else:
                        b = mid
                evaluations += 14
                upper = (a + b) / 2.0
    return curve, max(0.0, min(1.0, lower)), max(0.0, min(1.0, upper)), failures, evaluations


def _local_hessian_correlations(f_hat, observed, active, vectors, scales):
    """在最优点用数值 Hessian 估计端元相关（局部二次近似）。

    自由参数为前 n-1 个端元比例，最后一个由质量守恒确定。
    边界最优点先向内推 1e-3；Hessian 奇异（端元不可区分）时返回 None。
    """
    n = len(vectors)
    if n < 2:
        return None, "single_endmember"
    eps = 1e-3
    q = [f_hat[j] * (1.0 - n * eps) + eps for j in range(n - 1)]
    h = 5e-4

    def rss_q(free: Sequence[float]) -> float | None:
        if any(v < -1e-12 for v in free) or sum(free) > 1.0 + 1e-12:
            return None
        fractions = list(free) + [max(0.0, 1.0 - sum(free))]
        return _rss(fractions, observed, active, vectors, scales)

    size = n - 1
    hessian = [[0.0] * size for _ in range(size)]
    for i in range(size):
        for j in range(i, size):
            points = [0.0, 0.0, 0.0, 0.0]
            coords = [(1, 1), (1, -1), (-1, 1), (-1, -1)]
            ok = True
            for index, (si, sj) in enumerate(coords):
                free = list(q)
                free[i] += si * h
                free[j] += sj * h
                value = rss_q(free)
                if value is None:
                    ok = False
                    break
                points[index] = value
            if not ok:
                return None, "hessian_crosses_boundary"
            value = (points[0] - points[1] - points[2] + points[3]) / (4.0 * h * h)
            hessian[i][j] = hessian[j][i] = value

    # sigma2 因子对相关系数无影响（标量在归一化中消去），直接求逆 Hessian
    try:
        inv = _invert(hessian)
    except ValueError:
        return None, "hessian_singular_rank_deficient"
    # 质量守恒加边：f_last = 1 - sum(q)
    matrix = [[0.0] * n for _ in range(n)]
    for i in range(n - 1):
        for j in range(n - 1):
            matrix[i][j] = inv[i][j]
    for i in range(n - 1):
        matrix[i][n - 1] = matrix[n - 1][i] = -sum(inv[i])
    matrix[n - 1][n - 1] = sum(sum(row) for row in inv)

    stds = [math.sqrt(max(0.0, matrix[i][i])) for i in range(n)]
    if any(std <= 1e-12 for std in stds):
        return None, "hessian_singular_rank_deficient"
    corr = [[matrix[i][j] / (stds[i] * stds[j]) for j in range(n)] for i in range(n)]
    return [[max(-1.0, min(1.0, value)) for value in row] for row in corr], "local_hessian"


def _invert(matrix: Sequence[Sequence[float]]) -> list[list[float]]:
    """带部分主元的高斯-约当求逆；奇异或非正定时报错。"""
    n = len(matrix)
    work = [list(row) + [1.0 if i == j else 0.0 for j in range(n)]
            for i, row in enumerate(matrix)]
    scale = max((abs(matrix[i][j]) for i in range(n) for j in range(n)), default=0.0)
    pivot_floor = max(1e-14, scale * 1e-12)
    for column in range(n):
        pivot = max(range(column, n), key=lambda r: abs(work[r][column]))
        if abs(work[pivot][column]) <= pivot_floor:
            raise ValueError("singular_matrix")
        work[column], work[pivot] = work[pivot], work[column]
        divisor = work[column][column]
        work[column] = [v / divisor for v in work[column]]
        for row in range(n):
            if row == column:
                continue
            factor = work[row][column]
            if factor:
                work[row] = [a - factor * b for a, b in zip(work[row], work[column])]
    return [row[n:] for row in work]


def _separation_matrix(vectors: Sequence[Sequence[float]], active: Sequence[int],
                       sigmas: Sequence[float],
                       labels: Sequence[Mapping[str, Any]]) -> tuple[list[list[float]], list[str]]:
    """端元在观测噪声尺度下的两两标准化距离。

    相关矩阵受质量守恒约束（双端元恒为 -1）不能单独说明端元相似；
    距离 < 噪声尺度（阈值 3σ）才意味着数据无法区分这两个补给来源。
    """
    n = len(vectors)
    matrix = [[0.0] * n for _ in range(n)]
    warnings: list[str] = []
    for i in range(n):
        for j in range(i + 1, n):
            distance = math.sqrt(
                sum(((vectors[i][k] - vectors[j][k]) / sigmas[k]) ** 2 for k in active))
            matrix[i][j] = matrix[j][i] = round(distance, 8)
            if distance < 3.0:
                warnings.append(
                    f"endmembers_poorly_separated:{labels[i]['endmember_id']}~{labels[j]['endmember_id']}"
                    f":std_distance={round(distance, 3)}")
    return matrix, warnings


def _known_sigmas(observed, active, vectors, scales, sample, endmembers, weights):
    """已知噪声模型下每个示踪剂的合成标准差（自助扰动与剖面共用）。

    σ_k² = 观测分析误差² + Σ_j f_j·端元特征误差²（端元项以比例加权的插件近似）。
    返回 3 元数组（非活动示踪剂给占位 1.0）。
    """
    rel_error = float(sample.get("measurement_error") or 0.0)
    floor = float(sample.get("detection_limit") or 0.0)
    sigmas = [1.0, 1.0, 1.0]
    for k in active:
        obs_variance = max(abs(float(observed[k])) * rel_error, floor) ** 2
        endmember_variance = 0.0
        for j, e in enumerate(endmembers):
            sigma_jk = max(abs(vectors[j][k]), scales[k] * 0.02) * float(e.get("uncertainty") or 0.0)
            endmember_variance += weights[j] * sigma_jk ** 2
        total = math.sqrt(obs_variance + endmember_variance)
        sigmas[k] = total if total > 1e-12 else scales[k] * 1e-6
    return sigmas


def profile_intervals(*, sample: Mapping[str, Any],
                      endmembers: Sequence[Mapping[str, Any]],
                      point_result: Mapping[str, Any],
                      confidence_level: float, grid_points: int,
                      max_iterations: int, tolerance: float) -> dict[str, Any]:
    observed = _observed(sample)
    active = _active(observed)
    raw_vectors = _vectors(endmembers)
    base_scales = _scales(observed)
    n_members = len(endmembers)

    point = _point_map(point_result)
    labels = _labels(endmembers)
    f_point = [point.get(int(label["endmember_id"]), 1.0 / n_members)
               for label in labels]
    # 已知方差模型：σ 由分析误差与端元不确定性确定，无需用残差估计噪声尺度
    weight_scales = _known_sigmas(observed, active, raw_vectors, base_scales,
                                  sample, endmembers, f_point)
    vectors = raw_vectors

    # 无约束最优仅用于标定基准与 Hessian 中心；fraction_point 仍报告原始点估计
    optimum = solve_fractions(observed, vectors, active, weight_scales,
                              max_iterations, tolerance)
    f_hat = optimum["fractions"]
    rss0 = min(optimum["rss"], _rss(f_point, observed, active, vectors, weight_scales))

    warnings: list[str] = []
    degrees_of_freedom = len(active) - (n_members - 1)
    if degrees_of_freedom < 1:
        warnings.append(
            f"underidentified_model:tracers={len(active)} endmembers={n_members}")
    # χ²(1, confidence_level) = z((1+c)/2)²，单参数剖面的精确阈值
    z_value = normal_ppf(0.5 + confidence_level / 2.0)
    threshold = z_value * z_value

    intervals = []
    profiles = []
    total_failures = 0
    total_evaluations = 0
    for j, label in enumerate(labels):
        curve, lower, upper, failures, evaluations = _profile_curve(
            j, grid_points, observed, active, vectors, weight_scales,
            max_iterations, tolerance, rss0, threshold, f_hat[j])
        total_failures += failures
        total_evaluations += evaluations
        point_value = f_point[j]
        profiles.append({"endmember_id": label["endmember_id"], "name": label["name"],
                         "points": curve})
        if lower <= 1e-6 or upper >= 1.0 - 1e-6:
            warnings.append(f"interval_touches_boundary:endmember_id={label['endmember_id']}")
        intervals.append({
            "endmember_id": label["endmember_id"],
            "name": label["name"],
            "fraction_point": round(point_value, 8),
            "lower": round(lower, 8),
            "upper": round(upper, 8),
            "width": round(upper - lower, 8),
            "contains_point": bool(lower - 1e-6 <= point_value <= upper + 1e-6),
            "boundary": [b for b, hit in (("lower=0", lower <= 1e-6),
                                          ("upper=1", upper >= 1.0 - 1e-6)) if hit],
        })

    if total_failures:
        warnings.append(f"profile_grid_nonconverged:{total_failures}")

    matrix, corr_method = _local_hessian_correlations(
        f_hat, observed, active, vectors, weight_scales)
    separation, separation_warnings = _separation_matrix(
        vectors, active, weight_scales, labels)
    warnings.extend(separation_warnings)
    correlations: dict[str, Any]
    if matrix is None:
        correlations = {"method": "none", "reason": corr_method, "labels": labels,
                        "matrix": None, "endmember_std_distance": separation}
        warnings.append(f"correlation_unavailable:{corr_method}")
    else:
        correlations = {"method": corr_method, "labels": labels,
                        "matrix": [[round(value, 8) for value in row] for row in matrix],
                        "endmember_std_distance": separation}

    return {
        "method": PROFILE,
        "confidence_level": confidence_level,
        "summary": {
            "n_requested": grid_points * n_members,
            "n_samples": grid_points * n_members,
            "n_success": grid_points * n_members - total_failures,
            "n_failures": total_failures,
            "success_rate": round(
                (grid_points * n_members - total_failures) / (grid_points * n_members), 8),
        },
        "intervals": intervals,
        "correlations": correlations,
        "profiles": profiles,
        "diagnostics": {
            "warnings": warnings,
            "point_estimate_covered": all(item["contains_point"] for item in intervals),
            "profile_threshold_chi2_1": round(threshold, 8),
            "noise_model": "known_sigma:observed_analysis_error+mixture_weighted_endmember_uncertainty",
            "residual_chi2": round(rss0, 8),
            "degrees_of_freedom": degrees_of_freedom,
            "grid_points": grid_points,
            "refine_bisections": 14,
            "total_profile_evaluations": total_evaluations,
        },
    }
