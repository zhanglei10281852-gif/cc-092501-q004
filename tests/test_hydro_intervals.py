from __future__ import annotations


def _completed_inversion(client, suffix: str, model_version: str = "mix-test") -> dict:
    well = client.post("/api/hydro/wells", json={"code": f"W-{suffix}", "name": "北部监测井", "latitude": 35.1, "longitude": 116.2, "aquifer": "浅层孔隙含水层", "screen_depth_m": 42}).json()
    e1 = client.post("/api/hydro/endmembers", json={"name": f"山区降水-{suffix}", "isotope_d18o": -10, "isotope_d2h": -70, "solute_mg_l": 10, "uncertainty": 0.1, "version": "v1"}).json()
    e2 = client.post("/api/hydro/endmembers", json={"name": f"河流渗漏-{suffix}", "isotope_d18o": -5, "isotope_d2h": -35, "solute_mg_l": 50, "uncertainty": 0.2, "version": "v1"}).json()
    sample = client.post(f"/api/hydro/wells/{well['id']}/samples", json={"sample_code": f"S-{suffix}", "sampled_at": "2026-09-24T08:00:00+00:00", "isotope_d18o": -7.5, "isotope_d2h": -52.5, "solute_mg_l": 30, "detection_limit": 0.1, "measurement_error": 0.05}).json()
    task = client.post(f"/api/hydro/samples/{sample['id']}/inversions", json={"endmember_ids": [e1["id"], e2["id"]], "max_iterations": 1000, "tolerance": 1e-10, "model_version": model_version})
    assert task.status_code == 202, task.text
    done = client.post(f"/api/hydro/inversions/{task.json()['id']}/run?worker_id=test")
    assert done.status_code == 200 and done.json()["status"] == "done", done.text
    return done.json()


def _run_to_done(client, task_id: int, chunk_size: int | None = None) -> dict:
    query = f"/api/hydro/intervals/{task_id}/run?worker_id=test"
    if chunk_size is not None:
        query += f"&chunk_size={chunk_size}"
    for _ in range(50):
        response = client.post(query)
        assert response.status_code == 200, response.text
        if response.json()["status"] == "done":
            return response.json()
    raise AssertionError("区间任务未在预期次数内完成")


def test_bootstrap_interval_is_reproducible_and_resumable(client):
    inversion = _completed_inversion(client, "101")
    point_estimate = inversion["result_json"]
    payload = {"method": "parametric-bootstrap", "confidence_level": 0.95, "replicates": 60, "seed": 42, "model_version": "interval-1"}
    created = client.post(f"/api/hydro/inversions/{inversion['id']}/intervals", json=payload)
    assert created.status_code == 202, created.text
    duplicate = client.post(f"/api/hydro/inversions/{inversion['id']}/intervals", json=payload)
    assert duplicate.json()["id"] == created.json()["id"]
    task_id = created.json()["id"]
    partial = client.post(f"/api/hydro/intervals/{task_id}/run?worker_id=test&chunk_size=25")
    assert partial.json()["status"] == "queued" and partial.json()["progress"] == 25
    done = _run_to_done(client, task_id, chunk_size=25)
    result = client.get(f"/api/hydro/intervals/{task_id}").json()["result"]
    assert result["samples"] + result["failures"] == 60
    assert result["evaluations"] == 60
    assert len(result["intervals"]) == 2
    for entry in result["intervals"]:
        assert 0.0 <= entry["lower"] <= entry["point"] <= entry["upper"] <= 1.0
        assert entry["width"] > 0
    assert len(result["correlations"]) == 1
    correlation = result["correlations"][0]
    assert correlation["pearson"] <= -0.99
    assert "warning" in correlation
    # 重试已完成的任务不改变结果
    again = client.post(f"/api/hydro/intervals/{task_id}/run?worker_id=test")
    assert again.json()["result_json"] == done["result_json"]
    # 区间计算不覆盖原始点估计
    rerun = client.post(f"/api/hydro/inversions/{inversion['id']}/run?worker_id=test")
    assert rerun.json()["result_json"] == point_estimate


def test_chunked_and_single_run_give_identical_results(client):
    first = _completed_inversion(client, "201")
    second = _completed_inversion(client, "202")
    payload = {"method": "parametric-bootstrap", "confidence_level": 0.9, "replicates": 40, "seed": 7, "model_version": "interval-1"}
    task_a = client.post(f"/api/hydro/inversions/{first['id']}/intervals", json=payload).json()
    task_b = client.post(f"/api/hydro/inversions/{second['id']}/intervals", json=payload).json()
    single = _run_to_done(client, task_a["id"])
    chunked = _run_to_done(client, task_b["id"], chunk_size=7)
    import json
    result_single = json.loads(single["result_json"])
    result_chunked = json.loads(chunked["result_json"])
    numeric = lambda entries: [[item[key] for key in ("point", "lower", "upper", "width")] for item in entries]
    assert numeric(result_single["intervals"]) == numeric(result_chunked["intervals"])
    pearson = lambda entries: [item["pearson"] for item in entries]
    assert pearson(result_single["correlations"]) == pearson(result_chunked["correlations"])
    assert (result_single["samples"], result_single["failures"]) == (result_chunked["samples"], result_chunked["failures"])


def test_profile_intervals_and_confidence_level_comparison(client):
    inversion = _completed_inversion(client, "301")
    low = client.post(f"/api/hydro/inversions/{inversion['id']}/intervals", json={"method": "deterministic-profile", "confidence_level": 0.9, "replicates": 41, "model_version": "interval-1"}).json()
    high = client.post(f"/api/hydro/inversions/{inversion['id']}/intervals", json={"method": "deterministic-profile", "confidence_level": 0.99, "replicates": 41, "model_version": "interval-1"}).json()
    assert low["id"] != high["id"]
    done_low = _run_to_done(client, low["id"])
    done_high = _run_to_done(client, high["id"])
    import json
    result_low = json.loads(done_low["result_json"])
    result_high = json.loads(done_high["result_json"])
    assert result_low["samples"] + result_low["failures"] == 82
    for entry_low, entry_high in zip(result_low["intervals"], result_high["intervals"]):
        assert entry_low["lower"] <= entry_low["point"] <= entry_low["upper"]
        assert entry_high["lower"] <= entry_low["lower"]
        assert entry_high["upper"] >= entry_low["upper"]
    # 比较接口：按置信水平和模型版本过滤
    everything = client.get(f"/api/hydro/inversions/{inversion['id']}/intervals").json()["items"]
    assert len(everything) == 2
    only_high = client.get(f"/api/hydro/inversions/{inversion['id']}/intervals?confidence_level=0.99").json()["items"]
    assert len(only_high) == 1 and only_high[0]["result"]["confidence_level"] == 0.99
    other_version = client.get(f"/api/hydro/inversions/{inversion['id']}/intervals?model_version=interval-2").json()["items"]
    assert other_version == []


def test_seed_changes_bootstrap_task_and_result(client):
    inversion = _completed_inversion(client, "401")
    base = {"method": "parametric-bootstrap", "confidence_level": 0.95, "replicates": 30, "model_version": "interval-1"}
    first = client.post(f"/api/hydro/inversions/{inversion['id']}/intervals", json={**base, "seed": 1}).json()
    second = client.post(f"/api/hydro/inversions/{inversion['id']}/intervals", json={**base, "seed": 2}).json()
    assert first["id"] != second["id"]
    import json
    assert json.loads(_run_to_done(client, first["id"])["result_json"])["intervals"] != json.loads(_run_to_done(client, second["id"])["result_json"])["intervals"]


def test_interval_requires_completed_inversion(client):
    well = client.post("/api/hydro/wells", json={"code": "W-501", "name": "北部监测井", "latitude": 35.1, "longitude": 116.2, "aquifer": "浅层孔隙含水层", "screen_depth_m": 42}).json()
    e1 = client.post("/api/hydro/endmembers", json={"name": "山区降水-501", "isotope_d18o": -10, "isotope_d2h": -70, "solute_mg_l": 10, "uncertainty": 0.1, "version": "v1"}).json()
    e2 = client.post("/api/hydro/endmembers", json={"name": "河流渗漏-501", "isotope_d18o": -5, "isotope_d2h": -35, "solute_mg_l": 50, "uncertainty": 0.2, "version": "v1"}).json()
    sample = client.post(f"/api/hydro/wells/{well['id']}/samples", json={"sample_code": "S-501", "sampled_at": "2026-09-24T08:00:00+00:00", "isotope_d18o": -7.5, "isotope_d2h": -52.5, "solute_mg_l": 30, "detection_limit": 0.1, "measurement_error": 0.05}).json()
    task = client.post(f"/api/hydro/samples/{sample['id']}/inversions", json={"endmember_ids": [e1["id"], e2["id"]], "model_version": "mix-test"}).json()
    response = client.post(f"/api/hydro/inversions/{task['id']}/intervals", json={"replicates": 20})
    assert response.status_code == 422
    missing = client.post("/api/hydro/inversions/99999/intervals", json={"replicates": 20})
    assert missing.status_code == 404
