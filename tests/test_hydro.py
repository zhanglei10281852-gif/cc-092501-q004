from __future__ import annotations

import json


def create_well(client,code="W-001"):
    response=client.post("/api/hydro/wells",json={"code":code,"name":"北部监测井","latitude":35.1,"longitude":116.2,"aquifer":"浅层孔隙含水层","screen_depth_m":42})
    assert response.status_code==201,response.text
    return response.json()


def test_mixture_inversion_and_transport(client):
    well=create_well(client)
    e1=client.post("/api/hydro/endmembers",json={"name":"山区降水","isotope_d18o":-10,"isotope_d2h":-70,"solute_mg_l":10,"uncertainty":0.1,"version":"v1"}).json()
    e2=client.post("/api/hydro/endmembers",json={"name":"河流渗漏","isotope_d18o":-5,"isotope_d2h":-35,"solute_mg_l":50,"uncertainty":0.2,"version":"v1"}).json()
    sample=client.post(f"/api/hydro/wells/{well['id']}/samples",json={"sample_code":"S-001","sampled_at":"2026-09-24T08:00:00+00:00","isotope_d18o":-7.5,"isotope_d2h":-52.5,"solute_mg_l":30,"detection_limit":0.1,"measurement_error":0.05}).json()
    task=client.post(f"/api/hydro/samples/{sample['id']}/inversions",json={"endmember_ids":[e1['id'],e2['id']],"max_iterations":1000,"tolerance":1e-10,"model_version":"mix-test"})
    assert task.status_code==202,task.text
    done=client.post(f"/api/hydro/inversions/{task.json()['id']}/run?worker_id=test")
    assert done.status_code==200,done.text
    assert done.json()["status"]=="done"
    transport=client.post(f"/api/hydro/wells/{well['id']}/transport",json={"source_concentration":100,"distance_m":100,"velocity_m_day":2,"dispersion_m2_day":5,"decay_per_day":0.01,"duration_days":100,"step_days":5,"model_version":"ade-test"})
    assert transport.status_code==201,transport.text
    assert transport.json()["result_json"]


def test_missing_measurement_is_classified(client):
    well=create_well(client,"W-002")
    sample=client.post(f"/api/hydro/wells/{well['id']}/samples",json={"sample_code":"S-002","sampled_at":"2026-09-24T08:00:00+00:00","isotope_d18o":-7.5,"detection_limit":0.1,"measurement_error":0.05})
    assert sample.status_code==201
    assert sample.json()["quality_status"]=="incomplete"


def _setup_done_inversion(client, well_code="W-100", sample_code="S-100", close_endmembers=False):
    well=create_well(client,well_code)
    e1=client.post("/api/hydro/endmembers",json={"name":"山区降水","isotope_d18o":-10,"isotope_d2h":-70,"solute_mg_l":10,"uncertainty":0.15,"version":"v1"}).json()
    if close_endmembers:
        e2_payload={"name":"河流渗漏","isotope_d18o":-9.7,"isotope_d2h":-68.5,"solute_mg_l":10.6,"uncertainty":0.15,"version":"v1"}
    else:
        e2_payload={"name":"河流渗漏","isotope_d18o":-5,"isotope_d2h":-35,"solute_mg_l":50,"uncertainty":0.15,"version":"v1"}
    e2=client.post("/api/hydro/endmembers",json=e2_payload).json()
    sample=client.post(f"/api/hydro/wells/{well['id']}/samples",json={"sample_code":sample_code,"sampled_at":"2026-09-24T08:00:00+00:00","isotope_d18o":-7.5,"isotope_d2h":-52.5,"solute_mg_l":30,"detection_limit":0.1,"measurement_error":0.05}).json()
    task=client.post(f"/api/hydro/samples/{sample['id']}/inversions",json={"endmember_ids":[e1["id"],e2["id"]],"max_iterations":1000,"tolerance":1e-10,"model_version":"mix-test"})
    assert task.status_code==202,task.text
    done=client.post(f"/api/hydro/inversions/{task.json()['id']}/run?worker_id=test")
    assert done.status_code==200,done.text
    assert done.json()["status"]=="done"
    return well, e1, e2, sample, done.json()


def _run_interval(client, inversion_id, payload, worker="worker-a"):
    queued=client.post(f"/api/hydro/inversions/{inversion_id}/intervals",json=payload)
    assert queued.status_code==202,queued.text
    ran=client.post(f"/api/hydro/intervals/{queued.json()['id']}/run?worker_id={worker}")
    assert ran.status_code==200,ran.text
    return queued.json(), ran.json()


def test_bootstrap_interval_full_flow_and_shape(client):
    _,_,_,_,inversion=_setup_done_inversion(client)
    queued, ran = _run_interval(client, inversion["id"],
        {"method":"parametric-bootstrap","confidence_level":0.95,"random_seed":20260924,"n_bootstrap":200})
    assert queued["status"]=="queued"
    assert ran["status"]=="done"
    assert ran["attempts"]==1
    result=json.loads(ran["result_json"])
    assert result["method"]=="parametric-bootstrap"
    assert result["confidence_level"]==0.95
    summary=result["summary"]
    assert summary["n_samples"]==200
    assert summary["n_success"]+summary["n_failures"]==200
    assert summary["success_rate"]==1.0
    intervals=result["intervals"]
    assert len(intervals)==2
    for item in intervals:
        assert 0.0<=item["lower"]<=item["fraction_point"]<=item["upper"]<=1.0 or item["contains_point"]
        assert item["lower"]<=item["upper"]
        assert {"endmember_id","name","fraction_point","lower","upper","width"} <= set(item)
    matrix=result["correlations"]["matrix"]
    assert matrix[0][0]==1.0 and matrix[1][1]==1.0
    assert abs(matrix[0][1]-matrix[1][0])<1e-12
    # 点估计快照与原始任务保持不变
    snapshot=json.loads(ran["point_snapshot_json"])
    assert snapshot["result_json"]==json.loads(inversion["result_json"])
    # 结果中同时带回点估计，便于对照
    assert {p["endmember_id"] for p in result["point_estimate"]}=={i["endmember_id"] for i in intervals}


def test_interval_task_is_recoverable_and_deterministic(client):
    _,_,_,_,inversion=_setup_done_inversion(client,"W-101","S-101")
    payload={"method":"parametric-bootstrap","confidence_level":0.9,"random_seed":77,"n_bootstrap":150}
    first=client.post(f"/api/hydro/inversions/{inversion['id']}/intervals",json=payload)
    assert first.status_code==202
    first_task=first.json()
    # 同配置重复入队返回同一任务（幂等）
    again=client.post(f"/api/hydro/inversions/{inversion['id']}/intervals",json=payload)
    assert again.json()["id"]==first_task["id"]
    ran=client.post(f"/api/hydro/intervals/{first_task['id']}/run?worker_id=w1").json()
    result_one=ran["result_json"]
    # 已完成任务重跑不重算：attempts 不增、结果逐位一致
    rerun=client.post(f"/api/hydro/intervals/{first_task['id']}/run?worker_id=w2").json()
    assert rerun["attempts"]==ran["attempts"]==1
    assert rerun["result_json"]==result_one
    # 模拟失败后恢复：将任务置回 queued 再跑，结果仍与首次完全一致
    from app.database import get_connection
    get_connection().execute("UPDATE hydro_inversion_intervals SET status='queued' WHERE id=?",(first_task["id"],))
    recovered=client.post(f"/api/hydro/intervals/{first_task['id']}/run?worker_id=w3").json()
    assert recovered["attempts"]==2
    assert recovered["status"]=="done"
    assert recovered["result_json"]==result_one
    # 原始点估计结果没有被覆盖
    untouched=client.post(f"/api/hydro/inversions/{inversion['id']}/run?worker_id=test").json()
    assert untouched["result_json"]==inversion["result_json"]


def test_profile_interval_is_deterministic_and_levels_compare(client):
    _,_,_,_,inversion=_setup_done_inversion(client,"W-102","S-102")
    _, ran95 = _run_interval(client, inversion["id"],
        {"method":"deterministic-profile","confidence_level":0.95,"grid_points":25})
    r95=json.loads(ran95["result_json"])
    assert r95["method"]=="deterministic-profile"
    assert r95["summary"]["n_samples"]==25*2
    assert r95["summary"]["n_failures"]==0
    # 剖面法不含随机量：重跑逐位一致
    rerun=client.post(f"/api/hydro/intervals/{ran95['id']}/run?worker_id=w").json()
    assert rerun["result_json"]==ran95["result_json"]
    # 90% 区间不宽于 95%
    _, ran90 = _run_interval(client, inversion["id"],
        {"method":"deterministic-profile","confidence_level":0.90,"grid_points":25})
    r90=json.loads(ran90["result_json"])
    for lo, hi in zip(r90["intervals"], r95["intervals"]):
        assert lo["width"]<=hi["width"]+1e-9
    # 阈值随置信水平变化
    assert r90["diagnostics"]["profile_threshold_chi2_1"]<r95["diagnostics"]["profile_threshold_chi2_1"]
    # 区间列表接口支持按样本/版本/水平/方法比较
    listing=client.get(f"/api/hydro/intervals?sample_id={inversion['sample_id']}&model_version=mix-test&confidence_level=0.95&method=deterministic-profile")
    assert listing.status_code==200
    rows=listing.json()
    assert len(rows)==1 and rows[0]["id"]==ran95["id"]
    all_rows=client.get(f"/api/hydro/intervals?sample_id={inversion['sample_id']}").json()
    assert len(all_rows)==2


def test_different_seeds_and_levels_create_distinct_tasks(client):
    _,_,_,_,inversion=_setup_done_inversion(client,"W-103","S-103")
    t1=client.post(f"/api/hydro/inversions/{inversion['id']}/intervals",
                   json={"method":"parametric-bootstrap","confidence_level":0.95,"random_seed":1,"n_bootstrap":100}).json()
    t2=client.post(f"/api/hydro/inversions/{inversion['id']}/intervals",
                   json={"method":"parametric-bootstrap","confidence_level":0.95,"random_seed":2,"n_bootstrap":100}).json()
    t3=client.post(f"/api/hydro/inversions/{inversion['id']}/intervals",
                   json={"method":"parametric-bootstrap","confidence_level":0.90,"random_seed":1,"n_bootstrap":100}).json()
    assert len({t["task_key"] for t in (t1,t2,t3)})==3


def test_interval_requires_completed_point_estimate(client):
    well=create_well(client,"W-104")
    e1=client.post("/api/hydro/endmembers",json={"name":"山区降水","isotope_d18o":-10,"isotope_d2h":-70,"solute_mg_l":10,"uncertainty":0.1,"version":"v1"}).json()
    e2=client.post("/api/hydro/endmembers",json={"name":"河流渗漏","isotope_d18o":-5,"isotope_d2h":-35,"solute_mg_l":50,"uncertainty":0.1,"version":"v1"}).json()
    sample=client.post(f"/api/hydro/wells/{well['id']}/samples",json={"sample_code":"S-104","sampled_at":"2026-09-24T08:00:00+00:00","isotope_d18o":-7.5,"isotope_d2h":-52.5,"solute_mg_l":30,"detection_limit":0.1,"measurement_error":0.05}).json()
    queued=client.post(f"/api/hydro/samples/{sample['id']}/inversions",json={"endmember_ids":[e1["id"],e2["id"]],"model_version":"mix-test"})
    missing=client.post(f"/api/hydro/inversions/9999/intervals",json={"confidence_level":0.95})
    assert missing.status_code==404
    not_done=client.post(f"/api/hydro/inversions/{queued.json()['id']}/intervals",json={"confidence_level":0.95})
    assert not_done.status_code==422
    # 无效置信水平被拒绝
    done=client.post(f"/api/hydro/inversions/{queued.json()['id']}/run?worker_id=test")
    assert done.status_code==200
    bad=client.post(f"/api/hydro/inversions/{queued.json()['id']}/intervals",json={"confidence_level":0.4})
    assert bad.status_code==422


def test_close_endmembers_produce_wide_intervals_and_diagnostics(client):
    _,_,_,_,inversion=_setup_done_inversion(client,"W-105","S-105",close_endmembers=True)
    _, ran = _run_interval(client, inversion["id"],
        {"method":"parametric-bootstrap","confidence_level":0.95,"random_seed":3,"n_bootstrap":300})
    result=json.loads(ran["result_json"])
    warnings=";".join(result["diagnostics"]["warnings"])
    assert "poorly_separated" in warnings
    # 区间比端元可区分时更宽（任一触及单纯形边界或宽度显著）
    assert any(item["width"]>=0.5 or item["boundary"] for item in result["intervals"])
    distance=result["correlations"]["endmember_std_distance"]
    assert distance[0][1]<3.0


def test_model_versions_are_comparable_via_listing(client):
    _,_,_,sample,v1=_setup_done_inversion(client,"W-106","S-106")
    # 同一数据再做一个不同模型版本的点估计
    task2=client.post(f"/api/hydro/samples/{sample['id']}/inversions",
                      json={"endmember_ids":_endmember_ids(v1),"model_version":"mix-test-2","max_iterations":1000,"tolerance":1e-10})
    v2=client.post(f"/api/hydro/inversions/{task2.json()['id']}/run?worker_id=test").json()
    assert v2["model_version"]=="mix-test-2"
    _run_interval(client,v1["id"],{"method":"deterministic-profile","confidence_level":0.95,"grid_points":15})
    _run_interval(client,v2["id"],{"method":"deterministic-profile","confidence_level":0.95,"grid_points":15})
    rows=client.get(f"/api/hydro/intervals?sample_id={sample['id']}").json()
    assert {r["model_version"] for r in rows}=={"mix-test","mix-test-2"}
    only_v2=client.get("/api/hydro/intervals?model_version=mix-test-2").json()
    assert len(only_v2)==1


def _endmember_ids(inversion_row):
    return json.loads(inversion_row["input_json"])["endmember_ids"]
