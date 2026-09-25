from __future__ import annotations

from fastapi import APIRouter, HTTPException, Query

from app.hydro.schemas import EndmemberCreate, IntervalRequest, InversionRequest, SampleCreate, TransportRequest, WellCreate
from app.hydro.service import HydroService

router=APIRouter(prefix="/api/hydro",tags=["地下水科学计算"])

def service()->HydroService: return HydroService()

@router.post("/wells",status_code=201)
def create_well(payload:WellCreate):
    try: return service().create_well(payload.model_dump())
    except Exception as exc:
        if "UNIQUE" in str(exc).upper(): raise HTTPException(409,"井点编码已存在") from exc
        raise

@router.get("/wells/{well_id}")
def get_well(well_id:int):
    value=service().get_well(well_id)
    if value is None: raise HTTPException(404,"井点不存在")
    return value

@router.delete("/wells/{well_id}")
def delete_well(well_id:int):
    try: service().delete_well(well_id); return {"message":"井点已删除"}
    except KeyError as exc: raise HTTPException(404,"井点不存在") from exc

@router.post("/endmembers",status_code=201)
def create_endmember(payload:EndmemberCreate): return service().create_endmember(payload.model_dump())

@router.post("/wells/{well_id}/samples",status_code=201)
def add_sample(well_id:int,payload:SampleCreate):
    try: return service().add_sample(well_id,payload.model_dump())
    except KeyError as exc: raise HTTPException(404,"井点不存在") from exc

@router.post("/samples/{sample_id}/inversions",status_code=202)
def enqueue_inversion(sample_id:int,payload:InversionRequest):
    try: return service().enqueue_inversion(sample_id,payload.model_dump())
    except KeyError as exc: raise HTTPException(404,"样本不存在") from exc
    except ValueError as exc: raise HTTPException(422,str(exc)) from exc

@router.post("/inversions/{task_id}/run")
def run_inversion(task_id:int,worker_id:str=Query(...,min_length=1)):
    try: return service().run_inversion(task_id,worker_id)
    except KeyError as exc: raise HTTPException(404,"任务不存在") from exc
    except ValueError as exc: raise HTTPException(422,str(exc)) from exc

@router.post("/inversions/{inversion_id}/intervals",status_code=202)
def enqueue_interval(inversion_id:int,payload:IntervalRequest):
    try: return service().enqueue_interval(inversion_id,payload.model_dump())
    except KeyError as exc: raise HTTPException(404,"反演任务不存在") from exc
    except ValueError as exc: raise HTTPException(422,str(exc)) from exc

@router.post("/intervals/{task_id}/run")
def run_interval(task_id:int,worker_id:str=Query(...,min_length=1),chunk_size:int|None=Query(default=None,ge=1,le=5000)):
    try: return service().run_interval(task_id,worker_id,chunk_size)
    except KeyError as exc: raise HTTPException(404,"区间任务不存在") from exc
    except ValueError as exc: raise HTTPException(422,str(exc)) from exc

@router.get("/intervals/{task_id}")
def get_interval(task_id:int):
    value=service().get_interval(task_id)
    if value is None: raise HTTPException(404,"区间任务不存在")
    return value

@router.get("/inversions/{inversion_id}/intervals")
def list_intervals(inversion_id:int,confidence_level:float|None=Query(default=None),model_version:str|None=Query(default=None)):
    return {"items":service().list_intervals(inversion_id,confidence_level,model_version)}

@router.post("/wells/{well_id}/transport",status_code=201)
def run_transport(well_id:int,payload:TransportRequest):
    try: return service().run_transport(well_id,payload.model_dump())
    except KeyError as exc: raise HTTPException(404,"井点不存在") from exc
