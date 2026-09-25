# 地下水同位素与污染迁移计算服务

这是一个面向水文地质研究团队和环境监管人员的模块化后端，集中管理地下水井、同位素观测、补给端元、混合源反演、污染物迁移、计算任务、参数版本、结果置信区间、用户权限、会话和审计。项目使用 FastAPI 与 SQLite，所有运行数据保存在单个本地数据库文件中，不依赖另行部署的数据库、缓存或消息队列。

## 主要模块

- 井点与样本：登记井点坐标、含水层、采样批次和实验室测量结果。
- 同位素计算：处理稳定同位素、溶质浓度、检测限和质量守恒约束，反演多个补给端元比例。
- 区间估计：对已完成的点估计追加固定随机种子的参数自助法或确定性剖面法，输出每个端元的置信区间、样本数、失败次数与端元相关性/可区分性诊断；区间任务可恢复、可复现，且不覆盖点估计。
- 污染迁移：计算一维平流、弥散和一阶衰减，提供到达时间和浓度曲线。
- 任务与审计：保存参数版本、计算输入摘要、置信区间、失败重试和结果差异。
- 身份与权限：用户、角色、细粒度权限、会话令牌、账号停用和会话撤销。
- 审计记录：关键身份操作留痕，并对口令和令牌等敏感字段做过滤。
- 后台任务：使用 SQLite 保存待执行任务，支持去重、租约、重试和完成回执。

## 运行环境

- Python 3.11
- SQLite 3，由 Python 标准库提供
- Linux、macOS 或 Windows

## 安装

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install -e ".[dev]"
```

配置项均以 `TOWNSHIP_` 开头。可以复制 `.env.example` 后按需设置，默认数据库位于 `./data/township.db`。

## 初始化与检查

```bash
python -m app.cli init-db
python -m app.cli check-db
```

## 启动服务

```bash
uvicorn app.main:app --host 0.0.0.0 --port 8432
```

健康检查：

```bash
curl -sS http://127.0.0.1:8432/api/system/health
```

首次部署可创建唯一的初始管理员：

```bash
curl -sS -X POST http://127.0.0.1:8432/api/auth/bootstrap   -H 'Content-Type: application/json'   -d '{"username":"admin","password":"Admin!23456","client_label":"initial-setup"}'
```

之后通过 `/api/auth/login` 获取会话令牌，并在管理接口请求头中使用 `Authorization: Bearer <token>`。

## 混合反演区间估计

点估计只给出一个最优比例，端元特征接近时会掩盖很大的不确定性。区间估计作为独立的可恢复任务挂在已完成的点估计任务之后，结果写入独立的区间任务表，**不会修改或覆盖原始点估计**：

1. 先完成点估计：`POST /api/hydro/samples/{sample_id}/inversions` 入队，`POST /api/hydro/inversions/{task_id}/run` 执行。
2. 入队区间任务：`POST /api/hydro/inversions/{task_id}/intervals`，选择方法与配置，再由 `POST /api/hydro/intervals/{id}/run?worker_id=...` 执行；任务中断后可直接重跑。
3. 查询与比较：`GET /api/hydro/intervals/{id}` 取单个结果，`GET /api/hydro/intervals?sample_id=&model_version=&confidence_level=&method=` 按样本、模型版本、置信水平或方法过滤比较。

支持两种方法：

- `parametric-bootstrap`（默认）：按样本分析误差（`measurement_error`、`detection_limit`）与端元 `uncertainty` 的正态模型做参数自助，`random_seed` 固定时同一配置在任意机器、任意重试次数下逐位复现；返回 type-7 经验分位数区间和端元比例的 Pearson 相关矩阵。
- `deterministic-profile`：不使用随机数，做单参数似然剖面，粗网格扫描后二分细化，阈值取精确的 χ²(1)=z((1+c)/2)²，并用数值 Hessian 给出局部端元相关矩阵。

结果结构对两种方法统一：

- `intervals[]`：每个端元的 `fraction_point`、`lower`、`upper`、`width`、是否包含点估计、是否触及单纯形边界；
- `summary`：`n_samples`（样本/网格点数）、`n_success`、`n_failures`、`success_rate`；
- `correlations`：相关矩阵（含标签顺序）与噪声尺度下的端元标准化距离矩阵（距离 < 3σ 时产生 `endmembers_poorly_separated` 告警，提示端元不可区分）；
- `diagnostics`：尾部分辨率不足、区间触界、恒定比例、欠定模型、点估计未被覆盖等告警。

幂等与可复现保证：区间任务的 `task_key` 由点估计任务、完整区间配置（方法、置信水平、种子、次数/网格点数）和点估计结果快照共同决定。重复入队命中同一任务行；已完成任务重跑直接返回原结果，`attempts` 不增加；失败后重跑由于计算是纯确定性函数（自助法以整数种子派生每个复制的随机子流），结果与首次完全一致。不同置信水平、不同种子或不同模型版本产生独立任务，可并列比较。

## 测试

```bash
python -m pytest
```

测试覆盖身份初始化、登录、用户与角色维护、权限计算、账号停用后的会话撤销、审计脱敏、井点样本、同位素约束、迁移计算、任务恢复和数据库时间格式，以及区间估计的两种方法、幂等入队、失败恢复、逐位可复现、点估计不被覆盖和跨置信水平/模型版本比较。

## 编译检查

```bash
python -m compileall -q app tests
```

## API 冒烟

```bash
python -m app.cli smoke
```

该命令在进程内启动应用并检查服务根路径与健康接口，适合部署前快速确认路由和数据库初始化是否正常。

## 目录结构

```text
app/
  api/             用户、角色、审计、认证和系统接口
  core/            时钟、安全、异常和分页能力
  repositories/    SQLite 查询与持久化读取
  routers/         灾情、事件、公告、部门和信访业务接口
  hydro/            地下水、同位素反演和污染迁移服务
  schemas/         管理接口输入模型
  services/        身份、审计和后台任务领域服务
  cli.py           初始化、检查和冒烟入口
  database.py      SQLite 连接、事务、表结构与基础权限
tests/             核心、管理接口和原有业务回归测试
tools/             本地维护脚本
```

## 数据一致性

SQLite 连接默认启用外键、WAL、busy timeout 与同步写入策略。需要跨多张表更新的管理操作在即时事务中执行，失败会整体回滚。会话令牌只保存摘要；用户停用会撤销仍有效的会话。审计事件保存操作者、动作、资源、结果和前后状态，但不会保存明文密码或令牌。
