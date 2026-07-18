# PDBMutator 开发计划

## 一、项目概述

PDBMutator 是一个 Windows 本地桌面程序，用于对 PDB 文件中的蛋白质残基进行单点突变。它使用 Gemmi 解析 PDB 文件，使用 PDBFixer 执行突变操作，并提供 PySide6 图形界面。

## 二、目录结构

```
point_mutation/
├── requirements.txt              # pip 依赖
├── environment.yml               # conda 环境配置
├── README.md                     # 项目说明文档
├── PDBMutator.spec               # PyInstaller 打包配置
├── DEVELOPMENT_PLAN.md           # 本开发计划文件
│
├── src/
│   ├── __init__.py
│   ├── main.py                   # 程序入口
│   │
│   ├── core/                     # 业务逻辑层
│   │   ├── __init__.py
│   │   ├── models.py             # 数据模型 (dataclass)
│   │   ├── parser.py             # PDB 解析器 (Gemmi)
│   │   ├── engine.py             # 突变引擎接口与实现
│   │   ├── validator.py          # 输入验证与冲突检测
│   │   └── logger.py             # JSON 日志记录
│   │
│   ├── gui/                      # GUI 层
│   │   ├── __init__.py
│   │   ├── main_window.py        # 主窗口
│   │   ├── widgets.py            # 自定义组件
│   │   └── worker.py             # 后台工作线程
│   │
│   └── utils/                    # 工具模块
│       ├── __init__.py
│       └── exceptions.py         # 自定义异常与中文提示
│
└── tests/                        # 测试模块
    ├── __init__.py
    ├── test_parser.py            # 解析器测试
    ├── test_engine.py            # 突变引擎测试
    └── test_validator.py         # 验证器测试
```

## 三、模块详细设计

### 1. 数据模型 (`src/core/models.py`)

使用 dataclass 定义：

- `ResidueInfo`: 残基信息（model, chain, residue_number, insertion_code, residue_name, alt_loc）
- `MutationRecord`: 突变记录（输入文件、输出文件、原残基信息、新氨基酸、时间、警告）
- `MutationResult`: 突变结果（成功/失败、输出路径、警告列表、冲突列表）
- `AtomInfo`: 原子信息（名称、坐标、元素、是否主链原子）

### 2. PDB 解析器 (`src/core/parser.py`)

使用 Gemmi 库：

- `PdbParser` 类
  - `parse(path)`: 解析 PDB 文件，返回结构信息
  - `get_models()`: 获取所有 model 编号
  - `get_chains(model_id)`: 获取指定 model 的所有 chain
  - `get_residues(model_id, chain_id)`: 获取指定 chain 的所有残基
  - `get_residue_info(model_id, chain_id, residue_number, insertion_code)`: 获取单个残基详细信息
  - `check_backbone_atoms(residue)`: 检查主链原子 (N, CA, C, O) 是否完整
  - `check_alt_loc(residue)`: 检查是否存在 altLoc
  - `check_disulfide(residue, structure)`: 检查是否参与二硫键
  - `extract_atoms(residue)`: 提取残基原子信息

### 3. 突变引擎 (`src/core/engine.py`)

- `BaseMutationEngine` (抽象基类): 定义突变接口
  - `mutate(structure, target, new_residue_name) -> MutationResult`
- `PdbFixerMutationEngine`: 使用 PDBFixer 的实现
  - 保留主链原子 (N, CA, C, O)
  - 删除原侧链中不属于新残基的原子
  - 补充新残基缺失的重原子
  - 不添加氢原子（保持 PDB 格式兼容）

### 4. 验证器 (`src/core/validator.py`)

- `MutationValidator` 类
  - `validate_target(structure, target)`: 综合验证目标残基
    - 检查残基是否存在
    - 检查是否为标准 20 种氨基酸
    - 检查主链原子是否完整
    - 检查是否存在 altLoc
    - 检查是否参与二硫键
  - `check_clashes(structure, target, new_residue_name, threshold=2.0)`: 原子碰撞检查
    - 计算新残基重原子与周围原子距离
    - 距离 < 2.0 Å 标记为严重冲突

### 5. 日志模块 (`src/core/logger.py`)

- `MutationLogger` 类
  - `log(mutation_record)`: 将突变记录写入 JSON 文件
  - 日志文件名格式: `mutation_YYYYMMDD_HHMMSS.json`
  - 包含: 输入文件、输出文件、model、chain、residue_number、insertion_code、原氨基酸、新氨基酸、时间、警告列表

### 6. 异常处理 (`src/utils/exceptions.py`)

自定义异常类，所有异常消息为中文：

- `PdbParseError`: PDB 文件解析失败
- `ResidueNotFoundError`: 目标残基不存在
- `NonStandardResidueError`: 非标准氨基酸
- `MissingBackboneError`: 主链原子缺失
- `AltLocPresentError`: 存在 altLoc 需要处理
- `DisulfideBondError`: 参与二硫键
- `MutationError`: 突变执行失败
- `ClashError`: 原子碰撞冲突

### 7. GUI 层 (`src/gui/`)

- `MainWindow`: 主窗口
  - 菜单栏/工具栏
  - 打开 PDB 按钮
  - Model 选择框 (QComboBox)
  - Chain 选择框 (QComboBox)
  - Residue 表格 (QTableWidget): 显示残基编号、名称、主链状态、altLoc、二硫键
  - 新氨基酸选择框 (QComboBox, 20种标准氨基酸)
  - 执行突变按钮
  - 保存 PDB 按钮
  - 日志窗口 (QTextEdit, 只读)
- `Worker`: QThread 子类，后台执行突变操作，避免界面卡死

### 8. 入口 (`src/main.py`)

- 创建 QApplication
- 实例化 MainWindow
- 显示窗口并进入事件循环

## 四、开发阶段

### 第一阶段：项目骨架与依赖安装
- [x] 创建目录结构
- [x] 编写 requirements.txt、environment.yml
- [x] 安装依赖

### 第二阶段：核心数据模型与异常
- [ ] 实现 `src/core/models.py`
- [ ] 实现 `src/utils/exceptions.py`
- [ ] 运行基础测试

### 第三阶段：PDB 解析器
- [ ] 实现 `src/core/parser.py`
- [ ] 编写 `tests/test_parser.py`
- [ ] 运行解析器测试

### 第四阶段：验证器
- [ ] 实现 `src/core/validator.py`
- [ ] 编写 `tests/test_validator.py`
- [ ] 运行验证器测试

### 第五阶段：突变引擎
- [ ] 实现 `src/core/engine.py` (BaseMutationEngine + PdbFixerMutationEngine)
- [ ] 编写 `tests/test_engine.py`
- [ ] 运行引擎测试

### 第六阶段：日志模块
- [ ] 实现 `src/core/logger.py`

### 第七阶段：GUI 界面
- [ ] 实现 `src/gui/widgets.py`
- [ ] 实现 `src/gui/worker.py`
- [ ] 实现 `src/gui/main_window.py`
- [ ] 实现 `src/main.py`

### 第八阶段：集成测试与打包
- [ ] 完整功能测试
- [ ] 编写 README.md
- [ ] 编写 PDBMutator.spec
- [ ] 最终验证

## 五、技术决策

1. **Python 版本**: 当前环境为 Python 3.12，与 3.11 高度兼容，直接使用
2. **PDBFixer 使用方式**: 通过其 Python API 调用，而非命令行
3. **Gemmi 使用方式**: 使用 Gemmi 的 Python 绑定读取 PDB 结构
4. **原子操作策略**: 
   - 主链原子 (N, CA, C, O) 完全保留
   - 使用 PDBFixer 的 `apply_mutation` 方法执行突变
   - 手动删除/添加原子以精确控制
5. **碰撞检测**: 基于原子坐标的欧几里得距离计算
6. **线程模型**: 使用 QThread + 信号槽机制避免 GUI 卡死
7. **路径处理**: 全程使用 pathlib.Path，支持中文路径

## 六、20种标准氨基酸三字母代码

Ala, Arg, Asn, Asp, Cys, Gln, Glu, Gly, His, Ile, Leu, Lys, Met, Phe, Pro, Ser, Thr, Trp, Tyr, Val
