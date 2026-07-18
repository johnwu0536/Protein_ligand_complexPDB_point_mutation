# PDBMutator

蛋白质点突变工具 - 基于 PySide6 的 Windows 桌面应用程序。

## 功能特点

- 使用 **Gemmi** 解析 PDB 文件，提取 model、chain、residue 信息
- 使用 **PDBFixer** 执行 20 种标准氨基酸之间的单点突变
- 保留蛋白质主链原子，精确控制侧链原子替换
- 原子碰撞检测（默认阈值 2.0 Å）
- 检查目标残基是否存在、是否为标准氨基酸、主链完整性、altLoc、二硫键
- JSON 格式的突变日志记录
- 后台线程执行突变，界面不卡死
- 所有异常转换为中文提示
- 支持中文路径和文件名

## 系统要求

- Windows 10/11
- Python 3.11 或 3.12

## 安装

### 使用 pip

```bash
pip install -r requirements.txt
```

### 使用 conda

```bash
conda env create -f environment.yml
conda activate pdbutator
```

## 运行

```bash
python -m src.main
```

或

```bash
python src/main.py
```

## 使用说明

1. 点击 **"打开PDB文件"** 按钮选择 PDB 文件
2. 选择 **Model** 和 **Chain**
3. 在残基表格中选择目标残基
4. 在 **"突变至"** 下拉框中选择目标氨基酸
5. 点击 **"执行突变"** 按钮
6. 选择输出文件路径，等待突变完成

## 项目结构

```
point_mutation/
├── src/
│   ├── main.py              # 程序入口
│   ├── core/                # 业务逻辑层
│   │   ├── models.py        # 数据模型 (dataclass)
│   │   ├── parser.py        # PDB 解析器 (Gemmi)
│   │   ├── engine.py        # 突变引擎 (BaseMutationEngine + PDBFixer)
│   │   ├── validator.py     # 输入验证与碰撞检测
│   │   └── logger.py        # JSON 日志记录
│   ├── gui/                 # GUI 层 (PySide6)
│   │   ├── main_window.py   # 主窗口
│   │   ├── widgets.py       # 自定义组件
│   │   └── worker.py        # 后台工作线程
│   └── utils/
│       └── exceptions.py    # 自定义异常 (中文提示)
├── tests/                   # 单元测试
│   ├── test_parser.py
│   ├── test_engine.py
│   └── test_validator.py
├── requirements.txt         # pip 依赖
├── environment.yml          # conda 环境配置
└── PDBMutator.spec          # PyInstaller 打包配置
```

## 运行测试

```bash
python -m pytest tests/ -v
```

## 打包为可执行文件

```bash
pip install pyinstaller
pyinstaller PDBMutator.spec --clean
```

打包后的可执行文件位于 `dist/PDBMutator.exe`（约 178 MB，包含所有依赖）。

> **注意**: 如果遇到 `PermissionError`，请先删除旧的 `dist/` 和 `build/` 目录再重新打包。
> 可使用 `--distpath` 和 `--workpath` 参数指定不同的输出目录以避免文件锁定问题。

## 扩展开发

`BaseMutationEngine` 抽象基类定义了突变接口，后续可以扩展：

- **OpenMM 引擎**: 使用 OpenMM 进行能量最小化后突变
- **MODELLER 引擎**: 使用 MODELLER 进行同源建模
- **PyRosetta 引擎**: 使用 PyRosetta 进行构象优化

## 许可证

MIT License
