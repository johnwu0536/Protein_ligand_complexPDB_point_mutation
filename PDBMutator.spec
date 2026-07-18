# -*- mode: python ; coding: utf-8 -*-
"""
PDBMutator PyInstaller 打包配置。
"""

import sys
from pathlib import Path

block_cipher = None

# 项目根目录（spec 文件所在目录）
root_dir = Path('.')

# 收集 PySide6 插件目录
import PySide6
pyside6_dir = Path(PySide6.__file__).parent

# 收集 OpenMM 数据文件
import openmm.app
openmm_app_dir = Path(openmm.app.__file__).parent
openmm_data_dir = openmm_app_dir / 'data'

# 收集 PDBFixer 模板目录
import pdbfixer
pdbfixer_dir = Path(pdbfixer.__file__).parent
pdbfixer_templates_dir = pdbfixer_dir / 'templates'

# 手动递归收集 OpenMM data 目录下的所有文件
# pdbfile.py 使用 os.path.join(os.path.dirname(__file__), 'data', 'pdbNames.xml')
# 所以 data 目录需要放在 openmm/app/data/ 下
openmm_data_files = []
for f in openmm_data_dir.rglob('*'):
    if f.is_file():
        rel_path = f.relative_to(openmm_app_dir)
        openmm_data_files.append((str(f), str(Path('openmm/app') / rel_path.parent)))

# 手动递归收集 PDBFixer templates 目录下的所有文件
# pdbfixer.py 使用 os.path.join(os.path.dirname(__file__), 'templates')
# 所以 templates 目录需要放在 pdbfixer/templates/ 下
pdbfixer_template_files = []
for f in pdbfixer_templates_dir.rglob('*'):
    if f.is_file():
        rel_path = f.relative_to(pdbfixer_dir)
        pdbfixer_template_files.append((str(f), str(Path('pdbfixer') / rel_path.parent)))

a = Analysis(
    ['src\\main.py'],
    pathex=[str(root_dir)],
    binaries=[],
    datas=[
        (str(pyside6_dir / 'plugins'), 'PySide6/plugins'),
    ] + openmm_data_files + pdbfixer_template_files,
    hiddenimports=[
        # 项目模块
        'src',
        'src.core',
        'src.core.models',
        'src.core.parser',
        'src.core.validator',
        'src.core.engine',
        'src.core.logger',
        'src.gui',
        'src.gui.main_window',
        'src.gui.widgets',
        'src.gui.worker',
        'src.utils',
        'src.utils.exceptions',
        # Gemmi
        'gemmi',
        # PDBFixer & OpenMM
        'pdbfixer',
        'openmm',
        'openmm.app',
        'openmm.vec3',
    ],
    hookspath=[],
    hooksconfig={},
    runtime_hooks=[],
    excludes=[
        'tkinter',
        'matplotlib',
        'PIL',
        'pandas',
        'scipy',
        'notebook',
        'jupyter',
        'ipython',
        'setuptools',
        'pip',
        'pytest',
        'wheel',
    ],
    win_no_prefer_redirects=False,
    win_private_assemblies=False,
    cipher=block_cipher,
    noarchive=False,
)

pyz = PYZ(a.pure, a.zipped_data, cipher=block_cipher)

exe = EXE(
    pyz,
    a.scripts,
    a.binaries,
    a.zipfiles,
    a.datas,
    [],
    name='PDBMutator',
    debug=False,
    bootloader_ignore_signals=False,
    strip=False,
    upx=True,
    upx_exclude=[],
    runtime_tmpdir=None,
    console=False,          # 不显示控制台窗口
    disable_windowed_traceback=False,
    argv_emulation=False,
    target_arch=None,
    codesign_identity=None,
    entitlements_file=None,
    icon=None,              # 可添加图标文件: 'PDBMutator.ico'
)
