"""
PDBMutator 数据模型模块。
使用 dataclass 定义核心数据结构。
"""

from __future__ import annotations

from dataclasses import dataclass, field
from datetime import datetime
from enum import Enum
from pathlib import Path
from typing import Optional


# 20种标准氨基酸的三字母代码
STANDARD_AMINO_ACIDS: list[str] = [
    "ALA", "ARG", "ASN", "ASP", "CYS",
    "GLN", "GLU", "GLY", "HIS", "ILE",
    "LEU", "LYS", "MET", "PHE", "PRO",
    "SER", "THR", "TRP", "TYR", "VAL",
]

# 氨基酸三字母到中文名称的映射
AMINO_ACID_NAMES_CN: dict[str, str] = {
    "ALA": "丙氨酸", "ARG": "精氨酸", "ASN": "天冬酰胺", "ASP": "天冬氨酸",
    "CYS": "半胱氨酸", "GLN": "谷氨酰胺", "GLU": "谷氨酸", "GLY": "甘氨酸",
    "HIS": "组氨酸", "ILE": "异亮氨酸", "LEU": "亮氨酸", "LYS": "赖氨酸",
    "MET": "甲硫氨酸", "PHE": "苯丙氨酸", "PRO": "脯氨酸", "SER": "丝氨酸",
    "THR": "苏氨酸", "TRP": "色氨酸", "TYR": "酪氨酸", "VAL": "缬氨酸",
}

# 主链原子名称
BACKBONE_ATOMS: list[str] = ["N", "CA", "C", "O"]


@dataclass(frozen=True)
class AtomInfo:
    """原子信息。"""
    name: str                    # 原子名称，如 "CA", "CB", "CG"
    element: str                 # 元素符号，如 "C", "N", "O", "S"
    x: float                     # X 坐标
    y: float                     # Y 坐标
    z: float                     # Z 坐标
    is_backbone: bool = False    # 是否为主链原子
    alt_loc: str = ""            # 交替位置标识


@dataclass(frozen=True)
class ResidueInfo:
    """残基信息。"""
    model: int                   # Model 编号
    chain: str                   # Chain ID
    residue_number: int          # 残基编号
    insertion_code: str = ""     # 插入码
    residue_name: str = ""       # 残基名称（三字母代码）
    alt_loc: str = ""            # 交替位置标识
    atoms: tuple[AtomInfo, ...] = field(default_factory=tuple)  # 原子列表

    @property
    def is_standard(self) -> bool:
        """是否为20种标准氨基酸。"""
        return self.residue_name in STANDARD_AMINO_ACIDS

    @property
    def has_backbone(self) -> bool:
        """主链原子是否完整。"""
        backbone_names = set(a.name for a in self.atoms if a.is_backbone)
        return all(b in backbone_names for b in BACKBONE_ATOMS)

    @property
    def missing_backbone_atoms(self) -> list[str]:
        """缺失的主链原子列表。"""
        backbone_names = set(a.name for a in self.atoms if a.is_backbone)
        return [b for b in BACKBONE_ATOMS if b not in backbone_names]

    @property
    def has_alt_loc(self) -> bool:
        """是否存在交替位置。"""
        return any(a.alt_loc for a in self.atoms)

    @property
    def alt_loc_values(self) -> list[str]:
        """获取所有 altLoc 值。"""
        locs = set(a.alt_loc for a in self.atoms if a.alt_loc)
        return sorted(locs)

    def display_name(self) -> str:
        """返回可读的残基名称。"""
        cn = AMINO_ACID_NAMES_CN.get(self.residue_name, self.residue_name)
        return f"{self.residue_name}({cn})"


@dataclass
class MutationRecord:
    """突变记录，用于 JSON 日志。"""
    input_file: str              # 输入文件路径
    output_file: str             # 输出文件路径
    model: int                   # Model 编号
    chain: str                   # Chain ID
    residue_number: int          # 残基编号
    insertion_code: str          # 插入码
    original_residue: str        # 原氨基酸三字母代码
    new_residue: str             # 新氨基酸三字母代码
    timestamp: str               # 时间戳
    warnings: list[str] = field(default_factory=list)  # 警告列表
    clashes: list[dict] = field(default_factory=list)  # 冲突列表
    success: bool = True         # 是否成功


@dataclass
class MutationResult:
    """突变结果。"""
    success: bool                # 是否成功
    output_path: Optional[Path] = None  # 输出文件路径
    warnings: list[str] = field(default_factory=list)  # 警告列表
    clashes: list[tuple[str, str, float]] = field(default_factory=list)  # 冲突列表
    error_message: Optional[str] = None  # 错误消息


@dataclass
class StructureInfo:
    """PDB 结构概要信息。"""
    path: Path                   # 文件路径
    models: list[int] = field(default_factory=list)  # model 编号列表
    chains: dict[int, list[str]] = field(default_factory=dict)  # model -> chain 列表
    residues: dict[tuple[int, str], list[ResidueInfo]] = field(default_factory=dict)  # (model, chain) -> 残基列表
    disulfide_bonds: list[tuple[str, int, str, int]] = field(default_factory=list)  # (chain1, seq1, chain2, seq2)
