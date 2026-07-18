"""
PDBMutator 验证器模块。
负责输入验证和原子碰撞检测。
"""

from __future__ import annotations

import logging
import math
from typing import Optional

from src.core.models import (
    BACKBONE_ATOMS,
    ResidueInfo,
    STANDARD_AMINO_ACIDS,
    StructureInfo,
)
from src.core.parser import PdbParser
from src.utils.exceptions import (
    AltLocPresentError,
    ClashError,
    DisulfideBondError,
    MissingBackboneError,
    NonStandardResidueError,
    ResidueNotFoundError,
)

logger = logging.getLogger(__name__)


class MutationValidator:
    """突变验证器，负责检查目标残基的合法性。"""

    # 非键合重原子距离阈值（Å）
    DEFAULT_CLASH_THRESHOLD: float = 2.0

    def __init__(self, parser: PdbParser) -> None:
        """
        初始化验证器。

        Args:
            parser: 已加载 PDB 文件的解析器实例
        """
        self._parser = parser

    def validate_target(
        self,
        model_id: int,
        chain_id: str,
        residue_number: int,
        insertion_code: str = "",
        check_disulfide: bool = True,
    ) -> ResidueInfo:
        """
        综合验证目标残基。

        Args:
            model_id: Model 编号
            chain_id: Chain ID
            residue_number: 残基编号
            insertion_code: 插入码
            check_disulfide: 是否检查二硫键

        Returns:
            ResidueInfo: 验证通过的残基信息

        Raises:
            ResidueNotFoundError: 残基不存在
            NonStandardResidueError: 非标准氨基酸
            MissingBackboneError: 主链原子缺失
            AltLocPresentError: 存在 altLoc
            DisulfideBondError: 参与二硫键
        """
        # 1. 检查残基是否存在
        residue = self._parser.get_residue(
            model_id, chain_id, residue_number, insertion_code
        )
        if residue is None:
            raise ResidueNotFoundError(
                model_id, chain_id, residue_number, insertion_code
            )

        # 2. 检查是否为标准氨基酸
        if not residue.is_standard:
            raise NonStandardResidueError(residue.residue_name)

        # 3. 检查主链原子是否完整
        if not residue.has_backbone:
            raise MissingBackboneError(
                residue.residue_name,
                residue.residue_number,
                residue.missing_backbone_atoms,
            )

        # 4. 检查是否存在 altLoc
        if residue.has_alt_loc:
            raise AltLocPresentError(
                residue.residue_name,
                residue.residue_number,
                residue.alt_loc_values,
            )

        # 5. 检查是否参与二硫键
        if check_disulfide and residue.residue_name == "CYS":
            disulfide_partner = self._parser.check_disulfide(
                chain_id, residue_number
            )
            if disulfide_partner is not None:
                partner_str = f"Chain {disulfide_partner[0]}, Residue {disulfide_partner[1]}"
                raise DisulfideBondError(
                    residue.residue_name,
                    residue.residue_number,
                    partner_str,
                )

        return residue

    def check_clashes(
        self,
        target: ResidueInfo,
        new_residue_name: str,
        threshold: float = DEFAULT_CLASH_THRESHOLD,
    ) -> list[tuple[str, str, float]]:
        """
        检查突变后新残基与周围原子的碰撞。

        通过计算新残基重原子与周围所有非键合重原子的距离，
        将距离小于阈值的标记为冲突。

        Args:
            target: 目标残基信息
            new_residue_name: 新氨基酸名称
            threshold: 距离阈值（Å），默认 2.0

        Returns:
            list[tuple[str, str, float]]: 冲突列表，每项为
                (原子1名称, 原子2名称, 距离)
        """
        clashes: list[tuple[str, str, float]] = []

        # 获取目标残基周围的所有原子（同一 chain 中相邻残基）
        # 这里简化处理：检查目标残基与同一 chain 中所有其他残基的碰撞
        model_id = target.model
        chain_id = target.chain
        all_residues = self._parser.get_residues(model_id, chain_id)

        # 获取新残基的侧链原子名称（排除主链原子）
        # 注意：这里我们使用目标残基的侧链原子作为近似
        # 实际突变后侧链原子会变化，但位置未知
        # 因此我们检查目标残基现有侧链原子与周围原子的距离
        sidechain_atoms = [
            a for a in target.atoms if not a.is_backbone
        ]

        for other_res in all_residues:
            if other_res.residue_number == target.residue_number:
                continue  # 跳过自身

            for sc_atom in sidechain_atoms:
                for other_atom in other_res.atoms:
                    # 跳过氢原子（仅检查重原子）
                    if other_atom.element == "H":
                        continue

                    dist = self._calc_distance(sc_atom, other_atom)
                    if dist < threshold:
                        clashes.append((
                            f"{target.residue_name}{target.residue_number}:{sc_atom.name}",
                            f"{other_res.residue_name}{other_res.residue_number}:{other_atom.name}",
                            dist,
                        ))

        return clashes

    def _calc_distance(
        self,
        atom1: "AtomInfo",
        atom2: "AtomInfo",
    ) -> float:
        """计算两个原子之间的欧几里得距离。"""
        dx = atom1.x - atom2.x
        dy = atom1.y - atom2.y
        dz = atom1.z - atom2.z
        return math.sqrt(dx * dx + dy * dy + dz * dz)

    def validate_and_check(
        self,
        model_id: int,
        chain_id: str,
        residue_number: int,
        new_residue_name: str,
        insertion_code: str = "",
        clash_threshold: float = DEFAULT_CLASH_THRESHOLD,
    ) -> tuple[ResidueInfo, list[tuple[str, str, float]]]:
        """
        综合验证并检查碰撞。

        Args:
            model_id: Model 编号
            chain_id: Chain ID
            residue_number: 残基编号
            new_residue_name: 新氨基酸名称
            insertion_code: 插入码
            clash_threshold: 碰撞检测阈值

        Returns:
            tuple[ResidueInfo, list]: (残基信息, 冲突列表)

        Raises:
            各种验证异常
        """
        # 验证目标残基
        residue = self.validate_target(
            model_id, chain_id, residue_number, insertion_code
        )

        # 检查新氨基酸是否合法
        if new_residue_name not in STANDARD_AMINO_ACIDS:
            raise NonStandardResidueError(new_residue_name)

        # 检查碰撞
        clashes = self.check_clashes(residue, new_residue_name, clash_threshold)

        return residue, clashes
