"""
PDBMutator PDB 解析器模块。
使用 Gemmi 库解析 PDB 文件并提取结构信息。
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Optional

import gemmi

from src.core.models import (
    AtomInfo,
    BACKBONE_ATOMS,
    ResidueInfo,
    STANDARD_AMINO_ACIDS,
    StructureInfo,
)
from src.utils.exceptions import PdbParseError

logger = logging.getLogger(__name__)


class PdbParser:
    """PDB 文件解析器，封装 Gemmi 操作。"""

    def __init__(self) -> None:
        self._structure: Optional[gemmi.Structure] = None
        self._path: Optional[Path] = None

    @property
    def structure(self) -> gemmi.Structure:
        """获取 Gemmi Structure 对象。"""
        if self._structure is None:
            raise RuntimeError("请先调用 parse() 方法加载PDB文件")
        return self._structure

    @property
    def path(self) -> Path:
        """获取文件路径。"""
        if self._path is None:
            raise RuntimeError("请先调用 parse() 方法加载PDB文件")
        return self._path

    def parse(self, path: Path | str) -> StructureInfo:
        """
        解析 PDB 文件。

        Args:
            path: PDB 文件路径（支持中文路径）

        Returns:
            StructureInfo: 结构概要信息

        Raises:
            PdbParseError: 解析失败时抛出
        """
        path = Path(path)
        if not path.exists():
            raise PdbParseError(path, "文件不存在")

        try:
            self._structure = gemmi.read_structure(str(path))
            self._path = path
        except Exception as e:
            raise PdbParseError(path, str(e))

        return self._extract_structure_info()

    def _extract_structure_info(self) -> StructureInfo:
        """从 Gemmi Structure 提取概要信息。"""
        info = StructureInfo(path=self._path)

        for model in self.structure:
            model_num = model.num if model.num else 1
            try:
                model_num = int(model_num)
            except (ValueError, TypeError):
                model_num = 1

            info.models.append(model_num)
            info.chains[model_num] = []

            for chain in model:
                chain_id = chain.name
                info.chains[model_num].append(chain_id)

                residues = []
                for residue in chain:
                    if residue.name.strip() in ("HOH", "WAT"):
                        continue

                    res_info = self._extract_residue_info(
                        residue, model_num, chain_id
                    )
                    if res_info is not None:
                        residues.append(res_info)

                info.residues[(model_num, chain_id)] = residues

        # 提取二硫键信息
        try:
            for conn in self.structure.connections:
                if conn.type == gemmi.ConnectionType.Disulf:
                    if len(conn.partner1) >= 2 and len(conn.partner2) >= 2:
                        c1 = conn.partner1[0].chain_name
                        s1 = conn.partner1[1].seq_id.num
                        c2 = conn.partner2[0].chain_name
                        s2 = conn.partner2[1].seq_id.num
                        info.disulfide_bonds.append((c1, s1, c2, s2))
        except Exception:
            logger.warning("提取二硫键信息时出现异常", exc_info=True)

        return info

    def _extract_residue_info(
        self,
        residue: gemmi.Residue,
        model_num: int,
        chain_id: str,
    ) -> Optional[ResidueInfo]:
        """从 Gemmi Residue 提取残基信息。"""
        res_name = residue.name.strip()
        if not res_name:
            return None

        seq_id = residue.seqid
        res_num = seq_id.num
        ins_code = seq_id.icode.strip() if seq_id.icode else ""

        # 收集原子信息
        atoms = []
        alt_locs = set()
        for atom in residue:
            atom_name = atom.name.strip()
            element = atom.element.name if atom.element else ""
            if not element:
                element = atom_name[0] if atom_name else ""

            # 判断是否为主链原子
            is_backbone = atom_name in BACKBONE_ATOMS

            alt_loc = atom.altloc.strip() if atom.altloc and atom.altloc != '\x00' else ""
            if alt_loc:
                alt_locs.add(alt_loc)

            atom_info = AtomInfo(
                name=atom_name,
                element=element,
                x=atom.pos.x,
                y=atom.pos.y,
                z=atom.pos.z,
                is_backbone=is_backbone,
                alt_loc=alt_loc,
            )
            atoms.append(atom_info)

        if not atoms:
            return None

        alt_loc_str = ""
        if alt_locs:
            alt_loc_str = ",".join(sorted(alt_locs))

        return ResidueInfo(
            model=model_num,
            chain=chain_id,
            residue_number=res_num,
            insertion_code=ins_code,
            residue_name=res_name,
            alt_loc=alt_loc_str,
            atoms=tuple(atoms),
        )

    def get_models(self) -> list[int]:
        """获取所有 model 编号。"""
        models = []
        for model in self.structure:
            model_num = model.num if model.num else 1
            try:
                model_num = int(model_num)
            except (ValueError, TypeError):
                model_num = 1
            models.append(model_num)
        return models

    def get_chains(self, model_id: int) -> list[str]:
        """获取指定 model 的所有 chain ID。"""
        chains = []
        for model in self.structure:
            model_num = model.num if model.num else 1
            try:
                model_num = int(model_num)
            except (ValueError, TypeError):
                model_num = 1
            if model_num == model_id:
                for chain in model:
                    chains.append(chain.name)
                break
        return chains

    def get_residues(self, model_id: int, chain_id: str) -> list[ResidueInfo]:
        """获取指定 model 和 chain 的所有残基。"""
        residues = []
        for model in self.structure:
            model_num = model.num if model.num else 1
            try:
                model_num = int(model_num)
            except (ValueError, TypeError):
                model_num = 1
            if model_num == model_id:
                for chain in model:
                    if chain.name == chain_id:
                        for residue in chain:
                            if residue.name.strip() in ("HOH", "WAT"):
                                continue
                            res_info = self._extract_residue_info(
                                residue, model_num, chain_id
                            )
                            if res_info is not None:
                                residues.append(res_info)
                        break
                break
        return residues

    def get_residue(
        self,
        model_id: int,
        chain_id: str,
        residue_number: int,
        insertion_code: str = "",
    ) -> Optional[ResidueInfo]:
        """获取指定残基的详细信息。"""
        for model in self.structure:
            model_num = model.num if model.num else 1
            try:
                model_num = int(model_num)
            except (ValueError, TypeError):
                model_num = 1
            if model_num == model_id:
                for chain in model:
                    if chain.name == chain_id:
                        for residue in chain:
                            seq_id = residue.seqid
                            res_num = seq_id.num
                            ins_code = seq_id.icode.strip() if seq_id.icode else ""
                            if res_num == residue_number and ins_code == insertion_code:
                                return self._extract_residue_info(
                                    residue, model_num, chain_id
                                )
                        break
                break
        return None

    def get_gemmi_residue(
        self,
        model_id: int,
        chain_id: str,
        residue_number: int,
        insertion_code: str = "",
    ) -> Optional[gemmi.Residue]:
        """获取 Gemmi 原生 Residue 对象（用于 PDBFixer 操作）。"""
        for model in self.structure:
            model_num = model.num if model.num else 1
            try:
                model_num = int(model_num)
            except (ValueError, TypeError):
                model_num = 1
            if model_num == model_id:
                for chain in model:
                    if chain.name == chain_id:
                        for residue in chain:
                            seq_id = residue.seqid
                            res_num = seq_id.num
                            ins_code = seq_id.icode.strip() if seq_id.icode else ""
                            if res_num == residue_number and ins_code == insertion_code:
                                return residue
                        break
                break
        return None

    def check_disulfide(
        self,
        chain_id: str,
        residue_number: int,
    ) -> Optional[tuple[str, int]]:
        """检查残基是否参与二硫键。返回配对残基的 (chain, seq_num) 或 None。"""
        for conn in self.structure.connections:
            if conn.type == gemmi.ConnectionType.Disulf:
                if len(conn.partner1) >= 2 and len(conn.partner2) >= 2:
                    c1 = conn.partner1[0].chain_name
                    s1 = conn.partner1[1].seq_id.num
                    c2 = conn.partner2[0].chain_name
                    s2 = conn.partner2[1].seq_id.num
                    if c1 == chain_id and s1 == residue_number:
                        return (c2, s2)
                    if c2 == chain_id and s2 == residue_number:
                        return (c1, s1)
        return None
