"""
PDBMutator 突变引擎模块。
定义抽象基类 BaseMutationEngine 和 PDBFixer 实现。
"""

from __future__ import annotations

import logging
import tempfile
from abc import ABC, abstractmethod
from pathlib import Path
from typing import Optional

import gemmi

from src.core.models import (
    BACKBONE_ATOMS,
    MutationResult,
    ResidueInfo,
    STANDARD_AMINO_ACIDS,
)
from src.core.parser import PdbParser
from src.utils.exceptions import MutationError

logger = logging.getLogger(__name__)


class BaseMutationEngine(ABC):
    """
    突变引擎抽象基类。

    定义突变接口，后续可扩展为 OpenMM、MODELLER、PyRosetta 等引擎。
    """

    @abstractmethod
    def mutate(
        self,
        parser: PdbParser,
        target: ResidueInfo,
        new_residue_name: str,
        output_path: Path,
    ) -> MutationResult:
        """
        执行单点突变。

        Args:
            parser: 已加载 PDB 文件的解析器
            target: 目标残基信息
            new_residue_name: 新氨基酸三字母代码
            output_path: 输出 PDB 文件路径

        Returns:
            MutationResult: 突变结果
        """
        ...


class PdbFixerMutationEngine(BaseMutationEngine):
    """
    使用 PDBFixer 的突变引擎实现。

    策略：
    1. 保留主链原子 (N, CA, C, O)
    2. 删除原侧链中不属于新残基的原子
    3. 补充新残基缺失的重原子
    4. 不添加氢原子
    """

    def mutate(
        self,
        parser: PdbParser,
        target: ResidueInfo,
        new_residue_name: str,
        output_path: Path,
    ) -> MutationResult:
        """
        执行单点突变。

        Args:
            parser: 已加载 PDB 文件的解析器
            target: 目标残基信息
            new_residue_name: 新氨基酸三字母代码
            output_path: 输出 PDB 文件路径

        Returns:
            MutationResult: 突变结果
        """
        warnings: list[str] = []

        try:
            # 使用 PDBFixer 执行突变
            result = self._mutate_with_pdbfixer(
                parser, target, new_residue_name, output_path
            )
            return result
        except Exception as e:
            logger.error("PDBFixer 突变失败", exc_info=True)
            return MutationResult(
                success=False,
                error_message=f"PDBFixer 突变失败: {str(e)}",
                warnings=warnings,
            )

    @staticmethod
    def _setup_openmm_data_path() -> None:
        """
        设置 OpenMM 数据文件路径。

        在 PyInstaller 打包环境中，OpenMM 的数据文件（如 pdbNames.xml）
        会被解压到临时目录，需要设置 OPENMM_DATA_PATH 环境变量指向该目录。
        """
        import os
        import sys

        # 检查是否在 PyInstaller 打包环境中
        if getattr(sys, 'frozen', False):
            # 打包后的临时目录
            base_dir = Path(sys._MEIPASS)
            data_dir = base_dir / 'openmm' / 'app' / 'data'
            if data_dir.exists():
                os.environ['OPENMM_DATA_PATH'] = str(data_dir)
                logger.info(f"设置 OPENMM_DATA_PATH = {data_dir}")
            else:
                logger.warning(f"OpenMM 数据目录不存在: {data_dir}")

    def _mutate_with_pdbfixer(
        self,
        parser: PdbParser,
        target: ResidueInfo,
        new_residue_name: str,
        output_path: Path,
    ) -> MutationResult:
        """使用 PDBFixer 执行突变的核心逻辑。"""
        warnings: list[str] = []

        try:
            from pdbfixer import PDBFixer
            from openmm.app import PDBFile
        except ImportError:
            raise MutationError(
                "未安装 PDBFixer 或 OpenMM。请执行: pip install pdbfixer openmm"
            )

        # 设置 OpenMM 数据文件路径（支持 PyInstaller 打包环境）
        self._setup_openmm_data_path()


        # 将 Gemmi 结构保存到临时文件，供 PDBFixer 读取
        with tempfile.NamedTemporaryFile(
            suffix=".pdb", delete=False, mode="w"
        ) as tmp_file:
            tmp_path = tmp_file.name
            try:
                parser.structure.write_pdb(tmp_path)
            except Exception as e:
                Path(tmp_path).unlink(missing_ok=True)
                raise MutationError(f"无法写入临时PDB文件: {e}")

        try:
            # 使用 PDBFixer 加载
            fixer = PDBFixer(filename=tmp_path)

            # 应用突变
            # PDBFixer.applyMutation 参数:
            #   chainId, residueIndex, newResidueName
            # residueIndex 是 chain 中残基的索引（从0开始）

            # 找到目标残基在 chain 中的索引
            chain_index = -1
            residue_index = -1
            for i, chain in enumerate(fixer.topology.chains()):
                if chain.id == target.chain:
                    chain_index = i
                    for j, residue in enumerate(chain.residues()):
                        if (residue.id == str(target.residue_number) or
                            residue.index == target.residue_number):
                            residue_index = j
                            break
                    break

            if chain_index == -1:
                raise MutationError(f"在 PDBFixer 拓扑中未找到 chain '{target.chain}'")
            if residue_index == -1:
                raise MutationError(
                    f"在 PDBFixer 拓扑中未找到残基 {target.residue_number}"
                )

            # 执行突变
            # PDBFixer.applyMutations(mutations, chain_id)
            # mutations 格式: ["ALA-133-GLY"] 表示将残基133从ALA突变为GLY
            mutation_str = f"{target.residue_name}-{target.residue_number}-{new_residue_name}"
            fixer.applyMutations([mutation_str], target.chain)

            # 写入输出文件
            with open(str(output_path), "w") as f:
                PDBFile.writeFile(
                    fixer.topology,
                    fixer.positions,
                    f,
                    keepIds=True,
                )

            logger.info(
                f"突变成功: {target.residue_name}{target.residue_number} "
                f"-> {new_residue_name}, 输出: {output_path}"
            )

            return MutationResult(
                success=True,
                output_path=output_path,
                warnings=warnings,
            )

        except Exception as e:
            error_msg = str(e)
            logger.error(f"PDBFixer 突变失败: {error_msg}", exc_info=True)
            return MutationResult(
                success=False,
                error_message=f"PDBFixer 突变失败: {error_msg}",
                warnings=warnings,
            )
        finally:
            # 清理临时文件
            Path(tmp_path).unlink(missing_ok=True)
