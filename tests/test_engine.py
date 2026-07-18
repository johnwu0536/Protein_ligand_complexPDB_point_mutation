"""
突变引擎单元测试。
"""

from __future__ import annotations

from pathlib import Path

import pytest

from src.core.engine import BaseMutationEngine, PdbFixerMutationEngine
from src.core.models import MutationResult, ResidueInfo, AtomInfo
from src.core.parser import PdbParser


class TestBaseMutationEngine:
    """BaseMutationEngine 测试类。"""

    def test_abstract_class_cannot_instantiate(self):
        """测试抽象类不能直接实例化。"""
        with pytest.raises(TypeError):
            BaseMutationEngine()  # type: ignore

    def test_concrete_class_can_instantiate(self):
        """测试具体类可以实例化。"""
        engine = PdbFixerMutationEngine()
        assert isinstance(engine, BaseMutationEngine)
        assert isinstance(engine, PdbFixerMutationEngine)


class TestPdbFixerMutationEngine:
    """PdbFixerMutationEngine 测试类。"""

    @pytest.fixture
    def engine(self) -> PdbFixerMutationEngine:
        return PdbFixerMutationEngine()

    @pytest.fixture
    def parser_with_ala(self, tmp_path: Path) -> PdbParser:
        """创建包含 ALA 残基的解析器。"""
        content = """\
ATOM      1  N   ALA A   1       1.000   2.000   3.000  1.00  0.00           N
ATOM      2  CA  ALA A   1       2.000   2.000   3.000  1.00  0.00           C
ATOM      3  C   ALA A   1       3.000   2.000   3.000  1.00  0.00           C
ATOM      4  O   ALA A   1       4.000   2.000   3.000  1.00  0.00           O
ATOM      5  CB  ALA A   1       1.500   1.000   3.000  1.00  0.00           C
END
"""
        pdb_path = tmp_path / "ala.pdb"
        pdb_path.write_text(content, encoding="utf-8")
        parser = PdbParser()
        parser.parse(pdb_path)
        return parser

    def test_mutate_returns_mutation_result(self, engine, parser_with_ala, tmp_path):
        """测试 mutate 返回 MutationResult 类型。"""
        target = parser_with_ala.get_residue(1, "A", 1)
        assert target is not None

        output_path = tmp_path / "output.pdb"
        result = engine.mutate(parser_with_ala, target, "GLY", output_path)

        assert isinstance(result, MutationResult)

    def test_mutate_does_not_overwrite_input(self, engine, parser_with_ala, tmp_path):
        """测试突变不会覆盖输入文件。"""
        target = parser_with_ala.get_residue(1, "A", 1)
        assert target is not None

        input_path = parser_with_ala.path
        output_path = tmp_path / "output.pdb"

        engine.mutate(parser_with_ala, target, "GLY", output_path)

        # 输入文件应该保持不变
        assert input_path.exists()
        assert output_path.exists()
        assert input_path != output_path

    def test_mutate_with_invalid_target(self, engine, parser_with_ala, tmp_path):
        """测试无效目标的突变。"""
        # 创建一个不存在的残基
        target = ResidueInfo(
            model=1,
            chain="A",
            residue_number=999,
            residue_name="ALA",
        )
        output_path = tmp_path / "output.pdb"
        result = engine.mutate(parser_with_ala, target, "GLY", output_path)

        # 应该返回失败结果而不是抛出异常
        assert result.success is False
        assert result.error_message is not None
