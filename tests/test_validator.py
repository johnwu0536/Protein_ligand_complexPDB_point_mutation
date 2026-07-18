"""
验证器单元测试。
"""

from __future__ import annotations

from pathlib import Path

import pytest

from src.core.models import ResidueInfo, AtomInfo
from src.core.parser import PdbParser
from src.core.validator import MutationValidator
from src.utils.exceptions import (
    AltLocPresentError,
    DisulfideBondError,
    MissingBackboneError,
    NonStandardResidueError,
    ResidueNotFoundError,
)


class TestMutationValidator:
    """MutationValidator 测试类。"""

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

    @pytest.fixture
    def parser_with_nonstandard(self, tmp_path: Path) -> PdbParser:
        """创建包含非标准残基的解析器。"""
        content = """\
ATOM      1  N   ALA A   1       1.000   2.000   3.000  1.00  0.00           N
ATOM      2  CA  ALA A   1       2.000   2.000   3.000  1.00  0.00           C
ATOM      3  C   ALA A   1       3.000   2.000   3.000  1.00  0.00           C
ATOM      4  O   ALA A   1       4.000   2.000   3.000  1.00  0.00           O
ATOM      5  CB  ALA A   1       1.500   1.000   3.000  1.00  0.00           C
ATOM      6  N   MSE A   2       5.000   2.000   3.000  1.00  0.00           N
ATOM      7  CA  MSE A   2       6.000   2.000   3.000  1.00  0.00           C
ATOM      8  C   MSE A   2       7.000   2.000   3.000  1.00  0.00           C
ATOM      9  O   MSE A   2       8.000   2.000   3.000  1.00  0.00           O
END
"""
        pdb_path = tmp_path / "nonstandard.pdb"
        pdb_path.write_text(content, encoding="utf-8")
        parser = PdbParser()
        parser.parse(pdb_path)
        return parser

    @pytest.fixture
    def parser_with_missing_backbone(self, tmp_path: Path) -> PdbParser:
        """创建主链原子缺失的解析器。"""
        content = """\
ATOM      1  N   ALA A   1       1.000   2.000   3.000  1.00  0.00           N
ATOM      2  CA  ALA A   1       2.000   2.000   3.000  1.00  0.00           C
ATOM      3  CB  ALA A   1       1.500   1.000   3.000  1.00  0.00           C
END
"""
        pdb_path = tmp_path / "missing_bb.pdb"
        pdb_path.write_text(content, encoding="utf-8")
        parser = PdbParser()
        parser.parse(pdb_path)
        return parser

    def test_validate_valid_target(self, parser_with_ala):
        """测试验证有效的目标残基。"""
        validator = MutationValidator(parser_with_ala)
        residue = validator.validate_target(1, "A", 1)

        assert residue is not None
        assert residue.residue_name == "ALA"
        assert residue.residue_number == 1

    def test_validate_nonexistent_residue(self, parser_with_ala):
        """测试验证不存在的残基。"""
        validator = MutationValidator(parser_with_ala)
        with pytest.raises(ResidueNotFoundError):
            validator.validate_target(1, "A", 999)

    def test_validate_nonexistent_chain(self, parser_with_ala):
        """测试验证不存在的 chain。"""
        validator = MutationValidator(parser_with_ala)
        with pytest.raises(ResidueNotFoundError):
            validator.validate_target(1, "Z", 1)

    def test_validate_nonstandard_residue(self, parser_with_nonstandard):
        """测试验证非标准残基。"""
        validator = MutationValidator(parser_with_nonstandard)
        with pytest.raises(NonStandardResidueError):
            validator.validate_target(1, "A", 2)

    def test_validate_missing_backbone(self, parser_with_missing_backbone):
        """测试验证主链原子缺失。"""
        validator = MutationValidator(parser_with_missing_backbone)
        with pytest.raises(MissingBackboneError):
            validator.validate_target(1, "A", 1)

    def test_validate_and_check_valid(self, parser_with_ala):
        """测试综合验证。"""
        validator = MutationValidator(parser_with_ala)
        residue, clashes = validator.validate_and_check(
            1, "A", 1, "GLY"
        )
        assert residue.residue_name == "ALA"
        assert isinstance(clashes, list)

    def test_validate_and_check_invalid_new_residue(self, parser_with_ala):
        """测试无效的新氨基酸。"""
        validator = MutationValidator(parser_with_ala)
        with pytest.raises(NonStandardResidueError):
            validator.validate_and_check(1, "A", 1, "XYZ")

    def test_clash_detection(self, tmp_path: Path):
        """测试碰撞检测。"""
        # 创建两个距离很近的残基
        content = """\
ATOM      1  N   ALA A   1       1.000   2.000   3.000  1.00  0.00           N
ATOM      2  CA  ALA A   1       2.000   2.000   3.000  1.00  0.00           C
ATOM      3  C   ALA A   1       3.000   2.000   3.000  1.00  0.00           C
ATOM      4  O   ALA A   1       4.000   2.000   3.000  1.00  0.00           O
ATOM      5  CB  ALA A   1       2.500   2.000   3.000  1.00  0.00           C
ATOM      6  N   GLY A   2       5.000   2.000   3.000  1.00  0.00           N
ATOM      7  CA  GLY A   2       6.000   2.000   3.000  1.00  0.00           C
ATOM      8  C   GLY A   2       7.000   2.000   3.000  1.00  0.00           C
ATOM      9  O   GLY A   2       8.000   2.000   3.000  1.00  0.00           O
END
"""
        pdb_path = tmp_path / "clash_test.pdb"
        pdb_path.write_text(content, encoding="utf-8")
        parser = PdbParser()
        parser.parse(pdb_path)

        validator = MutationValidator(parser)
        target = parser.get_residue(1, "A", 1)
        assert target is not None

        clashes = validator.check_clashes(target, "GLY", threshold=5.0)
        assert len(clashes) > 0
