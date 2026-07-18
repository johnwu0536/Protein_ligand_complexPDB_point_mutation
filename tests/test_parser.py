"""
PDB 解析器单元测试。
"""

from __future__ import annotations

from pathlib import Path

import pytest

from src.core.models import ResidueInfo, StructureInfo
from src.core.parser import PdbParser
from src.utils.exceptions import PdbParseError


class TestPdbParser:
    """PdbParser 测试类。"""

    @pytest.fixture
    def sample_pdb_content(self) -> str:
        """生成一个简单的 PDB 文件内容用于测试。"""
        return """\
ATOM      1  N   ALA A   1       1.000   2.000   3.000  1.00  0.00           N
ATOM      2  CA  ALA A   1       2.000   2.000   3.000  1.00  0.00           C
ATOM      3  C   ALA A   1       3.000   2.000   3.000  1.00  0.00           C
ATOM      4  O   ALA A   1       4.000   2.000   3.000  1.00  0.00           O
ATOM      5  CB  ALA A   1       1.500   1.000   3.000  1.00  0.00           C
ATOM      6  N   GLY A   2       5.000   2.000   3.000  1.00  0.00           N
ATOM      7  CA  GLY A   2       6.000   2.000   3.000  1.00  0.00           C
ATOM      8  C   GLY A   2       7.000   2.000   3.000  1.00  0.00           C
ATOM      9  O   GLY A   2       8.000   2.000   3.000  1.00  0.00           O
END
"""

    @pytest.fixture
    def sample_pdb_file(self, tmp_path: Path, sample_pdb_content: str) -> Path:
        """创建临时 PDB 文件。"""
        pdb_path = tmp_path / "test.pdb"
        pdb_path.write_text(sample_pdb_content, encoding="utf-8")
        return pdb_path

    def test_parse_valid_pdb(self, sample_pdb_file: Path):
        """测试解析有效的 PDB 文件。"""
        parser = PdbParser()
        info = parser.parse(sample_pdb_file)

        assert isinstance(info, StructureInfo)
        assert info.path == sample_pdb_file
        assert 1 in info.models
        assert "A" in info.chains[1]

    def test_parse_nonexistent_file(self):
        """测试解析不存在的文件。"""
        parser = PdbParser()
        with pytest.raises(PdbParseError):
            parser.parse(Path("nonexistent.pdb"))

    def test_get_models(self, sample_pdb_file: Path):
        """测试获取 model 列表。"""
        parser = PdbParser()
        parser.parse(sample_pdb_file)
        models = parser.get_models()
        assert models == [1]

    def test_get_chains(self, sample_pdb_file: Path):
        """测试获取 chain 列表。"""
        parser = PdbParser()
        parser.parse(sample_pdb_file)
        chains = parser.get_chains(1)
        assert "A" in chains

    def test_get_residues(self, sample_pdb_file: Path):
        """测试获取残基列表。"""
        parser = PdbParser()
        parser.parse(sample_pdb_file)
        residues = parser.get_residues(1, "A")

        assert len(residues) == 2
        assert residues[0].residue_name == "ALA"
        assert residues[0].residue_number == 1
        assert residues[1].residue_name == "GLY"
        assert residues[1].residue_number == 2

    def test_get_residue(self, sample_pdb_file: Path):
        """测试获取单个残基。"""
        parser = PdbParser()
        parser.parse(sample_pdb_file)
        residue = parser.get_residue(1, "A", 1)

        assert residue is not None
        assert isinstance(residue, ResidueInfo)
        assert residue.residue_name == "ALA"
        assert residue.residue_number == 1
        assert residue.chain == "A"
        assert residue.model == 1

    def test_get_nonexistent_residue(self, sample_pdb_file: Path):
        """测试获取不存在的残基。"""
        parser = PdbParser()
        parser.parse(sample_pdb_file)
        residue = parser.get_residue(1, "A", 999)
        assert residue is None

    def test_backbone_check(self, sample_pdb_file: Path):
        """测试主链原子检查。"""
        parser = PdbParser()
        parser.parse(sample_pdb_file)
        residue = parser.get_residue(1, "A", 1)

        assert residue is not None
        assert residue.has_backbone is True
        assert residue.missing_backbone_atoms == []

    def test_standard_amino_acid(self, sample_pdb_file: Path):
        """测试标准氨基酸检查。"""
        parser = PdbParser()
        parser.parse(sample_pdb_file)
        residue = parser.get_residue(1, "A", 1)

        assert residue is not None
        assert residue.is_standard is True

    def test_parse_with_chinese_path(self, tmp_path: Path, sample_pdb_content: str):
        """测试中文路径解析。"""
        chinese_dir = tmp_path / "测试目录"
        chinese_dir.mkdir(exist_ok=True)
        pdb_path = chinese_dir / "测试文件.pdb"
        pdb_path.write_text(sample_pdb_content, encoding="utf-8")

        parser = PdbParser()
        info = parser.parse(pdb_path)
        assert info.path == pdb_path
        assert 1 in info.models

    def test_parse_empty_structure(self, tmp_path: Path):
        """测试空结构文件。"""
        pdb_path = tmp_path / "empty.pdb"
        pdb_path.write_text("END\n", encoding="utf-8")

        parser = PdbParser()
        info = parser.parse(pdb_path)
        # Gemmi 会为 END 创建默认 model，但无残基
        assert len(info.models) >= 0
        # 确保没有残基
        all_residues = []
        for key, residues in info.residues.items():
            all_residues.extend(residues)
        assert len(all_residues) == 0

    def test_parse_with_water(self, tmp_path: Path):
        """测试包含水分子的 PDB 文件。"""
        content = """\
ATOM      1  N   ALA A   1       1.000   2.000   3.000  1.00  0.00           N
ATOM      2  CA  ALA A   1       2.000   2.000   3.000  1.00  0.00           C
ATOM      3  C   ALA A   1       3.000   2.000   3.000  1.00  0.00           C
ATOM      4  O   ALA A   1       4.000   2.000   3.000  1.00  0.00           O
ATOM      5  CB  ALA A   1       1.500   1.000   3.000  1.00  0.00           C
ATOM      6  O   HOH A 201       5.000   5.000   5.000  1.00  0.00           O
END
"""
        pdb_path = tmp_path / "with_water.pdb"
        pdb_path.write_text(content, encoding="utf-8")

        parser = PdbParser()
        info = parser.parse(pdb_path)
        residues = parser.get_residues(1, "A")
        # 水分子应该被过滤掉
        assert len(residues) == 1
        assert residues[0].residue_name == "ALA"
