"""
PDBMutator 自定义异常模块。
所有异常消息均为中文，确保用户可理解。
"""

from pathlib import Path
from typing import Optional


class PDBMutatorError(Exception):
    """所有 PDBMutator 异常的基类。"""

    def __init__(self, message: str, detail: Optional[str] = None) -> None:
        self.message = message
        self.detail = detail
        super().__init__(self._format_message())

    def _format_message(self) -> str:
        if self.detail:
            return f"{self.message}\n详细信息: {self.detail}"
        return self.message


class PdbParseError(PDBMutatorError):
    """PDB 文件解析失败。"""

    def __init__(self, path: Path, detail: Optional[str] = None) -> None:
        message = f"无法解析PDB文件: {path}"
        super().__init__(message, detail)


class ResidueNotFoundError(PDBMutatorError):
    """目标残基不存在。"""

    def __init__(
        self,
        model: int,
        chain: str,
        residue_number: int,
        insertion_code: str = "",
        detail: Optional[str] = None,
    ) -> None:
        ic = f", insertion_code='{insertion_code}'" if insertion_code else ""
        message = f"未找到目标残基: model={model}, chain='{chain}', residue_number={residue_number}{ic}"
        super().__init__(message, detail)


class NonStandardResidueError(PDBMutatorError):
    """非标准氨基酸。"""

    def __init__(self, residue_name: str, detail: Optional[str] = None) -> None:
        message = f"残基 '{residue_name}' 不是20种标准氨基酸之一，无法执行突变"
        super().__init__(message, detail)


class MissingBackboneError(PDBMutatorError):
    """主链原子缺失。"""

    def __init__(
        self,
        residue_name: str,
        residue_number: int,
        missing_atoms: list[str],
        detail: Optional[str] = None,
    ) -> None:
        atoms_str = ", ".join(missing_atoms)
        message = (
            f"残基 {residue_name}({residue_number}) 缺失主链原子: {atoms_str}。"
            f"主链原子 (N, CA, C, O) 必须完整才能执行突变"
        )
        super().__init__(message, detail)


class AltLocPresentError(PDBMutatorError):
    """存在 altLoc 需要处理。"""

    def __init__(
        self,
        residue_name: str,
        residue_number: int,
        alt_locs: list[str],
        detail: Optional[str] = None,
    ) -> None:
        locs_str = ", ".join(alt_locs)
        message = (
            f"残基 {residue_name}({residue_number}) 存在交替位置 (altLoc): {locs_str}。"
            f"请先处理 altLoc 后再执行突变"
        )
        super().__init__(message, detail)


class DisulfideBondError(PDBMutatorError):
    """参与二硫键。"""

    def __init__(
        self,
        residue_name: str,
        residue_number: int,
        partner: str,
        detail: Optional[str] = None,
    ) -> None:
        message = (
            f"残基 {residue_name}({residue_number}) 参与二硫键 (与 {partner})。"
            f"突变半胱氨酸会破坏二硫键，请确认操作"
        )
        super().__init__(message, detail)


class MutationError(PDBMutatorError):
    """突变执行失败。"""

    def __init__(self, detail: Optional[str] = None) -> None:
        message = "执行突变操作时发生错误"
        super().__init__(message, detail)


class ClashError(PDBMutatorError):
    """原子碰撞冲突。"""

    def __init__(
        self,
        clashes: list[tuple[str, str, float]],
        threshold: float = 2.0,
        detail: Optional[str] = None,
    ) -> None:
        clash_msgs = []
        for atom1, atom2, dist in clashes[:5]:  # 最多显示5个冲突
            clash_msgs.append(f"  {atom1} <-> {atom2}: {dist:.2f} Å")
        if len(clashes) > 5:
            clash_msgs.append(f"  ... 还有 {len(clashes) - 5} 个冲突")
        clash_str = "\n".join(clash_msgs)
        message = (
            f"检测到 {len(clashes)} 个原子碰撞 (阈值: {threshold} Å):\n{clash_str}\n"
            f"建议选择其他突变位置或氨基酸"
        )
        super().__init__(message, detail)


class InvalidInputError(PDBMutatorError):
    """输入参数无效。"""

    def __init__(self, field: str, detail: Optional[str] = None) -> None:
        message = f"输入参数 '{field}' 无效"
        super().__init__(message, detail)
