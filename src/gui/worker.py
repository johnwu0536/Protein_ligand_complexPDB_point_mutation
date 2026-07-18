"""
PDBMutator 后台工作线程模块。
使用 QThread 在后台执行突变操作，避免 GUI 卡死。
"""

from __future__ import annotations

import logging
from datetime import datetime
from pathlib import Path
from typing import Optional

from PySide6.QtCore import QThread, Signal

from src.core.engine import BaseMutationEngine, PdbFixerMutationEngine
from src.core.logger import MutationLogger
from src.core.models import MutationRecord, MutationResult, ResidueInfo
from src.core.parser import PdbParser
from src.core.validator import MutationValidator

logger = logging.getLogger(__name__)


class MutationWorker(QThread):
    """后台突变工作线程。"""

    # 信号定义
    started = Signal(str)  # 开始消息
    progress = Signal(str)  # 进度消息
    finished = Signal(object)  # MutationResult
    error = Signal(str)  # 错误消息（中文）
    warning = Signal(str)  # 警告消息（中文）

    def __init__(
        self,
        parser: PdbParser,
        target: ResidueInfo,
        new_residue_name: str,
        output_path: Path,
        engine: Optional[BaseMutationEngine] = None,
        logger_instance: Optional[MutationLogger] = None,
        parent=None,
    ) -> None:
        """
        初始化工作线程。

        Args:
            parser: PDB 解析器
            target: 目标残基
            new_residue_name: 新氨基酸
            output_path: 输出路径
            engine: 突变引擎（默认使用 PDBFixer）
            logger_instance: 日志记录器
            parent: 父对象
        """
        super().__init__(parent)
        self._parser = parser
        self._target = target
        self._new_residue_name = new_residue_name
        self._output_path = output_path
        self._engine = engine or PdbFixerMutationEngine()
        self._logger = logger_instance or MutationLogger()

    def run(self) -> None:
        """执行后台任务。"""
        try:
            self.started.emit("开始执行突变...")

            # 1. 验证目标残基
            self.progress.emit("正在验证目标残基...")
            validator = MutationValidator(self._parser)
            warnings: list[str] = []

            try:
                residue = validator.validate_target(
                    self._target.model,
                    self._target.chain,
                    self._target.residue_number,
                    self._target.insertion_code,
                )
            except Exception as e:
                self.error.emit(str(e))
                return

            # 2. 检查碰撞
            self.progress.emit("正在检查原子碰撞...")
            clashes = validator.check_clashes(
                residue, self._new_residue_name
            )

            if clashes:
                clash_msg = (
                    f"检测到 {len(clashes)} 个原子碰撞 "
                    f"(阈值: {MutationValidator.DEFAULT_CLASH_THRESHOLD} Å)"
                )
                self.warning.emit(clash_msg)
                warnings.append(clash_msg)

            # 3. 执行突变
            self.progress.emit(
                f"正在执行突变: {residue.residue_name}{residue.residue_number} "
                f"-> {self._new_residue_name}..."
            )

            result = self._engine.mutate(
                self._parser,
                residue,
                self._new_residue_name,
                self._output_path,
            )

            # 4. 记录日志
            if result.success:
                record = MutationRecord(
                    input_file=str(self._parser.path),
                    output_file=str(self._output_path),
                    model=self._target.model,
                    chain=self._target.chain,
                    residue_number=self._target.residue_number,
                    insertion_code=self._target.insertion_code,
                    original_residue=residue.residue_name,
                    new_residue=self._new_residue_name,
                    timestamp=datetime.now().isoformat(),
                    warnings=warnings + result.warnings,
                    clashes=[
                        {"atom1": c[0], "atom2": c[1], "distance": round(c[2], 3)}
                        for c in clashes
                    ],
                    success=True,
                )
                try:
                    log_path = self._logger.log(record)
                    self.progress.emit(f"日志已保存: {log_path}")
                except Exception as e:
                    self.warning.emit(f"保存日志失败: {e}")

                self.progress.emit("突变完成！")
            else:
                # 记录失败日志
                record = MutationRecord(
                    input_file=str(self._parser.path),
                    output_file=str(self._output_path),
                    model=self._target.model,
                    chain=self._target.chain,
                    residue_number=self._target.residue_number,
                    insertion_code=self._target.insertion_code,
                    original_residue=residue.residue_name,
                    new_residue=self._new_residue_name,
                    timestamp=datetime.now().isoformat(),
                    warnings=warnings + result.warnings,
                    success=False,
                )
                try:
                    self._logger.log(record)
                except Exception:
                    pass

                self.error.emit(result.error_message or "未知错误")
                return

            self.finished.emit(result)

        except Exception as e:
            logger.error("后台工作线程异常", exc_info=True)
            self.error.emit(f"执行突变时发生未知错误: {str(e)}")
