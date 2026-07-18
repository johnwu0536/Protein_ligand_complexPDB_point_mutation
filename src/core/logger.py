"""
PDBMutator 日志模块。
负责将突变记录保存为 JSON 格式的日志文件。
"""

from __future__ import annotations

import json
import logging
from datetime import datetime
from pathlib import Path
from typing import Optional

from src.core.models import MutationRecord

logger = logging.getLogger(__name__)


class MutationLogger:
    """突变日志记录器，将突变记录保存为 JSON 文件。"""

    def __init__(self, log_dir: Optional[Path | str] = None) -> None:
        """
        初始化日志记录器。

        Args:
            log_dir: 日志文件保存目录。默认为当前目录下的 logs 文件夹。
        """
        if log_dir is None:
            self._log_dir = Path.cwd() / "logs"
        else:
            self._log_dir = Path(log_dir)
        self._log_dir.mkdir(parents=True, exist_ok=True)

    @property
    def log_dir(self) -> Path:
        """获取日志目录。"""
        return self._log_dir

    def log(self, record: MutationRecord) -> Path:
        """
        记录一条突变日志。

        Args:
            record: 突变记录

        Returns:
            Path: 日志文件路径
        """
        # 生成日志文件名
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        log_filename = f"mutation_{timestamp}.json"
        log_path = self._log_dir / log_filename

        # 转换为可序列化的字典
        log_data = {
            "input_file": str(record.input_file),
            "output_file": str(record.output_file),
            "model": record.model,
            "chain": record.chain,
            "residue_number": record.residue_number,
            "insertion_code": record.insertion_code,
            "original_residue": record.original_residue,
            "new_residue": record.new_residue,
            "timestamp": record.timestamp,
            "warnings": record.warnings,
            "clashes": record.clashes,
            "success": record.success,
        }

        # 写入 JSON 文件
        try:
            with open(str(log_path), "w", encoding="utf-8") as f:
                json.dump(log_data, f, ensure_ascii=False, indent=2)
            logger.info(f"突变日志已保存: {log_path}")
        except Exception as e:
            logger.error(f"保存突变日志失败: {e}", exc_info=True)
            raise

        return log_path

    def get_all_logs(self) -> list[Path]:
        """获取所有日志文件路径。"""
        return sorted(self._log_dir.glob("mutation_*.json"))

    def read_log(self, log_path: Path) -> Optional[dict]:
        """读取并解析日志文件。"""
        try:
            with open(str(log_path), "r", encoding="utf-8") as f:
                return json.load(f)
        except Exception as e:
            logger.error(f"读取日志文件失败: {log_path}, {e}", exc_info=True)
            return None
