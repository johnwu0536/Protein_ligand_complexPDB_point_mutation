"""
PDBMutator 程序入口。
"""

from __future__ import annotations

import logging
import sys
from pathlib import Path

from PySide6.QtWidgets import QApplication

from src.gui.main_window import MainWindow


def setup_logging() -> None:
    """配置日志系统。"""
    log_dir = Path.cwd() / "logs"
    log_dir.mkdir(parents=True, exist_ok=True)

    log_file = log_dir / "pdbutator.log"

    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
        handlers=[
            logging.FileHandler(str(log_file), encoding="utf-8"),
            logging.StreamHandler(sys.stdout),
        ],
    )


def main() -> None:
    """程序主入口。"""
    setup_logging()
    logger = logging.getLogger(__name__)
    logger.info("PDBMutator 启动")

    app = QApplication(sys.argv)
    app.setApplicationName("PDBMutator")
    app.setOrganizationName("PDBMutator")

    window = MainWindow()
    window.show()

    try:
        exit_code = app.exec()
        logger.info(f"PDBMutator 退出，退出码: {exit_code}")
        sys.exit(exit_code)
    except Exception as e:
        logger.critical(f"程序异常退出: {e}", exc_info=True)
        sys.exit(1)


if __name__ == "__main__":
    main()
