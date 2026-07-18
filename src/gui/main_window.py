"""
PDBMutator 主窗口模块。
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Optional

from PySide6.QtWidgets import (
    QApplication,
    QComboBox,
    QFileDialog,
    QGroupBox,
    QHBoxLayout,
    QLabel,
    QMainWindow,
    QMessageBox,
    QPushButton,
    QSplitter,
    QVBoxLayout,
    QWidget,
)
from PySide6.QtCore import Qt, QTimer

from src.core.engine import PdbFixerMutationEngine
from src.core.logger import MutationLogger
from src.core.models import MutationResult, ResidueInfo, STANDARD_AMINO_ACIDS
from src.core.parser import PdbParser
from src.gui.widgets import AminoAcidComboBox, LogWindow, ResidueTableWidget
from src.gui.worker import MutationWorker

logger = logging.getLogger(__name__)


class MainWindow(QMainWindow):
    """PDBMutator 主窗口。"""

    def __init__(self) -> None:
        super().__init__()
        self._parser: Optional[PdbParser] = None
        self._current_residue: Optional[ResidueInfo] = None
        self._worker: Optional[MutationWorker] = None
        self._logger = MutationLogger()

        self._setup_ui()
        self._connect_signals()

    def _setup_ui(self) -> None:
        """初始化用户界面。"""
        self.setWindowTitle("PDBMutator - 蛋白质点突变工具")
        self.setMinimumSize(900, 700)

        # 中央部件
        central = QWidget()
        self.setCentralWidget(central)
        main_layout = QVBoxLayout(central)

        # 文件操作区域
        file_group = QGroupBox("文件操作")
        file_layout = QHBoxLayout(file_group)

        self.open_button = QPushButton("打开PDB文件")
        self.save_button = QPushButton("保存PDB文件")
        self.save_button.setEnabled(False)

        file_layout.addWidget(self.open_button)
        file_layout.addWidget(self.save_button)
        file_layout.addStretch()

        main_layout.addWidget(file_group)

        # 选择区域
        select_group = QGroupBox("选择参数")
        select_layout = QHBoxLayout(select_group)

        select_layout.addWidget(QLabel("Model:"))
        self.model_combo = QComboBox()
        select_layout.addWidget(self.model_combo)

        select_layout.addWidget(QLabel("Chain:"))
        self.chain_combo = QComboBox()
        select_layout.addWidget(self.chain_combo)

        select_layout.addWidget(QLabel("突变至:"))
        self.amino_acid_combo = AminoAcidComboBox()
        select_layout.addWidget(self.amino_acid_combo)

        self.mutate_button = QPushButton("执行突变")
        self.mutate_button.setEnabled(False)
        select_layout.addWidget(self.mutate_button)

        main_layout.addWidget(select_group)

        # 分割器：残基表格 + 日志
        splitter = QSplitter(Qt.Vertical)

        # 残基表格
        self.residue_table = ResidueTableWidget()
        splitter.addWidget(self.residue_table)

        # 日志窗口
        self.log_window = LogWindow()
        self.log_window.setMinimumHeight(150)
        splitter.addWidget(self.log_window)

        main_layout.addWidget(splitter)

        # 状态栏
        self.statusBar().showMessage("就绪 - 请打开PDB文件")

    def _connect_signals(self) -> None:
        """连接信号与槽。"""
        self.open_button.clicked.connect(self._on_open_pdb)
        self.save_button.clicked.connect(self._on_save_pdb)
        self.model_combo.currentIndexChanged.connect(self._on_model_changed)
        self.chain_combo.currentIndexChanged.connect(self._on_chain_changed)
        self.mutate_button.clicked.connect(self._on_mutate)
        self.residue_table.residue_selected.connect(self._on_residue_selected)

    def _on_open_pdb(self) -> None:
        """打开 PDB 文件。"""
        file_path, _ = QFileDialog.getOpenFileName(
            self,
            "选择PDB文件",
            "",
            "PDB文件 (*.pdb *.ent);;所有文件 (*.*)",
        )
        if not file_path:
            return

        path = Path(file_path)
        self.log_window.append_info(f"正在打开文件: {path}")

        try:
            self._parser = PdbParser()
            info = self._parser.parse(path)
            self.log_window.append_success(f"成功解析PDB文件: {path.name}")
            self.statusBar().showMessage(f"已加载: {path.name}")

            # 更新 model 选择框
            self.model_combo.blockSignals(True)
            self.model_combo.clear()
            for m in info.models:
                self.model_combo.addItem(f"Model {m}", m)
            self.model_combo.blockSignals(False)

            # 自动选择第一个 model
            if info.models:
                self.model_combo.setCurrentIndex(0)
                self._on_model_changed()

            self.save_button.setEnabled(True)

        except Exception as e:
            self.log_window.append_error(f"打开文件失败: {str(e)}")
            QMessageBox.critical(self, "错误", f"无法打开PDB文件:\n{str(e)}")

    def _on_model_changed(self) -> None:
        """Model 选择变化时更新 Chain 列表。"""
        if self._parser is None:
            return

        model_id = self.model_combo.currentData()
        if model_id is None:
            return

        chains = self._parser.get_chains(model_id)
        self.chain_combo.blockSignals(True)
        self.chain_combo.clear()
        for c in chains:
            self.chain_combo.addItem(f"Chain {c}", c)
        self.chain_combo.blockSignals(False)

        if chains:
            self.chain_combo.setCurrentIndex(0)
            self._on_chain_changed()

    def _on_chain_changed(self) -> None:
        """Chain 选择变化时更新残基表格。"""
        if self._parser is None:
            return

        model_id = self.model_combo.currentData()
        chain_id = self.chain_combo.currentData()

        if model_id is None or chain_id is None:
            return

        residues = self._parser.get_residues(model_id, chain_id)
        self.residue_table.populate(residues)
        self.log_window.append_info(
            f"Chain {chain_id}: 共 {len(residues)} 个残基"
        )

    def _on_residue_selected(self, data: dict) -> None:
        """残基选中时更新状态。"""
        if self._parser is None:
            return

        model_id = self.model_combo.currentData()
        chain_id = self.chain_combo.currentData()

        self._current_residue = self._parser.get_residue(
            model_id,
            chain_id,
            int(data["residue_number"]),
            data["insertion_code"],
        )

        if self._current_residue:
            self.mutate_button.setEnabled(True)
            self.statusBar().showMessage(
                f"已选择: {self._current_residue.display_name()} "
                f"(残基 {self._current_residue.residue_number})"
            )
        else:
            self.mutate_button.setEnabled(False)

    def _on_mutate(self) -> None:
        """执行突变操作。"""
        if self._parser is None or self._current_residue is None:
            return

        new_aa = self.amino_acid_combo.current_amino_acid()
        if not new_aa:
            QMessageBox.warning(self, "警告", "请选择目标氨基酸")
            return

        # 检查是否选择了相同的氨基酸
        if new_aa == self._current_residue.residue_name:
            reply = QMessageBox.question(
                self,
                "确认",
                f"目标残基已经是 {self._current_residue.display_name()}，"
                f"确定要突变为相同的氨基酸吗？",
                QMessageBox.Yes | QMessageBox.No,
            )
            if reply == QMessageBox.No:
                return

        # 选择输出文件路径
        output_path, _ = QFileDialog.getSaveFileName(
            self,
            "保存突变后的PDB文件",
            str(self._parser.path.stem) + f"_{new_aa}.pdb",
            "PDB文件 (*.pdb);;所有文件 (*.*)",
        )
        if not output_path:
            return

        # 禁用按钮，防止重复操作
        self._set_ui_enabled(False)
        self.log_window.append_info(
            f"开始突变: {self._current_residue.display_name()} -> "
            f"{new_aa}"
        )

        # 启动后台工作线程
        self._worker = MutationWorker(
            parser=self._parser,
            target=self._current_residue,
            new_residue_name=new_aa,
            output_path=Path(output_path),
            logger_instance=self._logger,
            parent=self,
        )

        # 连接信号
        self._worker.started.connect(self.log_window.append_info)
        self._worker.progress.connect(self.log_window.append_info)
        self._worker.warning.connect(self.log_window.append_warning)
        self._worker.error.connect(self._on_mutation_error)
        self._worker.finished.connect(self._on_mutation_finished)

        self._worker.start()

    def _on_mutation_finished(self, result: MutationResult) -> None:
        """突变完成回调。"""
        self._set_ui_enabled(True)

        if result.success and result.output_path:
            self.log_window.append_success(
                f"突变成功！输出文件: {result.output_path}"
            )
            self.statusBar().showMessage(
                f"突变完成: {result.output_path.name}"
            )

            # 显示警告
            for w in result.warnings:
                self.log_window.append_warning(w)

            QMessageBox.information(
                self,
                "成功",
                f"突变成功完成！\n\n输出文件: {result.output_path}",
            )
        else:
            self.log_window.append_error(
                result.error_message or "突变失败"
            )

    def _on_mutation_error(self, error_msg: str) -> None:
        """突变错误回调。"""
        self._set_ui_enabled(True)
        self.log_window.append_error(error_msg)
        QMessageBox.critical(self, "突变失败", error_msg)

    def _on_save_pdb(self) -> None:
        """保存当前 PDB 结构。"""
        if self._parser is None:
            return

        output_path, _ = QFileDialog.getSaveFileName(
            self,
            "保存PDB文件",
            "",
            "PDB文件 (*.pdb);;所有文件 (*.*)",
        )
        if not output_path:
            return

        try:
            self._parser.structure.write_pdb(str(output_path))
            self.log_window.append_success(f"文件已保存: {output_path}")
            self.statusBar().showMessage(f"已保存: {Path(output_path).name}")
        except Exception as e:
            self.log_window.append_error(f"保存失败: {str(e)}")
            QMessageBox.critical(self, "错误", f"保存文件失败:\n{str(e)}")

    def _set_ui_enabled(self, enabled: bool) -> None:
        """设置 UI 组件的启用状态。"""
        self.open_button.setEnabled(enabled)
        self.save_button.setEnabled(enabled)
        self.model_combo.setEnabled(enabled)
        self.chain_combo.setEnabled(enabled)
        self.amino_acid_combo.setEnabled(enabled)
        self.mutate_button.setEnabled(enabled and self._current_residue is not None)
        self.residue_table.setEnabled(enabled)

    def closeEvent(self, event) -> None:
        """窗口关闭事件。"""
        if self._worker and self._worker.isRunning():
            reply = QMessageBox.question(
                self,
                "确认退出",
                "突变操作正在执行中，确定要退出吗？",
                QMessageBox.Yes | QMessageBox.No,
            )
            if reply == QMessageBox.No:
                event.ignore()
                return
            self._worker.quit()
            self._worker.wait(3000)

        event.accept()
