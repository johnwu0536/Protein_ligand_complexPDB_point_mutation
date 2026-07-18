"""
PDBMutator 自定义 GUI 组件模块。
"""

from __future__ import annotations

from PySide6.QtWidgets import (
    QComboBox,
    QTableWidget,
    QTableWidgetItem,
    QTextEdit,
    QHeaderView,
)
from PySide6.QtCore import Qt, Signal
from PySide6.QtGui import QColor, QFont

from src.core.models import STANDARD_AMINO_ACIDS, AMINO_ACID_NAMES_CN


class AminoAcidComboBox(QComboBox):
    """氨基酸选择下拉框。"""

    def __init__(self, parent=None) -> None:
        super().__init__(parent)
        self._populate()

    def _populate(self) -> None:
        """填充20种标准氨基酸。"""
        self.clear()
        for aa in STANDARD_AMINO_ACIDS:
            cn_name = AMINO_ACID_NAMES_CN.get(aa, "")
            self.addItem(f"{aa} ({cn_name})", aa)

    def current_amino_acid(self) -> str:
        """获取当前选中的氨基酸三字母代码。"""
        return self.currentData()


class ResidueTableWidget(QTableWidget):
    """残基信息表格组件。"""

    # 列索引
    COL_RES_NUM = 0
    COL_INS_CODE = 1
    COL_RES_NAME = 2
    COL_BACKBONE = 3
    COL_ALT_LOC = 4
    COL_DISULFIDE = 5

    # 信号：选中残基变化
    residue_selected = Signal(dict)

    def __init__(self, parent=None) -> None:
        super().__init__(parent)
        self._setup_ui()

    def _setup_ui(self) -> None:
        """初始化表格界面。"""
        headers = [
            "残基编号", "插入码", "残基名称",
            "主链完整", "AltLoc", "二硫键"
        ]
        self.setColumnCount(len(headers))
        self.setHorizontalHeaderLabels(headers)

        # 设置表头样式
        header = self.horizontalHeader()
        header.setSectionResizeMode(QHeaderView.ResizeToContents)
        header.setStretchLastSection(True)

        # 设置选择模式
        self.setSelectionBehavior(QTableWidget.SelectRows)
        self.setSelectionMode(QTableWidget.SingleSelection)

        # 连接点击事件
        self.itemSelectionChanged.connect(self._on_selection_changed)

    def _on_selection_changed(self) -> None:
        """选中行变化时发射信号。"""
        rows = self.selectedIndexes()
        if not rows:
            return

        row = rows[0].row()
        data = {
            "residue_number": self._get_item_data(row, self.COL_RES_NUM),
            "insertion_code": self._get_item_data(row, self.COL_INS_CODE),
            "residue_name": self._get_item_data(row, self.COL_RES_NAME),
            "has_backbone": self._get_item_data(row, self.COL_BACKBONE) == "是",
            "alt_loc": self._get_item_data(row, self.COL_ALT_LOC),
            "disulfide": self._get_item_data(row, self.COL_DISULFIDE),
        }
        self.residue_selected.emit(data)

    def _get_item_data(self, row: int, col: int) -> str:
        """获取表格单元格数据。"""
        item = self.item(row, col)
        return item.text() if item else ""

    def populate(self, residues: list) -> None:
        """
        填充残基数据。

        Args:
            residues: ResidueInfo 列表
        """
        self.setRowCount(len(residues))

        for i, res in enumerate(residues):
            # 残基编号
            self._set_item(i, self.COL_RES_NUM, str(res.residue_number))

            # 插入码
            self._set_item(i, self.COL_INS_CODE, res.insertion_code or "")

            # 残基名称
            self._set_item(i, self.COL_RES_NAME, res.display_name())

            # 主链完整
            bb_status = "是" if res.has_backbone else "否"
            bb_item = self._set_item(i, self.COL_BACKBONE, bb_status)
            if not res.has_backbone:
                bb_item.setBackground(QColor(255, 200, 200))  # 浅红色

            # AltLoc
            alt_loc = res.alt_loc if res.alt_loc else ""
            alt_item = self._set_item(i, self.COL_ALT_LOC, alt_loc)
            if res.alt_loc:
                alt_item.setBackground(QColor(255, 255, 200))  # 浅黄色

            # 二硫键（暂不显示，需要额外信息）
            self._set_item(i, self.COL_DISULFIDE, "")

    def _set_item(self, row: int, col: int, text: str) -> QTableWidgetItem:
        """设置表格单元格。"""
        item = QTableWidgetItem(text)
        item.setFlags(item.flags() & ~Qt.ItemIsEditable)  # 只读
        self.setItem(row, col, item)
        return item


class LogWindow(QTextEdit):
    """日志显示窗口。"""

    def __init__(self, parent=None) -> None:
        super().__init__(parent)
        self.setReadOnly(True)
        self.setFont(QFont("Consolas", 9))

        # 颜色定义
        self._info_color = QColor(0, 0, 0)       # 黑色
        self._warning_color = QColor(200, 150, 0)  # 橙色
        self._error_color = QColor(200, 0, 0)     # 红色
        self._success_color = QColor(0, 150, 0)   # 绿色

    def append_info(self, message: str) -> None:
        """添加信息日志。"""
        self._append_colored(f"[信息] {message}", self._info_color)

    def append_warning(self, message: str) -> None:
        """添加警告日志。"""
        self._append_colored(f"[警告] {message}", self._warning_color)

    def append_error(self, message: str) -> None:
        """添加错误日志。"""
        self._append_colored(f"[错误] {message}", self._error_color)

    def append_success(self, message: str) -> None:
        """添加成功日志。"""
        self._append_colored(f"[成功] {message}", self._success_color)

    def _append_colored(self, message: str, color: QColor) -> None:
        """添加带颜色的文本。"""
        self.setTextColor(color)
        self.append(message)
        # 自动滚动到底部
        scrollbar = self.verticalScrollBar()
        if scrollbar:
            scrollbar.setValue(scrollbar.maximum())
