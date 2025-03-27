import sys
import os
import platform
import yaml
import subprocess
import re
from PySide6.QtWidgets import (
    QApplication, QWidget, QLabel, QLineEdit, QVBoxLayout,
    QPushButton, QCheckBox, QFileDialog, QMessageBox, QHBoxLayout,
    QPlainTextEdit, QComboBox, QProgressBar
)
from PySide6.QtCore import Qt, QProcess

class LitDiverGUI(QWidget):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("LitDiver - GUI Configuration")
        self.setFixedWidth(500)
        self.process = QProcess(self)
        self.process.readyReadStandardOutput.connect(self.read_stdout)
        self.process.readyReadStandardError.connect(self.read_stderr)
        self.process.finished.connect(self.process_finished)
        self.init_ui()
        self.total_keywords = 0
        self.keywords_done = 0

    def init_ui(self):
        layout = QVBoxLayout()

        # Email input
        self.email_label = QLabel("Email Address:")
        self.email_input = QLineEdit("")
        layout.addWidget(self.email_label)
        layout.addWidget(self.email_input)

        # Max Results
        self.max_results_label = QLabel("Max Results (up to 10000):")
        self.max_results_input = QLineEdit("100")
        layout.addWidget(self.max_results_label)
        layout.addWidget(self.max_results_input)

        # Date Range
        self.date_range_label = QLabel("Date Range (e.g., 2020:2024):")
        self.date_range_input = QLineEdit("")
        layout.addWidget(self.date_range_label)
        layout.addWidget(self.date_range_input)

        # Search Field Selector
        self.field_label = QLabel("Search Field (optional):")
        self.field_dropdown = QComboBox()
        self.field_dropdown.addItems(["All Fields", "ti", "ab", "tiab", "mh", "tw"])
        layout.addWidget(self.field_label)
        layout.addWidget(self.field_dropdown)

        # Keywords File
        self.keyword_file_label = QLabel("Select Keywords File:")
        self.keyword_file_input = QLineEdit()
        self.keyword_file_button = QPushButton("Browse...")
        self.keyword_file_button.clicked.connect(self.select_keyword_file)
        layout.addWidget(self.keyword_file_label)
        layout.addWidget(self.keyword_file_input)
        layout.addWidget(self.keyword_file_button)

        # Output Directory
        self.output_dir_label = QLabel("Select Output Directory:")
        self.output_dir_input = QLineEdit("output_litdiver")
        self.output_dir_button = QPushButton("Browse...")
        self.output_dir_button.clicked.connect(self.select_output_dir)
        layout.addWidget(self.output_dir_label)
        layout.addWidget(self.output_dir_input)
        layout.addWidget(self.output_dir_button)

        # PDF Checkbox
        self.pdf_checkbox = QCheckBox("Download PDFs")
        self.pdf_checkbox.setChecked(True)
        layout.addWidget(self.pdf_checkbox)

        # Run and Stop Buttons
        button_layout = QHBoxLayout()
        self.run_button = QPushButton("Save Config and Run")
        self.run_button.clicked.connect(self.save_and_run)
        self.stop_button = QPushButton("Stop")
        self.stop_button.clicked.connect(self.stop_process)
        self.stop_button.setEnabled(False)
        button_layout.addWidget(self.run_button)
        button_layout.addWidget(self.stop_button)
        layout.addLayout(button_layout)

        # Progress Bar
        self.progress_bar = QProgressBar()
        self.progress_bar.setValue(0)
        self.progress_bar.setTextVisible(True)
        layout.addWidget(self.progress_bar)

        # Log output area
        self.log_output = QPlainTextEdit()
        self.log_output.setReadOnly(True)
        self.log_output.setPlaceholderText("Logs will appear here...")
        layout.addWidget(self.log_output)

        self.setLayout(layout)

    def select_keyword_file(self):
        file_path, _ = QFileDialog.getOpenFileName(self, "Select Keywords File", "", "Text Files (*.txt)")
        if file_path:
            self.keyword_file_input.setText(file_path)

    def select_output_dir(self):
        dir_path = QFileDialog.getExistingDirectory(self, "Select Output Directory")
        if dir_path:
            self.output_dir_input.setText(dir_path)

    def save_and_run(self):
        try:
            max_results = int(self.max_results_input.text())
            if max_results <= 0 or max_results > 10000:
                raise ValueError
        except ValueError:
            QMessageBox.critical(self, "Error", "Max Results must be a number between 1 and 10000.")
            return

        keyword_file = self.keyword_file_input.text().strip()
        if not keyword_file or not os.path.exists(keyword_file):
            QMessageBox.warning(self, "Missing Input", "Please select a valid keyword file.")
            return

        email = self.email_input.text().strip()
        if not email:
            QMessageBox.warning(self, "Missing Input", "Please enter a valid email address.")
            return

        field_selection = self.field_dropdown.currentText()
        field = None if field_selection == "All Fields" else field_selection

        with open(keyword_file, "r") as f:
            self.total_keywords = sum(1 for line in f if line.strip())
            self.progress_bar.setMaximum(self.total_keywords)
            self.keywords_done = 0

        config = {
            "email": email,
            "max_results": max_results,
            "date_range": self.date_range_input.text(),
            "output_dir": self.output_dir_input.text(),
            "download_pdfs": self.pdf_checkbox.isChecked(),
            "field": field_selection if field else "all"
        }
        with open("config.yml", "w") as f:
            yaml.dump(config, f)

        self.log_output.clear()
        self.stop_button.setEnabled(True)
        self.run_button.setEnabled(False)
        self.progress_bar.setValue(0)

        cmd = ["main.py", "--keywords", keyword_file]
        if field:
            cmd.extend(["--field", field])

        self.process.start("python3", cmd)

    def read_stdout(self):
        output = self.process.readAllStandardOutput().data().decode()
        self.log_output.appendPlainText(output)

        match = re.search(r"\[(\d+)/(\d+)] Processing keyword", output)
        if match:
            self.keywords_done = int(match.group(1))
            self.progress_bar.setValue(self.keywords_done)

    def read_stderr(self):
        output = self.process.readAllStandardError().data().decode()
        self.log_output.appendPlainText(output)

    def process_finished(self):
        self.stop_button.setEnabled(False)
        self.run_button.setEnabled(True)
        self.log_output.appendPlainText("✅ LitDiver has finished.")
        self.open_output_dir(self.output_dir_input.text())

        summary = f"✅ Run Complete\nProcessed {self.keywords_done} of {self.total_keywords} keywords."
        QMessageBox.information(self, "Summary", summary)

    def stop_process(self):
        if self.process and self.process.state() == QProcess.Running:
            self.process.terminate()
            self.log_output.appendPlainText("❌ LitDiver was stopped by the user.")
            self.stop_button.setEnabled(False)
            self.run_button.setEnabled(True)

    def open_output_dir(self, path):
        if platform.system() == "Windows":
            os.startfile(path)
        elif platform.system() == "Darwin":
            subprocess.run(["open", path])
        else:
            subprocess.run(["xdg-open", path])

if __name__ == "__main__":
    app = QApplication(sys.argv)
    gui = LitDiverGUI()
    gui.show()
    sys.exit(app.exec())
