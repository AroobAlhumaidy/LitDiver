import sys
import os
import platform
import yaml
import subprocess
from PySide6.QtWidgets import (
    QApplication, QWidget, QLabel, QLineEdit, QVBoxLayout,
    QPushButton, QCheckBox, QFileDialog, QMessageBox, QHBoxLayout
)
from PySide6.QtCore import Qt, QTimer

class LitDiverGUI(QWidget):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("LitDiver - GUI Configuration")
        self.setFixedWidth(400)
        self.process = None
        self.init_ui()

    def init_ui(self):
        layout = QVBoxLayout()

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
        self.save_button = QPushButton("Save Config and Run")
        self.save_button.clicked.connect(self.save_and_run)
        self.stop_button = QPushButton("Stop")
        self.stop_button.clicked.connect(self.stop_process)
        self.stop_button.setEnabled(False)
        button_layout.addWidget(self.save_button)
        button_layout.addWidget(self.stop_button)
        layout.addLayout(button_layout)

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

        config = {
            "email": "your-email@example.com",
            "max_results": max_results,
            "date_range": self.date_range_input.text(),
            "output_dir": self.output_dir_input.text(),
            "download_pdfs": self.pdf_checkbox.isChecked(),
        }
        with open("config.yaml", "w") as f:
            yaml.dump(config, f)

        try:
            self.process = subprocess.Popen(["python3", "main.py", "--keywords", keyword_file])
            self.stop_button.setEnabled(True)
            QTimer.singleShot(2000, lambda: self.check_process(config["output_dir"]))
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Failed to run main.py: {e}")

    def check_process(self, output_dir):
        if self.process and self.process.poll() is None:
            QTimer.singleShot(1000, lambda: self.check_process(output_dir))
        else:
            self.stop_button.setEnabled(False)
            QMessageBox.information(self, "Success", "LitDiver has finished running.")
            self.open_output_dir(output_dir)

    def stop_process(self):
        if self.process and self.process.poll() is None:
            self.process.terminate()
            self.stop_button.setEnabled(False)
            QMessageBox.information(self, "Stopped", "LitDiver was stopped by the user.")

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