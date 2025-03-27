import sys
import os
import yaml
from PySide6.QtWidgets import (
    QApplication, QWidget, QLabel, QLineEdit, QPushButton,
    QCheckBox, QFileDialog, QVBoxLayout, QHBoxLayout, QMessageBox
)
from PySide6.QtCore import Qt

CONFIG_PATH = "config.yaml"

def load_config():
    if os.path.exists(CONFIG_PATH):
        with open(CONFIG_PATH, 'r') as f:
            return yaml.safe_load(f)
    return {
        'email': '',
        'max_results': 100,
        'date_range': '',
        'output_dir': 'output_litdiver',
        'download_pdfs': True
    }

def save_config(config):
    with open(CONFIG_PATH, 'w') as f:
        yaml.dump(config, f)

class LitDiverGUI(QWidget):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("LitDiver Configuration")
        self.config = load_config()
        self.setup_ui()

    def setup_ui(self):
        layout = QVBoxLayout()

        # Email
        self.email_input = QLineEdit(self.config.get('email', ''))
        layout.addLayout(self._form_row("Email:", self.email_input))

        # Max Results
        self.max_results_input = QLineEdit(str(self.config.get('max_results', 100)))
        layout.addLayout(self._form_row("Max Results (up to 10000):", self.max_results_input))

        # Date Range
        self.date_input = QLineEdit(self.config.get('date_range', ''))
        layout.addLayout(self._form_row("Date Range (e.g., 2020:2024):", self.date_input))

        # Output Directory with file picker
        self.output_input = QLineEdit(self.config.get('output_dir', 'output_litdiver'))
        output_browse = QPushButton("Browse")
        output_browse.clicked.connect(self.pick_output_dir)
        output_layout = QHBoxLayout()
        output_layout.addWidget(QLabel("Output Directory:"))
        output_layout.addWidget(self.output_input)
        output_layout.addWidget(output_browse)
        layout.addLayout(output_layout)

        # PDF Download checkbox
        self.pdf_checkbox = QCheckBox("Download PDFs")
        self.pdf_checkbox.setChecked(self.config.get('download_pdfs', True))
        layout.addWidget(self.pdf_checkbox)

        # Buttons
        button_layout = QHBoxLayout()
        save_btn = QPushButton("Save and Run")
        save_btn.clicked.connect(self.save_and_run)
        cancel_btn = QPushButton("Cancel")
        cancel_btn.clicked.connect(self.close)
        button_layout.addWidget(save_btn)
        button_layout.addWidget(cancel_btn)
        layout.addLayout(button_layout)

        self.setLayout(layout)

    def _form_row(self, label_text, widget):
        layout = QHBoxLayout()
        label = QLabel(label_text)
        label.setFixedWidth(180)
        layout.addWidget(label)
        layout.addWidget(widget)
        return layout

    def pick_output_dir(self):
        dir_path = QFileDialog.getExistingDirectory(self, "Select Output Directory")
        if dir_path:
            self.output_input.setText(dir_path)

    def save_and_run(self):
        try:
            config = {
                'email': self.email_input.text().strip(),
                'max_results': int(self.max_results_input.text().strip()),
                'date_range': self.date_input.text().strip(),
                'output_dir': self.output_input.text().strip(),
                'download_pdfs': self.pdf_checkbox.isChecked()
            }
            save_config(config)
            QMessageBox.information(self, "Saved", "Configuration saved. You can now run LitDiver from the terminal.")
            self.close()
        except ValueError:
            QMessageBox.critical(self, "Error", "Max results must be a number.")

if __name__ == "__main__":
    app = QApplication(sys.argv)
    gui = LitDiverGUI()
    gui.show()
    sys.exit(app.exec())