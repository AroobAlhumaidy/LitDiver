import PySimpleGUI as sg
import yaml
import os

CONFIG_FILE = "config.yaml"

# Load current config
def load_config():
    if not os.path.exists(CONFIG_FILE):
        return {
            'email': '',
            'max_results': 100,
            'date_range': '',
            'download_pdfs': True,
            'output_dir': 'output_litdiver'
        }
    with open(CONFIG_FILE, 'r') as f:
        return yaml.safe_load(f)

# Save updated config
def save_config(updated_config):
    with open(CONFIG_FILE, 'w') as f:
        yaml.dump(updated_config, f)

def main():
    config = load_config()

    layout = [
        [sg.Text("LitDiver - Configuration", font=("Helvetica", 16), justification='center')],
        [sg.Text("Email:"), sg.InputText(config.get('email'), key='email')],
        [sg.Text("Max Results (up to 10000):"), sg.InputText(config.get('max_results'), key='max_results')],
        [sg.Text("Date Range (e.g., 2020:2024):"), sg.InputText(config.get('date_range'), key='date_range')],
        [sg.Text("Output Directory:"), sg.InputText(config.get('output_dir'), key='output_dir')],
        [sg.Checkbox("Download PDFs", default=config.get('download_pdfs', True), key='download_pdfs')],
        [sg.Button("Save and Run"), sg.Button("Cancel")]
    ]

    window = sg.Window("LitDiver GUI", layout)

    while True:
        event, values = window.read()
        if event in (sg.WIN_CLOSED, "Cancel"):
            break
        elif event == "Save and Run":
            updated_config = {
                'email': values['email'],
                'max_results': int(values['max_results']),
                'date_range': values['date_range'],
                'download_pdfs': values['download_pdfs'],
                'output_dir': values['output_dir']
            }
            save_config(updated_config)
            sg.popup("Config saved! You can now run LitDiver from the terminal.")
            break

    window.close()

if __name__ == "__main__":
    main()