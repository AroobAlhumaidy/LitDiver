import PySimpleGUI as sg
import os
from subprocess import run

# --- GUI Layout ---
layout = [
    [sg.Text("LitDiver - Literature Fetching Tool", font=("Helvetica", 16), justification='center')],
    [sg.Text("Select Keywords File (.txt):"), sg.InputText(key="-KEYWORDS-"), sg.FileBrowse(file_types=(("Text Files", "*.txt"),))],
    [sg.Text("Search Field:"), sg.Combo(["All Fields", "Title (ti)", "Abstract (ab)", "Title/Abstract (tiab)", "MeSH Terms (mh)", "Text Words (tw)"], default_value="All Fields", key="-FIELD-")],
    [sg.Text("Max Results per Query:"), sg.InputText("1000", key="-MAX_RESULTS-", size=(10,1))],
    [sg.Text("Date Range (e.g., 2020:2023):"), sg.InputText("", key="-DATE_RANGE-", size=(20,1))],
    [sg.Checkbox("Download PDFs", default=True, key="-PDFS-")],
    [sg.Text("Output Directory:"), sg.InputText(key="-OUTPUT_DIR-"), sg.FolderBrowse()],
    [sg.Button("Start Fetching", size=(15,1), button_color=("white", "green")), sg.Button("Exit", size=(10,1))],
    [sg.Output(size=(100, 20), key="-OUTPUT-")]
]

window = sg.Window("LitDiver GUI", layout, resizable=True)

# --- Event Loop ---
while True:
    event, values = window.read()
    if event in (sg.WIN_CLOSED, "Exit"):
        break
    elif event == "Start Fetching":
        keyword_file = values["-KEYWORDS-"]
        field = values["-FIELD-"]
        max_results = values["-MAX_RESULTS-"]
        date_range = values["-DATE_RANGE-"]
        pdfs = values["-PDFS-"]
        output_dir = values["-OUTPUT_DIR-"] or "output_litdiver"

        # Normalize field code
        field_map = {
            "All Fields": None, "Title (ti)": "ti", "Abstract (ab)": "ab",
            "Title/Abstract (tiab)": "tiab", "MeSH Terms (mh)": "mh", "Text Words (tw)": "tw"
        }
        field_code = field_map.get(field)

        # Build command to run main.py with arguments
        command = ["python3", "main.py"]
        if field_code:
            command += ["--field", field_code]

        # Prepare config.yaml dynamically
        with open("config.yaml", "w") as cfg:
            cfg.write(f"""
email: your-email@example.com
max_results: {max_results}
date_range: '{date_range}'
download_pdfs: {str(pdfs).lower()}
output_dir: {output_dir}
""")

        print("\n--- Fetching Literature ---\n")
        run(command)
        print("\n✅ Done. Check the output directory.")

window.close()
