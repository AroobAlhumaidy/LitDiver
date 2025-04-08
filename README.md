# 🤿  LitDiver (LITerature DIVER) 🤿  "beta version 2.2"
![License](https://img.shields.io/badge/license-MIT-blue)

---

### 📖 Overview:
LitDiver automates the process of retrieving scientific literature from PubMed using custom keywords (search queries), downloading available PDFs (from PMC and open-access sources), exporting citation data as a RIS file for reference managers, and generating a detailed summary report. It is designed to streamline literature searches for both small- and large-scale research projects

💻 LitDiver began with PubMed, In the future, support for additional databases will be added.

---
## 🚀 Features:
1. **✅ Interactive GUI** 
   - User-friendly interface to configure and run searches  
   - Real-time progress tracking and log updates  
   - Output folder opens automatically at the end

2. **🧠 Rich Metadata Retrieval**. 

Fetches detailed article metadata including:  

| Keyword | PMID | Title | Journal Title | Publication Date | Abstract | Authors | Author Address | Publication Types | MeSH Terms | Country of Publication | Volume | Issue | Page Numbers | Language | PubMed Central ID | Grant Info | CAS Registry Numbers | Comments/Corrections | DOI | Source | Entry Date | Last Revision Date |
|---------|------|-------|----------------|------------------|----------|---------|----------------|--------------------|-------------|-------------------------|--------|--------|----------------|----------|--------------------|------------|------------------------|------------------------|------|--------|-------------|----------------------|

All results are saved in structured .csv formats, ready for analysis and filtering.

3. **📥 Download PDFs Automatically**  
   - Retrieves available PDFs via PubMed Central (PMC) or Open Access (Unpaywall).
   - PDFs saved in an organized folder per run.

4. **📑 Reference File Export (RIS Format)**  
   - Automatically creates .RIS files for each keyword and a combined deduplicated RIS, compatible with EndNote, Mendeley, and other referencing tools.

5. **⚙️ Customizable Config File (YAML)**  
   - Set:
     - Max results: per keyword (up to 10,000)
     - Date range filter: (e.g., `2018:2024`)
     - Output directoryy
     - Enable/disable PDF download

6. **🔍 Flexible Search Scope Control**  
   - Customize where PubMed searches: Title, Abstract, MeSH terms, etc.
   - Check out options with `--help`
     Example:
     ```bash
     python3 main.py --help
     ```
   - Default: Broad search across all fields.

7. **📊 Summary Report Generation**  
   - Generates a Markdown report with:
     - Total articles fetched
     - PDF download attempts
     - Top journals
   - File: `summary_report.md` in the output folder.

8. **🪵 Download Failure Logging**  
   - Failed PDF downloads logged with reasons for easy troubleshooting.  
   - File: `download_failures.log`

---

## 🖥️ Cross-Platform Support

- **Linux (Ubuntu, Debian):** `gui.py`
- **macOS:** `gui_MacOS.py`

---


### 📂 File structure:
| File/Folder         | Description                                |
|---------------------|--------------------------------------------|
| `gui.py`            | # GUI for Linux                |
| `gui_MacOS.py`            | # GUI for macOS                |
| `keywords.txt`      | List of search keywords (one per line)     |
| `config.yml`        | Config file (date range, results limit, etc.)   |
| `utils/`            | Contains modular scripts (fetch, export, pdfs, log, etc.) |
| `requirements.txt`  | Required Python packages                   |

---

### 💻 Installation 
```
git clone https://github.com/AroobAlhumaidy/LitDiver
```
---
### 🛠️ Requirements:
- Python 3.7+
- Install dependencies:
```bash
cd LitDiver
pip install -r requirements.txt
```
---
### 🚀 How It Works 
#### 1. Open LitDiver from the Terminal
##### • Linux (Ubuntu/Debian)
   1. Open the Terminal app
   2. Navigate to the folder where you downloaded LitDiver:
   ```bash
   cd /path/to/LitDiver
   ```
   3. Start the app:
   ```bash
   python3 gui.py
   ```
##### • macOS
   1. Open the Terminal app
   2. Navigate to the LitDiver folder
   ```sh
   cd /path/to/LitDiver
   ```
   3. Start the app:
   ```sh
   python3 gui_MacOS.py
   ```
#### 2. Fill in Your Search Settings
   - Select your keywords file (.txt, one keyword or phrase per line)
   - Choose an output folder where results will be saved
   - Set how many articles you want (example: 100, maximum is 10,000 )
   - (Optional) Add a date range (example: 2020:2024)
   - *(Optional)* Choose a **field to search in** (example: title, abstract, MeSH)  
     - **Default:** All fields LitDiver will search across the full PubMed record, including title, abstract, MeSH terms, and more
   - Keep "Download PDFs" checked to fetch open-access full-texts

#### 3. Click ‘Save Config and Run’
   - A live progress bar and log window will keep you updated
   - You can safely stop the run at any time

#### 4. Get Your Results
When the search is complete, LitDiver will:

- Open your **output folder** automatically  
- Save all available **PDFs**  
- Export a **metadata spreadsheet** (`.csv`)
- Export **citation files** in `.ris` format  
- Generate a **summary report** in `summary_report.md`  
- PDFs that could not be downloaded are logged with reasons in `download_failures.log`

---
### 📝 Notes
- Ensure your keywords.txt file is formatted with one keyword/search phrase per line.
- An internet connection is required for downloading articles and PDFs.
- No NCBI API key is required, but rate limits apply.

---
### 🧑‍💻 Contributing
Pull requests are welcome, Please open an issue first to discuss what you’d like to change.

