# 🤿  LitDiver (LITerature DIVER) 🤿  "beta version"
![License](https://img.shields.io/badge/license-MIT-blue)

---

### 📖 Overview:
LitDiver automates the process of retrieving literature from PubMed based on custom keywords, downloading available PDFs (PMC + Open Access), ris file for referincing tools, and generating a detailed summary report. Designed to streamlining papers search for small or large volume of scientific articles.

💻 LitDiver began with PubMed, In the future, support for additional databases will be added.

---
## 🚀 Features:
1. **📥 Download PDFs Automatically**  
   - Retrieves available PDFs via **PubMed Central (PMC)** or **Open Access (Unpaywall)**.
   - PDFs saved in an organized folder per run.

2. **📑 Reference File Export (RIS Format)**  
   - Automatically creates **.RIS files** for each keyword and a **combined deduplicated RIS**, compatible with **EndNote**, **Mendeley**, and other referencing tools.

3. **⚙️ Customizable Config File (YAML)**  
   - Set:
     - **Max results** per keyword (up to 10,000)
     - **Date range filter** (e.g., `2018:2024`)
     - **Output directory**
     - Enable/disable **PDF download**

4. **🔍 Flexible Search Scope Control**  
   - Customize **where PubMed searches**: Title, Abstract, MeSH terms, etc.
   - Use command-line option `--field` to set (`ti`, `ab`, `tiab`, `mh`, `tw`).  
     Example:
     ```bash
     python3 main.py --field tiab
     ```
   - Default: Broad search across all fields.

5. **📊 Summary Report Generation**  
   - Generates a Markdown report with:
     - Total articles fetched
     - PDF download attempts
     - Top journals
   - File: `summary_report.md` in the output folder.

6. **🪵 Download Failure Logging**  
   - Failed PDF downloads logged with reasons for easy troubleshooting.  
   - File: `download_failures.log`

---

### 📂 File structure:
| File/Folder         | Description                                |
|---------------------|--------------------------------------------|
| `main.py`           | Main script to run LitDiver                |
| `keywords.txt`      | List of search keywords (one per line)     |
| `config.yml`        | Config file (date range, results limit, etc.)   |
| `utils/`            | Contains modular scripts (fetch, export, pdfs, log, etc.) |
| `requirements.txt`  | Required Python packages                   |

---

### 🛠️ Requirements:
- Python 3.7+
- Install dependencies:
```bash
pip install -r requirements.txt
```

---

### 📑 Usage:
1. Edit `config.yml` with your preferences.
2. Add your keywords to `keywords.txt`.
3. Run:
```bash
python3 main.py
```
4. Results saved in `output_litdiver/` or the directory you specified in the `config.yml`

---

### 📊 Example Output:
- `results_<keyword>.csv` / `.xml` 
- Combined deduplicated XML
- `summary_report.md`
- `pdfs/` folder with downloaded PDFs
- `download_failures.log`