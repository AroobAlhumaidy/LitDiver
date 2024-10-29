
# 🤿  LitAutoDive (LITerature AUTOmation DIVEr) 🤿
This repository attempts to automate the process of literature retrieval and converting the results to CSV format. This setup is particularly useful for literature reviews, systematic reviews, and other academic research requiring access to a large volume of articles.

💻  Here, I started with Pubmed. utilizing its Edirect commandline tool for large-scale literature curation. Later-on I will add more databases 

The search strategy here is more inclusive, as it dose not specify time frame for the results (If you need one let me know 😬 ). Hence, it retrieves all available metadata. 


🌟 This script was tested with the manual method and both retrived the same papers .. you can find the example input and output in the directory "Example" 

All what you need to do is follow the steps: 

## Overview
The primary goal of these scripts is to automate the process of:
1. Retrieving PubMed articles based on keywords or phrases.
2. Converting the medline output to CSV for easier data handling and analysis.

## Pre-requisition
1. Ensure **NCBI Edirect** is installed. For installation instructions, see [NCBI Edirect Documentation](https://www.ncbi.nlm.nih.gov/books/NBK179288/).
2. Ensure Python 3.x is installed.

## QuickStart
Follow these steps to set up and run the scripts:

### Step 0: Get the codes 
download these two scripts (fetch_pubmed_lit.prod.v1.0.sh) and (convert.medline.csv.prod.v1.0.py) to your local machine, or just clone the repo: 

```bash
git clone https://github.com/AroobAlhumaidy/automate_litrature_search
cd automate_litrature_search
```

### Step 1: Prepare Keywords
Generate a list of keywords and save them in a text file named `mesh_keywords.txt` (or any name you prefer).

### Step2: Fetch the litrature
Run the fetching script using the following command 
    * replace mesh_keywords.txt with your keywords file name 

```bash
time bash fetch_pubmed_lit.Edirect.v1.4.sh mesh_keywords.txt
```
   *Note:* The `time` command is optional; it simply shows how long the script takes to complete.

#### Output
This script will generate a file for each keyword, these files are in medline format .. you can use the next step to convert these format to a regualr human readable csv files 📰

### Step 2: Convert to CSV

Use the Python script to convert MEDLINE files to a CSV format for easier reading and analysis. Run:
    * replace mesh_keywords.txt with your keywords file name 

```bash
time python3 convert.medline.csv.prod.v1.0.py mesh_keywords.txt
```

## Features

- **Automated Metadata Curation**: Retrieve all available PubMed records for each keyword without time restrictions.
- **CSV Conversion**: The final CSV file includes a keyword column and essential metadata columns. allowing you to easily merge or manipulate multiple files if needed.
- **Validated Accuracy**: Tested manually to confirm the accuracy of fetched papers.
- **Logging:** Errors encountered during the execution of PubMed searches are saved to `failed_keywords.txt`, allowing for easy retries.
- **Retry Mechanism:** Failed searches are reattempted by adjusting special characters.

## Notes
- This script does not filter results by date, pulling all available records.
- The keyword-based metadata retrieval provides comprehensive coverage across all specified terms.
