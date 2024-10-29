# This script aims to convert the medline format to csv format + adding the keyword in the first column 
# Log: this is prod version of the dev version: convert.medline.csv.v1.5.py
# Usage: time python3 convert.medline.csv.prod.v1.0.py mesh_keywords.txt
# By: Aroob ALhumaidy

import argparse
from Bio import Medline
import csv
import glob
import os
import re

def standardize(text):
    return re.sub(r'[^a-zA-Z0-9]', '', text).lower()

parser = argparse.ArgumentParser(description="Process Medline results files and add keywords from a specified file.")
parser.add_argument("keyword_file", type=str, help="Path to the keyword file (e.g., mesh_keywords.txt)")
args = parser.parse_args()

# Load keywords from input keyword sfile 
keywords_dict = {}
with open(args.keyword_file, "r") as keywords_file:
    for line in keywords_file:
        keyword = line.strip()
        # Standardize the keyword for matching
        standardized_keyword = standardize(keyword)
        keywords_dict[standardized_keyword] = keyword

# Fetch all files starting with "results_"
files = glob.glob("results_*.txt")
if not files:
    print("No files found starting with 'results_' in the current directory.")
else:
    for file_path in files:
        base_name = os.path.splitext(os.path.basename(file_path))[0].replace("results_", "")
        standardized_base_name = standardize(base_name)
        
        # Find the matching keyword in the keyword file 
        keyword_used = None
        for key in keywords_dict:
            if key in standardized_base_name:
                keyword_used = keywords_dict[key]
                break
        
        if not keyword_used:
            print(f"No matching keyword found for {file_path}. Skipping file.")
            continue

        print(f"Processing file: {file_path} with keyword: {keyword_used}")
        
        # Parse the fetched results from each Medline file
        with open(file_path, "r") as handle:
            records = Medline.parse(handle)
            records = list(records)

            if not records:
                print(f"No records found in {file_path}")
                continue
        
        # Output CSV file name
        csv_filename = os.path.splitext(file_path)[0] + ".csv"
        
        with open(csv_filename, "w", newline="", encoding="utf-8") as csvfile:
            writer = csv.writer(csvfile)
            writer.writerow(["Keyword", "PMID", "Title", "Journal", "Abstract", "Publication Date"])

            for record in records:
                pmid = record.get("PMID", "N/A")
                title = record.get("TI", "N/A")
                journal = record.get("JT", "N/A")
                abstract = record.get("AB", "N/A")
                pub_date = record.get("DP", "N/A")
                writer.writerow([keyword_used, pmid, title, journal, abstract, pub_date])

        print(f"Results saved to {csv_filename}")