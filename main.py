import os
from utils.fetch import fetch_pubmed, save_as_xml
from utils.export import save_as_csv, generate_summary_report
from utils.pdf_downloader import download_pdfs
from utils.log import log_failures
from utils.config_loader import load_config

from Bio import Entrez

# Placeholder main script to orchestrate LitDiver modular functions
def main():
    config = load_config()
    Entrez.email = config['email']

    keyword_file = "keywords.txt"
    if not os.path.exists(keyword_file):
        print(f"Keyword file '{keyword_file}' not found.")
        return

    with open(keyword_file, "r") as f:
        keywords = [line.strip() for line in f if line.strip()]

    output_dir = config['output_dir']
    os.makedirs(output_dir, exist_ok=True)

    all_records = []
    for idx, keyword in enumerate(keywords, 1):
        print(f"[{idx}/{len(keywords)}] Processing keyword: {keyword}")
        try:
            search_term = keyword
            if config['date_range']:
                search_term += f" AND ({config['date_range']}[PDAT])"

            id_list, records = fetch_pubmed(search_term, max_results=config['max_results'])
            if not id_list:
                print(f"No records found for keyword: {keyword}")
                continue
            save_as_xml(keyword, id_list, output_dir)
            save_as_csv(keyword, records, output_dir)
            if config['download_pdfs']:
                download_pdfs(records, output_dir)
            all_records.extend(records)
        except Exception as e:
            print(f"Error processing keyword '{keyword}': {e}")

    generate_summary_report(all_records, output_dir)
    log_failures(output_dir)

if __name__ == "__main__":
    main()
