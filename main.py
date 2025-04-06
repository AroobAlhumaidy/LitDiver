import os
import argparse
from utils.fetch import fetch_pubmed, save_as_xml
from utils.export import save_as_csv, generate_summary_report, save_combined_xml, save_combined_csv
from utils.pdf_downloader import download_pdfs
from utils.pdf_downloader import failure_log
from utils.log import log_failures
from utils.config_loader import load_config
from utils.export import save_combined_ris
from utils.export import save_as_ris


from Bio import Entrez

# Main script for LitDiver modular functions
def main():
    parser = argparse.ArgumentParser(description="LitDiver: Automated Literature Retrieval and PDF Download Tool")
    parser.add_argument("--field", type=str, default=None,
                        help="Optional field to search in: ti (Title), ab (Abstract), tiab (Title/Abstract), mh (MeSH), tw (Text Words). Default: all fields")
    parser.add_argument("--keywords", type=str, default=None,
                        help="Path to keyword file (e.g., MeSH.keywords.txt)")
    args = parser.parse_args()

    config = load_config()
    Entrez.email = config['email']

    keyword_file = args.keywords if args.keywords else "keywords.txt"
    if not os.path.exists(keyword_file):
        print(f"Keyword file '{keyword_file}' not found.")
        return

    with open(keyword_file, "r") as f:
        raw_keywords = [line.strip() for line in f if line.strip()]

    # Apply field tag if provided
    keywords = []
    if args.field and args.field.lower() != 'all':
        field_tag = args.field.lower()
        valid_fields = ['ti', 'ab', 'tiab', 'mh', 'tw']
        if field_tag not in valid_fields:
            print(f"Invalid field: {field_tag}. Valid options: {', '.join(valid_fields)}")
            return
        for kw in raw_keywords:
            parts = kw.split()
            tagged = [f"{p}[{field_tag}]" for p in parts]
            keywords.append(" ".join(tagged))
    else:
        keywords = raw_keywords

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
            for r in records:
                r['__keyword'] = keyword            
            all_records.extend(records)
        except Exception as e:
            print(f"Error processing keyword '{keyword}': {e}")
    if all_records:
        save_combined_xml(all_records, output_dir)
        save_combined_csv(all_records, output_dir)
        save_combined_ris(all_records, output_dir)

    generate_summary_report(all_records, output_dir)
    log_failures(output_dir, failure_log)
if __name__ == "__main__":
    main()