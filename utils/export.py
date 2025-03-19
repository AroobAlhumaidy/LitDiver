import os
import pandas as pd
from collections import Counter
from datetime import datetime
from Bio import Entrez
from tqdm import tqdm
import time

# Save fetched articles as CSV
def save_as_csv(keyword, records, output_dir):
    data = []
    for record in tqdm(records, desc=f"Processing records for: {keyword}"):
        row = [
            keyword,
            record.get('PMID', ''),
            record.get('TI', ''),
            record.get('JT', ''),
            record.get('DP', ''),
            record.get('AB', ''),
            "; ".join(record.get('AU', [])),
            record.get('AD', ''),
            "; ".join(record.get('PT', [])),
            "; ".join(record.get('MH', [])),
            record.get('PL', ''),
            record.get('VI', ''),
            record.get('IP', ''),
            record.get('PG', ''),
            record.get('LA', ''),
            record.get('PMC', ''),
            "; ".join(record.get('GR', [])),
            "; ".join(record.get('RN', [])),
            "; ".join(record.get('CI', [])),
            record.get('LID', ''),
            record.get('SO', ''),
            record.get('EDAT', ''),
            record.get('LR', '')
        ]
        data.append(row)

    columns = [
        'Keyword', 'PMID', 'Title', 'Journal Title', 'Publication Date', 'Abstract',
        'Authors', 'Author Address', 'Publication Types', 'MeSH Terms',
        'Country of Publication', 'Volume', 'Issue', 'Page Numbers', 'Language',
        'PubMed Central ID', 'Grant Info', 'CAS Registry Numbers',
        'Comments/Corrections', 'DOI', 'Source', 'Entry Date', 'Last Revision Date'
    ]

    df = pd.DataFrame(data, columns=columns)

    csv_filename = os.path.join(output_dir, f"results_{keyword.replace(' ', '_')}.csv")
    df.to_csv(csv_filename, index=False)

    print(f"Saved CSV to {csv_filename}")
    return df

# Generate Markdown dashboard summary report
def generate_summary_report(records, output_dir):
    total_records = len(records)
    pmc_downloads = sum(1 for r in records if r.get('PMC', ''))

    journals = [r.get('JT', '') for r in records if r.get('JT', '')]
    journal_counts = Counter(journals).most_common(5)

    report_lines = [
        f"# LitDiver Summary Report",
        f"**Date:** {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n",
        f"- **Total Articles Fetched:** {total_records}",
        f"- **PMC PDF Downloads Attempted:** {pmc_downloads}",
        f"\n## Top Journals:",
    ]
    for journal, count in journal_counts:
        report_lines.append(f"- {journal}: {count} articles")

    report_path = os.path.join(output_dir, "summary_report.md")
    with open(report_path, "w") as report_file:
        report_file.write("\n".join(report_lines))

    print(f"Summary report saved to {report_path}")

# Save combined deduplicated XML from all records
def save_combined_xml(all_records, output_dir):
    unique_pmids = list(set(record.get('PMID', '') for record in all_records if record.get('PMID', '')))

    print(f"\n📦 Fetching combined deduplicated XML for {len(unique_pmids)} unique PMIDs...")

    combined_xml = ""
    for start in tqdm(range(0, len(unique_pmids), 1000), desc="📄 Fetching Combined XML", unit="batch"):
        end = start + 1000
        batch_ids = unique_pmids[start:end]
        handle = Entrez.efetch(db="pubmed", id=batch_ids, rettype="abstract", retmode="xml")
        xml_data = handle.read()
        handle.close()
        combined_xml += xml_data
        time.sleep(0.5)

    xml_path = os.path.join(output_dir, "results_combined_deduplicated.xml")
    with open(xml_path, "w", encoding="utf-8") as xml_file:
        xml_file.write(combined_xml)

    print(f"Combined deduplicated XML saved to {xml_path}")