import os
import pandas as pd
from collections import Counter
from datetime import datetime

# Save fetched articles as CSV
def save_as_csv(keyword, records, output_dir):
    data = []
    for record in records:
        row = [
            keyword,
            record.get('PMID', ''),
            record.get('TI', ''),
            record.get('JT', ''),
            record.get('DP', ''),
            record.get('AB', ''),
        ]
        data.append(row)

    columns = ['Keyword', 'PMID', 'Title', 'Journal Title', 'Publication Date', 'Abstract']

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
