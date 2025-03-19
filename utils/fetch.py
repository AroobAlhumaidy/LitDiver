import time
import os
from tqdm import tqdm
from Bio import Entrez, Medline

# Fetch articles from PubMed
def fetch_pubmed(keyword, max_results=10000):
    handle = Entrez.esearch(db='pubmed', term=keyword, retmax=max_results)
    record = Entrez.read(handle)
    handle.close()

    total_hits = int(record["Count"])
    print(f"\n🔍 Total hits found: {total_hits} for keyword: {keyword}\n")

    id_list = record['IdList']
    print(f"📥 Fetching {len(id_list)} records in MEDLINE format...\n")

    all_records = []
    for start in tqdm(range(0, len(id_list), 1000),
                      desc="📖 Fetching MEDLINE",
                      unit="batch",
                      colour="green",
                      ncols=80,
                      bar_format='{l_bar}{bar}| {n_fmt}/{total_fmt} [{elapsed}<{remaining}]'):
        end = start + 1000
        batch_ids = id_list[start:end]
        handle = Entrez.efetch(db="pubmed", id=batch_ids, rettype="medline", retmode="text")
        records = list(Medline.parse(handle))
        handle.close()
        all_records.extend(records)
        time.sleep(0.5)

    return id_list, all_records

# Save fetched articles as XML
def save_as_xml(keyword, id_list, output_dir):
    filename = os.path.join(output_dir, f"results_{keyword.replace(' ', '_')}.xml")
    all_xml_data = ""
    for start in range(0, len(id_list), 1000):
        end = start + 1000
        batch_ids = id_list[start:end]
        handle = Entrez.efetch(db="pubmed", id=batch_ids, rettype="abstract", retmode="xml")
        xml_data = handle.read()
        handle.close()
        all_xml_data += xml_data
        time.sleep(0.5)

    with open(filename, "w", encoding="utf-8") as file:
        file.write(all_xml_data)

    print(f"XML results saved to {filename}")
