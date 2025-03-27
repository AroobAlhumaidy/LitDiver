import time
import os
from tqdm import tqdm
from Bio import Entrez, Medline
import shutil

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

    # Insert XSL stylesheet reference
    xml_with_style = '<?xml-stylesheet type="text/xsl" href="pubmed_style.xsl"?>\n' + all_xml_data

    with open(filename, "w", encoding="utf-8") as file:
        file.write(xml_with_style)

    # Copy the stylesheet into output dir
    shutil.copy("pubmed_style.xsl", output_dir)

    xml_dir = os.path.join(output_dir, "xml")
    os.makedirs(xml_dir, exist_ok=True)
    filename = os.path.join(xml_dir, f"results_{keyword.replace(' ', '_')}.xml")

    print(f"XML results saved to {filename}")