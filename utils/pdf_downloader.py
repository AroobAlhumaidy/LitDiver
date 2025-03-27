import os
import time
import requests

failure_log = []

# Download PDF from PMC if available
def download_pmc_pdf(pmc_id, output_dir, pmid, journal=''):
    pdf_url = f"https://www.ncbi.nlm.nih.gov/pmc/articles/PMC{pmc_id}/pdf/"
    headers = {
        'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/120.0 Safari/537.36'
    }
    try:
        response = requests.get(pdf_url, stream=True, headers=headers, timeout=10)
        if response.status_code == 200:
            pdf_path = os.path.join(output_dir, f"{pmid}_PMC{pmc_id}.pdf")
            with open(pdf_path, "wb") as f:
                for chunk in response.iter_content(chunk_size=8192):
                    f.write(chunk)
            print(f"Downloaded PMC PDF for PMID {pmid} to {pdf_path}")
        else:
            failure_log.append((pmid, f"PMC PDF download failed (status {response.status_code}) - Journal: {journal}"))
            print(f"Failed to download PMC PDF for PMID {pmid}, status code {response.status_code}")
    except Exception as e:
        failure_log.append((pmid, f"PMC PDF download error - {str(e)} - Journal: {journal}"))
        print(f"Error downloading PMC PDF for PMID {pmid}: {e}")

# Get Open Access PDF link via Unpaywall
def get_unpaywall_pdf(doi):
    email = "aroob.alhumaidy@gmail.com"
    url = f"https://api.unpaywall.org/v2/{doi}?email={email}"
    try:
        response = requests.get(url, timeout=10)
        if response.status_code == 200:
            data = response.json()
            best_oa = data.get('best_oa_location')
            oa_link = best_oa.get('url_for_pdf') if best_oa else None
            return oa_link
        else:
            print(f"Unpaywall API error for DOI {doi}, status: {response.status_code}")
            failure_log.append((doi, f"Unpaywall API status {response.status_code}"))
    except Exception as e:
        print(f"Unpaywall fetch failed for DOI {doi}: {e}")
        failure_log.append((doi, f"Unpaywall fetch error: {e}"))
    return None

# Download Open Access PDF via URL
def download_open_pdf(url, output_dir, pmid):
    try:
        response = requests.get(url, stream=True, timeout=15)
        if response.status_code == 200:
            pdf_path = os.path.join(output_dir, f"{pmid}_OA.pdf")
            with open(pdf_path, "wb") as f:
                for chunk in response.iter_content(chunk_size=8192):
                    f.write(chunk)
            print(f"Downloaded OA PDF for PMID {pmid} to {pdf_path}")
        else:
            print(f"Failed to download OA PDF for PMID {pmid}, status: {response.status_code}")
            failure_log.append((pmid, f"OA PDF download failed, status {response.status_code}"))
    except Exception as e:
        print(f"Error downloading OA PDF for PMID {pmid}: {e}")
        failure_log.append((pmid, f"OA PDF download error: {e}"))

# Process all records for PDF download
def download_pdfs(records, output_dir):
    pdf_dir = os.path.join(output_dir, "pdfs")
    os.makedirs(pdf_dir, exist_ok=True)

    for r in records:
        pmid = r.get('PMID', '')
        pmc_id = r.get('PMC', '').replace('PMC', '')
        journal = r.get('JT', '')

        if pmc_id:
            download_pmc_pdf(pmc_id, pdf_dir, pmid, journal)
        else:
            lid_value = r.get('LID', '').strip()
            if lid_value:
                doi_raw = lid_value.split()[0]
                doi = doi_raw.replace('[doi]', '').strip()
                if doi and '/' in doi:
                    pdf_url = get_unpaywall_pdf(doi)
                    if pdf_url:
                        download_open_pdf(pdf_url, pdf_dir, pmid)
                    else:
                        failure_log.append((pmid, f"No OA PDF found via Unpaywall - Journal: {journal}"))
            else:
                failure_log.append((pmid, f"No DOI found in LID field - Journal: {journal}"))