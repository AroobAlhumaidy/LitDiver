import os
from utils.pdf_downloader import failure_log

# Save download failure log
def log_failures(output_dir, failures):
    if failures:
        log_path = os.path.join(output_dir, "download_failures.log")
        with open(log_path, "w") as log_file:
            log_file.write("ID - Error Message\n")
            for entry in failures:
                log_file.write(f"{entry[0]} - {entry[1]}\n")
        print(f"Download failure log saved to {log_path}")
