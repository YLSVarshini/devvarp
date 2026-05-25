import os, requests, gzip, shutil
from pathlib import Path

BASE_DIR = Path.home() / "devvarp_project"
RAW_DIR  = BASE_DIR / "data" / "raw"
RAW_DIR.mkdir(parents=True, exist_ok=True)

def download_file(url, dest, name):
    if Path(dest).exists():
        print(f"[SKIP] {name} already exists")
        return
    print(f"[DOWN] {name}...")
    r = requests.get(url, stream=True, timeout=120)
    with open(dest, 'wb') as f:
        for chunk in r.iter_content(1024*1024):
            f.write(chunk)
    print(f"[DONE] {name}")

def decomp(src, dst, name):
    if Path(dst).exists():
        print(f"[SKIP] {name} already decompressed")
        return
    print(f"[DECOMP] {name}...")
    with gzip.open(src,'rb') as fi, open(dst,'wb') as fo:
        shutil.copyfileobj(fi, fo)
    print(f"[DONE] {name}")

download_file(
    "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz",
    RAW_DIR/"variant_summary.txt.gz", "ClinVar")
decomp(RAW_DIR/"variant_summary.txt.gz", RAW_DIR/"variant_summary.txt", "ClinVar")

download_file(
    "https://downloads.wenglab.org/Registry-V3/GRCh38-cCREs.bed",
    RAW_DIR/"GRCh38-cCREs.bed", "ENCODE cCREs")

download_file(
    "https://github.com/obophenotype/human-phenotype-ontology/releases/latest/download/genes_to_phenotype.txt",
    RAW_DIR/"genes_to_phenotype.txt", "HPO annotations")

download_file(
    "https://github.com/obophenotype/human-phenotype-ontology/releases/latest/download/hp.obo",
    RAW_DIR/"hp.obo", "HPO OBO")

print("\nDone. Check ~/devvarp_project/data/raw/ for all files.")
