"""
SkinAge Atlas: Data Download and Ingestion Module
Downloads verified GEO datasets programmatically and extracts expression matrices and metadata.
"""

import os
import sys
import ssl
import gzip
import tarfile
import urllib.request
import pandas as pd
import numpy as np
from datetime import datetime

# Disable SSL verification for NCBI servers when local cert chain is incomplete
CTX = ssl.create_default_context()
CTX.check_hostname = False
CTX.verify_mode = ssl.CERT_NONE

HEADERS = {'User-Agent': 'Mozilla/5.0 (Bioinformatics SkinAge Atlas Pipeline)'}

DATA_RAW = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'data', 'raw'))
DATA_PROCESSED = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'data', 'processed'))
DATA_META = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'data', 'metadata'))

for d in [DATA_RAW, DATA_PROCESSED, DATA_META]:
    os.makedirs(d, exist_ok=True)


def download_file(url: str, dest_path: str, force: bool = False) -> str:
    """Download a remote file if it does not already exist."""
    if os.path.exists(dest_path) and not force and os.path.getsize(dest_path) > 0:
        print(f"[CACHE] {os.path.basename(dest_path)} already exists.")
        return dest_path
    print(f"[DOWNLOAD] Downloading {url} -> {dest_path}")
    req = urllib.request.Request(url, headers=HEADERS)
    with urllib.request.urlopen(req, context=CTX) as resp, open(dest_path, 'wb') as out_f:
        while True:
            chunk = resp.read(1024 * 1024)
            if not chunk:
                break
            out_f.write(chunk)
    print(f"[DONE] Saved {os.path.basename(dest_path)} ({os.path.getsize(dest_path) / 1024:.1f} KB)")
    return dest_path


def download_gene_info() -> str:
    """Download NCBI human gene_info table for Ensembl-to-Symbol mapping."""
    url = "https://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Mammalia/Homo_sapiens.gene_info.gz"
    dest = os.path.join(DATA_RAW, "Homo_sapiens.gene_info.gz")
    return download_file(url, dest)


def build_ensembl_symbol_map() -> dict:
    """Build Ensembl ID -> Gene Symbol dictionary from NCBI gene_info."""
    gene_info_path = download_gene_info()
    mapping_cache = os.path.join(DATA_META, "ensembl_to_symbol.csv")
    if os.path.exists(mapping_cache):
        print(f"[CACHE] Loading cached Ensembl-to-Symbol map from {mapping_cache}")
        df = pd.read_csv(mapping_cache)
        return dict(zip(df['ensembl_id'], df['symbol']))

    print("[PROCESSING] Parsing NCBI gene_info to extract Ensembl -> Symbol map...")
    ensembl_map = {}
    with gzip.open(gene_info_path, 'rt', encoding='utf-8', errors='replace') as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            if len(parts) < 6:
                continue
            symbol = parts[2]
            dbxrefs = parts[5]
            if 'Ensembl:' in dbxrefs:
                for item in dbxrefs.split('|'):
                    if item.startswith('Ensembl:'):
                        ens_id = item.split(':')[1].strip()
                        ensembl_map[ens_id] = symbol

    map_df = pd.DataFrame(list(ensembl_map.items()), columns=['ensembl_id', 'symbol'])
    map_df.to_csv(mapping_cache, index=False)
    print(f"[DONE] Indexed {len(ensembl_map):,} Ensembl-to-Symbol associations.")
    return ensembl_map


def download_gse226189():
    """
    Download GSE226189 (Discovery cohort: 82 RNA-seq samples of primary skin fibroblasts).
    Extracts raw counts for each sample and compiles metadata.
    """
    tar_url = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE226nnn/GSE226189/suppl/GSE226189_RAW.tar"
    matrix_url = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE226nnn/GSE226189/matrix/GSE226189_series_matrix.txt.gz"

    tar_path = os.path.join(DATA_RAW, "GSE226189_RAW.tar")
    matrix_path = os.path.join(DATA_RAW, "GSE226189_series_matrix.txt.gz")

    download_file(tar_url, tar_path)
    download_file(matrix_url, matrix_path)

    # Parse metadata from matrix header
    print("[METADATA] Parsing GSE226189 metadata...")
    with gzip.open(matrix_path, 'rt', encoding='utf-8', errors='replace') as f:
        matrix_lines = f.readlines()

    sample_ids, sample_titles, sexes, ages = [], [], [], []
    for line in matrix_lines:
        if line.startswith('!Sample_geo_accession'):
            sample_ids = [s.strip('"\n\r') for s in line.split('\t')[1:]]
        elif line.startswith('!Sample_title'):
            sample_titles = [s.strip('"\n\r') for s in line.split('\t')[1:]]
        elif line.startswith('!Sample_characteristics_ch1'):
            parts = [s.strip('"\n\r') for s in line.split('\t')[1:]]
            if len(parts) > 0 and parts[0].startswith('Sex:'):
                sexes = [p.replace('Sex:', '').strip() for p in parts]
            elif len(parts) > 0 and parts[0].startswith('age (years):'):
                ages = [int(p.replace('age (years):', '').strip()) for p in parts]

    meta_df = pd.DataFrame({
        'sample_id': sample_ids,
        'title': sample_titles,
        'sex': sexes,
        'age': ages,
        'cohort': 'GSE226189_Discovery',
        'tissue': 'Primary Skin Fibroblasts',
        'technology': 'RNA-seq (Illumina NovaSeq 6000)',
        'sun_exposure': 'Unspecified'
    })
    meta_df['age_group'] = pd.cut(
        meta_df['age'],
        bins=[0, 35, 60, 100],
        labels=['Young (<35)', 'Middle (35-60)', 'Old (>60)']
    )
    meta_csv = os.path.join(DATA_META, "gse226189_metadata.csv")
    meta_df.to_csv(meta_csv, index=False)
    print(f"[DONE] Saved GSE226189 metadata: {len(meta_df)} samples to {meta_csv}")

    # Extract raw counts from TAR
    counts_csv = os.path.join(DATA_PROCESSED, "gse226189_raw_counts.csv")
    if os.path.exists(counts_csv):
        print(f"[CACHE] Raw counts matrix already compiled: {counts_csv}")
        return meta_df, pd.read_csv(counts_csv, index_col=0)

    print("[PROCESSING] Extracting and merging per-sample geneCOUNT files from GSE226189_RAW.tar...")
    counts_dict = {}
    with tarfile.open(tar_path, 'r') as tar:
        for member in tar.getmembers():
            if 'geneCOUNT.txt.gz' in member.name:
                # e.g., GSM7067566_SKIN_AGE22_M_GR1_Gt073_geneCOUNT.txt.gz
                gsm = member.name.split('_')[0]
                f = tar.extractfile(member)
                decompressed = gzip.decompress(f.read()).decode('utf-8', errors='replace')
                sample_counts = {}
                for idx, line in enumerate(decompressed.splitlines()):
                    if idx == 0 or not line.strip():
                        continue
                    gene_id, count_val = line.strip().split('\t')
                    sample_counts[gene_id] = int(count_val)
                counts_dict[gsm] = sample_counts

    counts_df = pd.DataFrame(counts_dict)
    # Sort columns by metadata sample_id order
    counts_df = counts_df[[s for s in meta_df['sample_id'] if s in counts_df.columns]]
    counts_df.index.name = 'ensembl_id'
    counts_df.to_csv(counts_csv)
    print(f"[DONE] Compiled GSE226189 raw counts: {counts_df.shape[0]} genes x {counts_df.shape[1]} samples -> {counts_csv}")
    return meta_df, counts_df


def download_gse38308():
    """
    Download GSE38308 (Validation cohort: 42 samples, 21 matched pairs of
    pre-auricular sun-exposed vs post-auricular sun-protected human facial skin).
    """
    proc_url = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE38nnn/GSE38308/suppl/GSE38308%5Fprocessed%5Fdata%5Ffrom%5Fcommon%5Fprobe%5Fsequence%5Fbetween%5FGPL6106%5Fand%5FGPL6884%2Etxt%2Egz"
    probe_url = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE38nnn/GSE38308/suppl/GSE38308%5Fprobe%5Finfo%5FGPL6884%2Etxt%2Egz"
    mat_url = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE38nnn/GSE38308/matrix/GSE38308-GPL6884_series_matrix.txt.gz"
    mat6106_url = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE38nnn/GSE38308/matrix/GSE38308-GPL6106_series_matrix.txt.gz"

    proc_path = os.path.join(DATA_RAW, "GSE38308_processed_data.txt.gz")
    probe_path = os.path.join(DATA_RAW, "GSE38308_probe_info_GPL6884.txt.gz")
    mat_path = os.path.join(DATA_RAW, "GSE38308-GPL6884_series_matrix.txt.gz")
    mat6106_path = os.path.join(DATA_RAW, "GSE38308-GPL6106_series_matrix.txt.gz")

    download_file(proc_url, proc_path)
    download_file(probe_url, probe_path)
    download_file(mat_url, mat_path)
    download_file(mat6106_url, mat6106_path)

    # Parse probe info to get probe -> symbol map
    print("[PROCESSING] Parsing GSE38308 probe mapping...")
    probe_map = {}
    with gzip.open(probe_path, 'rt', encoding='utf-8', errors='replace') as f:
        for idx, line in enumerate(f):
            if idx == 0:
                continue
            parts = line.strip().split('\t')
            if len(parts) >= 2:
                probe_id, symbol = parts[0], parts[1]
                if symbol and symbol != '':
                    probe_map[probe_id] = symbol

    # Parse metadata from matrix headers
    print("[METADATA] Parsing GSE38308 metadata...")
    ages_dict = {}
    for mp in [mat_path, mat6106_path]:
        with gzip.open(mp, 'rt', encoding='utf-8', errors='replace') as f:
            cur_titles = []
            for line in f:
                if line.startswith('!Sample_title'):
                    cur_titles = [s.strip('"\n\r') for s in line.split('\t')[1:]]
                elif line.startswith('!Sample_characteristics_ch1') and 'age:' in line:
                    cur_ages = [int(s.strip('"\n\r').replace('age:', '').strip()) for s in line.split('\t')[1:]]
                    for title, age in zip(cur_titles, cur_ages):
                        # title example: 'sun-exposed pre-auricular skin 10A'
                        donor_code = title.split()[-1]  # '10A'
                        donor_id = donor_code[:-1]      # '10'
                        ages_dict[donor_code] = (donor_id, age)

    # Read processed expression table
    print("[PROCESSING] Loading GSE38308 processed expression data...")
    expr_df = pd.read_csv(proc_path, sep='\t', compression='gzip', index_col=0)

    # Build metadata dataframe matching expr_df columns
    meta_records = []
    for col in expr_df.columns:
        donor_id, age = ages_dict.get(col, (col[:-1], np.nan))
        exposure = "Sun-Exposed" if col.endswith('A') else "Sun-Protected"
        site = "pre-auricular" if col.endswith('A') else "post-auricular"
        meta_records.append({
            'sample_id': col,
            'donor_id': f"Donor_{donor_id}",
            'age': age,
            'sex': 'Female',
            'sun_exposure': exposure,
            'anatomical_site': site,
            'tissue': 'Human Facial Skin',
            'cohort': 'GSE38308_Validation',
            'technology': 'Microarray (Illumina HumanWG-6 v3.0)'
        })
    meta_df = pd.DataFrame(meta_records)
    meta_csv = os.path.join(DATA_META, "gse38308_metadata.csv")
    meta_df.to_csv(meta_csv, index=False)
    print(f"[DONE] Saved GSE38308 metadata: {len(meta_df)} samples to {meta_csv}")

    # Map probe IDs to gene symbols and collapse duplicates by mean
    expr_df['symbol'] = expr_df.index.map(probe_map)
    expr_clean = expr_df.dropna(subset=['symbol']).groupby('symbol').mean()
    clean_csv = os.path.join(DATA_PROCESSED, "gse38308_expression_symbols.csv")
    expr_clean.to_csv(clean_csv)
    print(f"[DONE] Saved GSE38308 gene expression matrix: {expr_clean.shape[0]} genes x {expr_clean.shape[1]} samples -> {clean_csv}")
    return meta_df, expr_clean


if __name__ == '__main__':
    print(f"=== Starting SkinAge Atlas Data Download [{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}] ===")
    ensembl_map = build_ensembl_symbol_map()
    meta_disc, counts_disc = download_gse226189()
    meta_val, expr_val = download_gse38308()
    print("=== All datasets downloaded and parsed successfully! ===")
