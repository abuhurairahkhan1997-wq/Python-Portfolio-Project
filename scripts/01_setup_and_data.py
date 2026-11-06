# ==================================================================
# PART 1 — Environment, data download and provenance check
# Scanpy re-analysis of one healthy human PBMC sample
# Data:   GSM5227130 (PBMC_healthy [52BA23-00]) from GEO series GSE171555
# Source: Yu C, et al. Med 2021;2(6):755-772.e5. doi:10.1016/j.medj.2021.04.008
# Author: Abuhurairah Rahman Pathan
# ==================================================================
# Google Colab: run this part once, restart the session when asked
# (Runtime > Restart session), then run everything from the top again.
# Local: create a Python 3.12 environment, `pip install -r requirements.txt`,
# then run `python scripts/run_all.py` from the repository folder.
# ==================================================================

import sys, os, subprocess, shutil, tarfile, hashlib, urllib.request
from pathlib import Path
import importlib.metadata as md_version

# ---------- 0. Pinned environment ----------
# These are the exact versions that produced the results in this repository.
PINNED = {
    "numpy": "2.0.2",
    "pandas": "2.2.3",
    "scipy": "1.16.3",
    "scikit-learn": "1.6.1",
    "matplotlib": "3.10.9",
    "h5py": "3.16.0",
    "anndata": "0.12.19",
    "scanpy": "1.11.5",
    "numba": "0.68.1",
    "llvmlite": "0.50.0",
    "pynndescent": "0.5.13",
    "umap-learn": "0.5.7",
    "igraph": "1.0.0",
    "scikit-image": "0.25.2",
}

def installed(pkg):
    try:
        return md_version.version(pkg)
    except md_version.PackageNotFoundError:
        return None

mismatched = {p: (installed(p), v) for p, v in PINNED.items() if installed(p) != v}

if mismatched:
    print("Installing pinned versions for:")
    for p, (have, want) in mismatched.items():
        print(f"  {p}: {have} -> {want}")
    subprocess.run(
        [sys.executable, "-m", "pip", "install", "-q",
         *[f"{p}=={v}" for p, v in PINNED.items()]],
        check=True,
    )
    print("\nDone. Restart the session (Colab: Runtime > Restart session), "
          "then run everything from the top again.")
    raise SystemExit("Restart required after installing pinned packages.")

print("✓ All pinned package versions are installed.")
if sys.version_info[:2] != (3, 12):
    print(f"⚠ Python {sys.version.split()[0]} detected; the reference run used Python 3.12. "
          "Part 2's cluster check will tell you if results differ.")

# ---------- 1. Imports ----------
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scanpy as sc

try:
    from IPython.display import display
except ImportError:
    display = print

try:
    from google.colab import files
    IN_COLAB = True
except ImportError:
    files = None
    IN_COLAB = False

np.random.seed(42)
sc.settings.verbosity = 1
sc.set_figure_params(dpi=100, facecolor="white")

print("Python:", sys.version.split()[0])
for p in ["scanpy", "anndata", "numpy", "pandas", "scipy", "igraph"]:
    print(f"{p}: {md_version.version(p)}")

# ---------- 2. Project folders ----------
if IN_COLAB:
    ROOT = Path("/content")
else:
    ROOT = Path.cwd()
    if ROOT.name in {"notebooks", "scripts"}:
        ROOT = ROOT.parent
os.chdir(ROOT)

# ---------- 3. Dataset provenance ----------
SAMPLE_ID   = "GSM5227130"
SAMPLE_NAME = "PBMC_healthy [52BA23-00]"
GEO_SERIES  = "GSE171555"
FILE_PREFIX = "GSM5227130_52BA23-00_"
GEO_URL     = "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM5227nnn/GSM5227130/suppl/"
CITATION    = ("Yu C, Littleton S, Giroux NS, et al. Mucosal-associated invariant T cell "
               "responses differ by sex in COVID-19. Med. 2021;2(6):755-772.e5. "
               "doi:10.1016/j.medj.2021.04.008")

# MD5 checksums of the three processed GEO files for this sample
EXPECTED_MD5 = {
    "barcodes.tsv.gz": "542a03cbcc30842a55b7bc462cbb3f7f",
    "features.tsv.gz": "8361c1748151e98503e656c4e2d89a48",
    "matrix.mtx.gz":   "9b4cadcc807f8b3c68d7cf31a9f630ca",
}

DATA_DIR = ROOT / "data" / SAMPLE_ID
DATA_DIR.mkdir(parents=True, exist_ok=True)

def md5(path, chunk=1 << 20):
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(chunk), b""):
            h.update(block)
    return h.hexdigest()

def sample_files_ok():
    return all((DATA_DIR / f"{FILE_PREFIX}{e}").exists() for e in EXPECTED_MD5)

def download_from_geo():
    for ending in EXPECTED_MD5:
        name = f"{FILE_PREFIX}{ending}"
        target = DATA_DIR / name
        if target.exists():
            continue
        print(f"  downloading {name} ...")
        tmp = target.with_suffix(target.suffix + ".part")
        urllib.request.urlretrieve(GEO_URL + name, tmp)
        tmp.rename(target)

def extract_from_tar(archive_path):
    with tarfile.open(archive_path) as tar:
        for member in tar.getmembers():
            name = Path(member.name).name
            if member.isfile() and name.startswith(FILE_PREFIX):
                with tar.extractfile(member) as src, open(DATA_DIR / name, "wb") as dst:
                    shutil.copyfileobj(src, dst)

# ---------- 4. Get the data ----------
# Order: files already present -> download from GEO -> local .tar -> Colab upload
if sample_files_ok():
    print(f"\nSample files already present in {DATA_DIR.relative_to(ROOT)}/")
else:
    try:
        print("\nDownloading the three GSM5227130 files from GEO ...")
        download_from_geo()
    except Exception as err:
        print(f"GEO download failed ({err}). Falling back to a .tar archive.")
        tars = sorted(ROOT.glob("*.tar")) + sorted((ROOT / "data").glob("*.tar"))
        if tars:
            archive_path = tars[0]
        elif IN_COLAB:
            print("Click 'Choose Files' and select a .tar containing the GSM5227130 files.")
            uploaded = files.upload()
            if not uploaded:
                raise RuntimeError("No file was uploaded.")
            archive_path = ROOT / list(uploaded.keys())[0]
        else:
            raise FileNotFoundError(
                f"Download the three {FILE_PREFIX}* files from {GEO_URL} into {DATA_DIR} and re-run."
            )
        print(f"Extracting {archive_path.name}")
        extract_from_tar(archive_path)
    if not sample_files_ok():
        raise FileNotFoundError(f"Could not find all three {FILE_PREFIX}* files.")

# ---------- 5. Verify provenance by MD5 ----------
print("\nChecking file integrity:")
for ending, expected in EXPECTED_MD5.items():
    f = DATA_DIR / f"{FILE_PREFIX}{ending}"
    got = md5(f)
    if got != expected:
        raise ValueError(f"MD5 mismatch for {f.name}: {got} (expected {expected}). "
                         "Delete the file and re-run to download it again.")
    print(f"  ✓ {f.name}  {got}")

print(f"\n✓ Part 1 complete. Data verified in {DATA_DIR.relative_to(ROOT)}/")
