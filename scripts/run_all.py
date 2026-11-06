"""Run the full analysis from the command line (outside Jupyter or Colab).

Usage, from the repository folder:
    pip install -r requirements.txt
    python scripts/run_all.py
"""
import os
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")  # save figures without opening windows

SCRIPTS = Path(__file__).resolve().parent
os.chdir(SCRIPTS.parent)

namespace = {"__name__": "__main__"}
for part in ["01_setup_and_data.py", "02_analysis.py"]:
    print(f"\n========== Running {part} ==========\n")
    code = (SCRIPTS / part).read_text(encoding="utf-8")
    exec(compile(code, part, "exec"), namespace)
