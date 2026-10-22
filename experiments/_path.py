import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RESULTS = os.path.join(ROOT, "results")
os.makedirs(RESULTS, exist_ok=True)
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
