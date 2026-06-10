"""Ensure the repo root is importable during tests even without `pip install -e .`."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
