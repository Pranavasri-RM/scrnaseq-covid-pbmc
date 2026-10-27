import sys
from pathlib import Path

# Make the scripts in bin/ importable as modules (e.g. `import annotate`)
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "bin"))