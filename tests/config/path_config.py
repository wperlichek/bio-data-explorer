from pathlib import Path

TEST_ROOT = Path(__file__).resolve().parents[1]

DATA_DIR = TEST_ROOT / "data"

print(TEST_ROOT)
print(DATA_DIR)
