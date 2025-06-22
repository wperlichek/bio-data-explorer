from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[3]

DATA_DIR = PROJECT_ROOT / "data"
FASTA_PATH = DATA_DIR / "fasta"
FASTQ_PATH = DATA_DIR / "fastq"
SAM_BAM_PATH = DATA_DIR / "sam-bam"
VCF_PATH = DATA_DIR / "vcf"
