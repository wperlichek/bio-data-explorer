from pathlib import Path

data_directory_path = Path(__file__).resolve().parent.parent.parent / "data"


def parse_vcf_file(vcf_file_name: str = ""):
    print(vcf_file_name)
