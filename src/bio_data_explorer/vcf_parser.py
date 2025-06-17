import logging
from pathlib import Path
import cyvcf2  # type: ignore
from cyvcf2 import VCF  # type: ignore

data_directory_path = Path(__file__).resolve().parent.parent.parent / "data"

logger = logging.getLogger(__name__)


def parse_vcf_file(vcf_file_name: str = ""):
    print("cyvcf2 version:", cyvcf2.__version__)
    try:
        for variant in VCF(f"{data_directory_path}/{vcf_file_name}"):  # type: ignore
            print(variant.REF)  # type: ignore
    except Exception as e:
        logger.error(f"Could not open {vcf_file_name}: {e}")
