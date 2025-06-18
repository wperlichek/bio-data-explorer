import logging
from pathlib import Path
from typing import List
from cyvcf2 import VCF  # type: ignore

data_directory_path = Path(__file__).resolve().parent.parent.parent / "data"

logger = logging.getLogger(__name__)

MIN_QUAL_SCORE = 30.0
MIN_INFO_DP = 20


def show_low_confidence_variants(vcf_file_name: str = "") -> List[str]:
    try:
        vcf = VCF(f"{data_directory_path}/{vcf_file_name}")  # type: ignore
    except Exception as e:
        logger.error(f"Could not open {vcf_file_name}: {e}")

    low_confidence_variants: List[str] = []

    for variant in vcf:  # type: ignore
        # https://brentp.github.io/cyvcf2/
        # https://brentp.github.io/cyvcf2/docstrings.html#api
        filter = str(variant.FILTER)  # type: ignore
        qual = float(variant.QUAL)  # type: ignore
        info_dp = int(variant.INFO.get("DP"))  # type: ignore
        variant_identifer = f"{variant.CHROM}:{variant.POS}_{variant.REF}>{" ".join(variant.ALT)}"  # type: ignore
        if is_low_confidence_variant(filter, qual, info_dp, variant_identifer):  # type: ignore
            logger.warning(f"Filtering out variant {variant_identifer}")
            low_confidence_variants.append(variant_identifer)

    return low_confidence_variants


def is_low_confidence_variant(
    filter: str = "", qual: float = 0.0, info_dp: int = 0, variant_identifer: str = ""
) -> bool:
    if filter != "None":
        logger.warning(
            f"Variant {variant_identifer} did not PASS in filter, filter value: {filter}"
        )
        return True
    if qual < MIN_QUAL_SCORE:
        logger.warning(
            f"Variant {variant_identifer} qual score {qual} does not meet min qual threshold {MIN_QUAL_SCORE}"
        )
        return True
    if info_dp < MIN_INFO_DP:
        logger.warning(
            f"Variant {variant_identifer} depth (info_dp) {info_dp} does not meet min depth threshold {MIN_INFO_DP}"
        )
        return True
    return False
