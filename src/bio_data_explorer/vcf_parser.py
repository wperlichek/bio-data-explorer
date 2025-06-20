import logging
from typing import List
from cyvcf2 import VCF

logger = logging.getLogger(__name__)

MIN_QUAL_SCORE = 30.0
MIN_INFO_DP = 20
PASSING_FILTER_VAL = "None"


def show_low_confidence_variants(vcf_file: str = "") -> List[str]:
    """
    Returns a list of low confidence variants given vcf file
    Uses helper method _is_low_confidence_variant to determine confidence
    Creates unique variant identifier to help with logging
    """

    try:
        vcf = VCF(vcf_file)
    except Exception as e:
        logger.error(f"Could not open {vcf_file}: {e}")

    low_confidence_variants: List[str] = []

    for variant in vcf:
        # See https://brentp.github.io/cyvcf2/
        filter = str(variant.FILTER)
        qual = float(variant.QUAL)
        info_dp = int(variant.INFO.get("DP"))
        alts = " ".join(variant.ALT)
        variant_identifier = f"{variant.CHROM}:{variant.POS}_{variant.REF}>{alts}"
        if _is_low_confidence_variant(filter, qual, info_dp, variant_identifier):
            logger.warning(f"Found low confidence variant: {variant_identifier}")
            low_confidence_variants.append(variant_identifier)

    return low_confidence_variants


def _is_low_confidence_variant(
    filter: str = "", qual: float = 0.0, info_dp: int = 0, variant_identifier: str = ""
) -> bool:
    if filter != PASSING_FILTER_VAL:
        logger.warning(
            f"Variant {variant_identifier} did not PASS in filter, filter value: {filter}"
        )
        return True
    if qual < MIN_QUAL_SCORE:
        logger.warning(
            f"Variant {variant_identifier} qual score {qual} does not meet min qual threshold {MIN_QUAL_SCORE}"
        )
        return True
    if info_dp < MIN_INFO_DP:
        logger.warning(
            f"Variant {variant_identifier} depth (info_dp) {info_dp} does not meet min depth threshold {MIN_INFO_DP}"
        )
        return True
    return False
