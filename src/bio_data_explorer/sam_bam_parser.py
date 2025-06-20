import logging, os
from pathlib import Path
from typing import Dict
import pysam
from pysam import AlignmentFile

logger = logging.getLogger(__name__)

MIN_MAPPING_QUALITY_THRESHOLD = 30


class SamBamParsingError(Exception):
    pass


def open_alignment_file(sam_or_bam_file: str = "") -> AlignmentFile:
    """
    This function also creates a .bai index
    if sam_or_bam_file is .bam and an index doesn't already exist
    """
    mode = "r" if Path(sam_or_bam_file).suffix == ".sam" else "rb"
    try:
        alignment_file = AlignmentFile(sam_or_bam_file, mode)
        if mode == "rb" and not alignment_file.has_index():
            create_index_for_bam_file(sam_or_bam_file)
            alignment_file = AlignmentFile(sam_or_bam_file, mode)
    except Exception as e:
        logger.error(f"Could not open {sam_or_bam_file}: {e}")
        raise SamBamParsingError(e)
    return alignment_file


def create_index_for_bam_file(bam_file: str = ""):
    temp_sorted_bam_file = f"{bam_file}.tmp_sorted.bam"
    pysam.sort("-o", temp_sorted_bam_file, bam_file)
    os.replace(temp_sorted_bam_file, bam_file)
    pysam.index(bam_file)
    logger.info(f"Created index {bam_file}.bai")


def get_read_alignment_stats_summary(alignment_file: AlignmentFile):
    reads_count = 0
    mapped_reads = 0
    unmapped_reads = 0
    low_quality_reads = 0
    duplicate_reads = 0
    supplimentary_reads = 0

    for read in alignment_file:
        if read.is_mapped:
            mapped_reads += 1
            if read.mapping_quality < MIN_MAPPING_QUALITY_THRESHOLD:
                low_quality_reads += 1
                logger.warning(
                    f"Read {read.query_name} map quality of {read.mapping_quality} is less than threshold of {MIN_MAPPING_QUALITY_THRESHOLD}"
                )
            if read.is_duplicate:
                duplicate_reads += 1
            if read.is_supplementary:
                supplimentary_reads += 1
        else:
            unmapped_reads += 1

        reads_count += 1

    mapping_rate_percentage = int(round(mapped_reads / reads_count, 2) * 100)

    alignment_stats: Dict[str, int] = {}
    alignment_stats["reads_count"] = reads_count
    alignment_stats["mapped_reads"] = mapped_reads
    alignment_stats["unmapped_reads"] = unmapped_reads
    alignment_stats["mapping_rate_percentage"] = mapping_rate_percentage
    alignment_stats["low_quality_reads"] = low_quality_reads
    alignment_stats["duplicate_reads"] = duplicate_reads
    alignment_stats["supplimentary_reads"] = supplimentary_reads

    return alignment_stats
