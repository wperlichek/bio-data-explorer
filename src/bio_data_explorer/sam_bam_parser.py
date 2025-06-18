import logging
from pathlib import Path
from typing import Dict
import pysam
from pysam import AlignmentFile

data_directory_path = Path(__file__).resolve().parent.parent.parent / "data"

logger = logging.getLogger(__name__)


class SamBamParsingError(Exception):
    pass


def parse_sam_or_bam_file(sam_or_bam_file_name: str = "") -> AlignmentFile:
    logger.info(f"Parsing {sam_or_bam_file_name}...")
    mode = "r" if Path(sam_or_bam_file_name).suffix == ".sam" else "rb"
    try:
        alignment_file = AlignmentFile(
            f"{data_directory_path}/{sam_or_bam_file_name}", mode
        )
        logger.info(f"Succesfully parsed {sam_or_bam_file_name}")
        if mode == "rb":  # a .bam file
            create_bai_from_bam_file(sam_or_bam_file_name)
    except Exception as e:
        logger.error(f"Could not open {sam_or_bam_file_name}: {e}")
        raise SamBamParsingError(e)
    return alignment_file


def create_bai_from_bam_file(bam_file_name: str = ""):
    bai_index_file = f"{data_directory_path}/sorted_{bam_file_name}"
    logger.info(f"Creating index for {bam_file_name} at {bai_index_file}")
    pysam.sort(
        "-o",
        bai_index_file,
        f"{data_directory_path}/{bam_file_name}",
    )
    pysam.index(bai_index_file)
    logger.info(f"Created index {bai_index_file}")


def get_read_alignment_stats_summary(alignment_file: AlignmentFile):
    reads_count = 0
    mapped_reads = 0
    unmapped_reads = 0

    for read in alignment_file:
        if read.is_mapped:
            mapped_reads += 1
        else:
            unmapped_reads += 1

        reads_count += 1

    mapping_rate_percentage = int(round(mapped_reads / reads_count, 2) * 100)

    alignment_stats: Dict[str, int] = {}

    alignment_stats["reads_count"] = reads_count
    alignment_stats["mapped_reads"] = mapped_reads
    alignment_stats["unmapped_reads"] = unmapped_reads
    alignment_stats["mapping_rate_percentage"] = mapping_rate_percentage

    return alignment_stats
