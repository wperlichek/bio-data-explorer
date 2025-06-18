import logging
from pathlib import Path
import pysam

data_directory_path = Path(__file__).resolve().parent.parent.parent / "data"

logger = logging.getLogger(__name__)


class SamBamParsingError(Exception):
    pass


def parse_sam_or_bam_file(sam_or_sam_file_name: str = ""):
    mode = "r" if Path(sam_or_sam_file_name).suffix == ".sam" else "rb"
    try:
        sam_file = pysam.AlignmentFile(
            f"{data_directory_path}/{sam_or_sam_file_name}", mode
        )
    except Exception as e:
        logger.error(f"Could not open {sam_or_sam_file_name}: {e}")
        raise SamBamParsingError(e)

    for read in sam_file:
        print(read)

    sam_file.close()


def create_bai_from_bam_file(bam_file_name: str = ""):
    bai_index_file = f"{data_directory_path}/sorted_{bam_file_name}"
    logger.info(f"Creating index for {bam_file_name} at {bai_index_file}")
    pysam.sort(
        "-o",
        bai_index_file,
        f"{data_directory_path}/{bam_file_name}",
    )
    pysam.index(bai_index_file)
