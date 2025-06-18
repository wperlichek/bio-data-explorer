import logging
from pathlib import Path
import pysam

data_directory_path = Path(__file__).resolve().parent.parent.parent / "data"

logger = logging.getLogger(__name__)


class SamBamParsingError(Exception):
    pass


def parse_sam_file(sam_file_name: str = ""):
    try:
        sam_file = pysam.AlignmentFile(f"{data_directory_path}/{sam_file_name}", "r")
    except Exception as e:
        logger.error(f"Could not open {sam_file_name}: {e}")
        raise SamBamParsingError(e)

    for read in sam_file.fetch("chr1", 100, 120):
        print(read)
    sam_file.close()
    return


def parse_bam_file():
    return
