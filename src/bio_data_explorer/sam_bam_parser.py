import logging
from pathlib import Path
import pysam

data_directory_path = Path(__file__).resolve().parent.parent.parent / "data"

logger = logging.getLogger(__name__)


class SamBamParsingError(Exception):
    pass


def parse_sam_bam_file(sam_file_name: str = ""):
    mode = "r" if Path(sam_file_name).suffix == ".sam" else "rb"
    try:
        sam_file = pysam.AlignmentFile(f"{data_directory_path}/{sam_file_name}", mode)
    except Exception as e:
        logger.error(f"Could not open {sam_file_name}: {e}")
        raise SamBamParsingError(e)

    for read in sam_file:
        print(read)

    sam_file.close()
    return
