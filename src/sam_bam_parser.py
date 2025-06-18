import logging
from pathlib import Path
import pysam

data_directory_path = Path(__file__).resolve().parent.parent.parent / "data"

logger = logging.getLogger(__name__)


def parse_sam_file(file_name: str = ""):
    sam_file = pysam.AlignmentFile(f"{data_directory_path}/{file_name}", "r")
    for read in sam_file.fetch("chr1", 100, 120):
        print(read)
    sam_file.close()
    return


def parse_bam_file():
    return
