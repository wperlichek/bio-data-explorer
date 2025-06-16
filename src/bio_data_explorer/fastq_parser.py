from typing import List
from Bio import SeqIO
from pathlib import Path
import logging, gzip
from .gene_explorer import Gene

data_directory_path = Path(__file__).resolve().parent.parent.parent / "data"

FILE_TYPE_FASTQ = "fastq"

logger = logging.getLogger(__name__)


class FastqParsingError(Exception):
    pass


def parse_fastq_file(fastq_file_name: str = "") -> List[Gene]:
    try:
        with gzip.open(fastq_file_name) as fastq_file:
            try:
                # https://biopython.org/wiki/SeqIO
                for record in SeqIO.parse(fastq_file, FILE_TYPE_FASTQ):  # type: ignore
                    logger.info(f"Parsed record: {record.name}")  # type: ignore
            except Exception as e:
                logger.error(f"Could not parse fastq file: {e}")
            genes: List[Gene] = []
            return genes
    except FileNotFoundError as e:
        logger.error(f"Could not open {fastq_file_name}: {e.strerror}")
        raise FastqParsingError(e)
