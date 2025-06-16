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
        with gzip.open(
            f"{data_directory_path}/{fastq_file_name}", "rt", encoding="utf-8"
        ) as fastq_file:
            genes: List[Gene] = []
            try:
                # https://biopython.org/wiki/SeqIO
                for record in SeqIO.parse(fastq_file, FILE_TYPE_FASTQ):  # type: ignore
                    genes.append(Gene(record.id, record.description, record.seq, record.letter_annotations["phred_quality"]))  # type: ignore
            except Exception as e:
                logger.error(f"Could not parse fastq file: {e}")
            return genes
    except FileNotFoundError as e:
        logger.error(f"Could not open {fastq_file_name}: {e.strerror}")
        raise FastqParsingError(e)
