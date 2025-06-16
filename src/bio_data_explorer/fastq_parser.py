from typing import List
from Bio import SeqIO
from pathlib import Path
import logging, gzip
from .gene_explorer import Gene

data_directory_path = Path(__file__).resolve().parent.parent.parent / "data"

FILE_TYPE_FASTQ = "fastq"

UNKNOWN_BASES_THRESHOLD_PERCENTAGE_TO_OMIT_READ = 5

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
                    seq = str(record.seq)  # type: ignore
                    if not should_discard_read_due_to_high_unknown_base_count(seq):
                        genes.append(Gene(record.id, record.description, seq, record.letter_annotations["phred_quality"]))  # type: ignore
                    else:
                        logger.warning(
                            f"Discarding sequence {seq} due to high % of unknown bases, max threshold is {UNKNOWN_BASES_THRESHOLD_PERCENTAGE_TO_OMIT_READ}%"
                        )
            except Exception as e:
                logger.error(f"Could not parse fastq file: {e}")
            return genes
    except FileNotFoundError as e:
        logger.error(f"Could not open {fastq_file_name}: {e.strerror}")
        raise FastqParsingError(e)


def should_discard_read_due_to_high_unknown_base_count(sequence: str = "") -> bool:
    if sequence == "":
        logger.warning(
            "Must provide a sequence to check its unknown base count, returning False without analysis"
        )
        return False
    unknown_base_count = 0
    for ch in sequence.upper():
        if ch == "N":
            unknown_base_count += 1
    percentage_unknowns = round(unknown_base_count / len(sequence), 3) * 100
    logger.info(
        f"Percentage of unknown bases in sequence {sequence}: {percentage_unknowns}%"
    )
    return percentage_unknowns >= UNKNOWN_BASES_THRESHOLD_PERCENTAGE_TO_OMIT_READ
