from typing import List
from Bio import SeqIO
from pathlib import Path
import logging, gzip
from .gene_explorer import Gene

data_directory_path = Path(__file__).resolve().parent.parent.parent / "data"

FILE_TYPE_FASTQ = "fastq"

UNKNOWN_BASES_THRESHOLD_PERCENTAGE_TO_OMIT_READ = 5

MIN_PHRED_QUALITY_TO_KEEP_WHILE_TRIMMING_END = 20

MIN_SEQUENCE_READ_LEN = 50

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
                    phred_qualities: List[int] = record.letter_annotations["phred_quality"]  # type: ignore
                    if valid_sequence(seq, phred_qualities):  # type: ignore
                        genes.append(Gene(record.id, record.description, seq, phred_qualities))  # type: ignore
            except Exception as e:
                logger.error(f"Could not parse fastq file: {e}")
            return genes
    except FileNotFoundError as e:
        logger.error(f"Could not open {fastq_file_name}: {e.strerror}")
        raise FastqParsingError(e)


def valid_sequence(sequence: str = "", phred_qualities: List[int] = []) -> bool:

    if len(sequence) < MIN_SEQUENCE_READ_LEN:
        logger.warning(
            f"Discarding sequence {sequence} due to it not being >= min sequence length {MIN_SEQUENCE_READ_LEN}"
        )
        return False

    if should_discard_read_due_to_high_unknown_base_count(sequence):
        logger.warning(
            f"Discarding sequence {sequence} due to high % of unknown bases, max threshold is {UNKNOWN_BASES_THRESHOLD_PERCENTAGE_TO_OMIT_READ}%"
        )
        return False

    trimmed_sequence = quality_trim_sequence_end(sequence, phred_qualities)

    if len(trimmed_sequence) < MIN_SEQUENCE_READ_LEN:
        logger.warning(
            f"Discarding trimmed sequence {sequence} due to it not being >= min sequence length {MIN_SEQUENCE_READ_LEN}"
        )
        return False

    return True


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


def quality_trim_sequence_end(
    sequence: str = "", phred_qualities: List[int] = []
) -> str:
    first_bad_read_idx = -1
    for idx in range(len(phred_qualities)):
        if phred_qualities[idx] < MIN_PHRED_QUALITY_TO_KEEP_WHILE_TRIMMING_END:
            first_bad_read_idx = idx
            break
    return sequence if first_bad_read_idx == -1 else sequence[0:first_bad_read_idx]
