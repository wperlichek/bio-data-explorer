from typing import List
from Bio import SeqIO
import logging, gzip
from .gene_explorer import Gene

UNKNOWN_BASES_THRESHOLD_PERCENTAGE_TO_OMIT_READ = 5

MIN_PHRED_QUALITY_TO_KEEP_WHILE_TRIMMING_END = 20

MIN_SEQUENCE_READ_LEN = 50

logger = logging.getLogger(__name__)


class FastqParsingError(Exception):
    pass


def parse_fastq_file(fastq_file_name: str = "") -> List[Gene]:
    """
    Uses Biopython to parse the fastq file
    Excludes low-quality reads based on sequence length and base composition
    """

    try:
        with gzip.open(fastq_file_name, "rt", encoding="utf-8") as fastq_file:
            try:
                genes: List[Gene] = []
                # See https://biopython.org/wiki/SeqIO for seq_record fields
                for seq_record in SeqIO.parse(fastq_file, "fastq"):
                    seq = str(seq_record.seq)
                    phred_qualities: List[int] = seq_record.letter_annotations[
                        "phred_quality"
                    ]
                    if _valid_sequence(seq, phred_qualities):
                        trimmed_seq = _quality_trim_sequence_end(seq)
                        genes.append(
                            Gene(
                                seq_record.id,
                                seq_record.description,
                                trimmed_seq,
                                phred_qualities[0 : len(trimmed_seq)],
                            )
                        )
                return genes
            except Exception as e:
                logger.error(f"Could not parse {fastq_file_name}: {e}")
                raise (FastqParsingError(e))
    except FileNotFoundError as e:
        logger.error(f"Could not open {fastq_file_name}: {e.strerror}")
        raise FastqParsingError(e)


def _valid_sequence(sequence: str = "", phred_qualities: List[int] = []) -> bool:
    if len(sequence) < MIN_SEQUENCE_READ_LEN:
        logger.warning(
            f"Discarding sequence {sequence} because its length is < min sequence length {MIN_SEQUENCE_READ_LEN}"
        )
        return False

    if _should_discard_read_due_to_high_unknown_base_count(sequence):
        logger.warning(
            f"Discarding sequence {sequence} because high % of unknown bases, max threshold is {UNKNOWN_BASES_THRESHOLD_PERCENTAGE_TO_OMIT_READ}%"
        )
        return False

    trimmed_sequence = _quality_trim_sequence_end(sequence, phred_qualities)

    if len(trimmed_sequence) < MIN_SEQUENCE_READ_LEN:
        logger.warning(
            f"Discarding trimmed sequence {sequence} because its length is < min sequence length {MIN_SEQUENCE_READ_LEN}"
        )
        return False

    return True


def _should_discard_read_due_to_high_unknown_base_count(sequence: str = "") -> bool:
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
    return percentage_unknowns >= UNKNOWN_BASES_THRESHOLD_PERCENTAGE_TO_OMIT_READ


def _quality_trim_sequence_end(
    sequence: str = "", phred_qualities: List[int] = []
) -> str:
    first_bad_read_idx = -1
    for idx in range(len(phred_qualities)):
        if phred_qualities[idx] < MIN_PHRED_QUALITY_TO_KEEP_WHILE_TRIMMING_END:
            first_bad_read_idx = idx
            break
    final_sequence = (
        sequence if first_bad_read_idx == -1 else sequence[0:first_bad_read_idx]
    )
    if len(final_sequence) < len(sequence):
        logger.warning(f"Trimmed sequence {sequence} down to {final_sequence}")
    return final_sequence
