from typing import Optional, Dict
from Bio import SeqIO
from pathlib import Path
import logging, gzip

data_directory_path = Path(__file__).resolve().parent.parent.parent / "data"

FILE_TYPE_FASTQ = "fastq"

logger = logging.getLogger(__name__)

class FastqParsingError(Exception):
    pass

def parse_fastq_file(fastq_file_name: str = "") -> Optional[Dict[str, str]]:
    try:
        with gzip.open(fastq_file_name) as fastq_file:
            try:
                # https://biopython.org/wiki/SeqIO
                for record in SeqIO.parse(fastq_file, FILE_TYPE_FASTQ): # type: ignore
                    logger.info(f"Parsed record: {record.name}") # type: ignore
            except Exception as e:
                logger.error(f"Could not parse fastq file: {e}")
            sequence_to_qualities: Dict[str, str] = {}
            return sequence_to_qualities
    except FileNotFoundError as e:
        logger.error(f"Could not open {fastq_file_name}: {e.strerror}")
        raise FastqParsingError(e)
