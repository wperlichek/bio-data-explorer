from Bio.Blast import NCBIWWW, NCBIXML, Record
from enum import Enum
from typing import Iterator, Optional
from unittest.mock import patch
import logging

logger = logging.getLogger(__name__)


class BlastProgram(Enum):
    BLASTN = "blastn"


class BlastDatabase(Enum):
    NT = "nt"


def make_blast_call(
    program: Optional[BlastProgram] = None,
    database: Optional[BlastDatabase] = None,
    sequence: str = "",
) -> Optional[Iterator[Record]]:
    if sequence == "":
        logger.warning("Must provide sequence to make a BLAST call")
        return None
    if program is None:
        program = BlastProgram.BLASTN
    if database is None:
        database = BlastDatabase.NT

    try:
        result_handle = NCBIWWW.qblast(program.value, database.value, sequence)
    except Exception as e:
        logger.error(f"Problem during external call to {NCBIWWW.NCBI_BLAST_URL}: {e}")
        return None

    try:
        # See https://biopython.org/docs/1.76/api/Bio.Blast.Record.html
        blast_record = NCBIXML.parse(result_handle)
    except Exception as e:
        logger.error(f"Problem parsing result handle {result_handle.getvalue}: {e}")
        return None

    return blast_record
