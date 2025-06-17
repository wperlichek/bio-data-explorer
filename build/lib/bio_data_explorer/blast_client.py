from Bio.Blast import NCBIWWW, NCBIXML
from enum import Enum
from typing import List, Optional, Any
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
) -> Optional[List[Any]]:
    if sequence == "":
        logger.warning("Must provide sequence to make a BLAST call")
        return None
    if program is None:
        program = BlastProgram.BLASTN
    if database is None:
        database = BlastDatabase.NT
    try:
        result_handle = NCBIWWW.qblast(program.value, database.value, sequence)  # type: ignore[reportUnknownMemberType]
    except Exception as e:
        logger.error(f"Problem during external call to {NCBIWWW.NCBI_BLAST_URL}: {e}")
        return None

    try:
        # https://biopython.org/docs/1.76/api/Bio.Blast.Record.html
        blast_records = NCBIXML.parse(result_handle)  # type: ignore
    except Exception as e:
        logger.error(f"Problem parsing result handle {result_handle.getvalue}: {e}")
        return None

    return blast_records  # type: ignore
