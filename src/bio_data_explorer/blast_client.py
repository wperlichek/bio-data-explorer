from Bio.Blast import NCBIWWW, NCBIXML
from enum import Enum
from typing import Optional
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
) -> Optional[int]:
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
        blast_records = NCBIXML.parse(result_handle)  # type: ignore
    except Exception as e:
        logger.error(f"Problem parsing result handle {result_handle.getvalue}: {e}")
        return None

    count_hits = 0
    for record in blast_records:  # type: ignore
        for alignment in record.alignments:  # type: ignore
            count_hits += 1
            print(f"  Hit: {alignment.title}")  # type: ignore

    return count_hits
