from Bio.Blast import NCBIWWW, NCBIXML
from enum import Enum
from typing import Optional
import logging

logger = logging.getLogger(__name__)


class BlastProgram(Enum):
    BLASTN = "blastn"


class BlastDatabase(Enum):
    NT = "nt"
    CORE_NOT = "Core_nt"


def make_blast_call(
    program: Optional[BlastProgram] = None,
    database: Optional[BlastDatabase] = None,
    sequence: str = "",
) -> Optional[str]:
    if sequence == "":
        logger.warning("Must provide sequence to make a BLAST call")
        return None
    if program is None:
        program = BlastProgram.BLASTN
    if database is None:
        database = BlastDatabase.NT
    result_handle = NCBIWWW.qblast(program.value, database.value, sequence)  # type: ignore[reportUnknownMemberType]
    blast_records = NCBIXML.parse(result_handle) # type: ignore
    for record in blast_records: # type: ignore
        for alignment in record.alignments: # type: ignore
            print(f"  Hit: {alignment.title}") # type: ignore
    return ""