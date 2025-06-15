from Bio.Blast import NCBIWWW
from enum import Enum
from typing import Optional, cast
from io import StringIO
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
) -> str:
    if sequence == "":
        logger.warning("Must provide sequence to make a BLAST call")
    if program is None:
        program = BlastProgram.BLASTN
    if database is None:
        database = BlastDatabase.NT
    
    return cast(StringIO, NCBIWWW.qblast(program.value, database.value, sequence))  # type: ignore[reportUnknownMemberType]