from Bio.Blast import NCBIWWW
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
) -> None:
    if sequence == "":
        logger.warning("Must provide sequence to make a BLAST call")
    if program is None:
        program = BlastProgram.BLASTN
    if database is None:
        database = BlastDatabase.NT

    result = cast(StringIO, NCBIWWW.qblast(program.value, database.value, sequence))  # type: ignore[reportUnknownMemberType]
