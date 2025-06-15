from Bio import Blast
from enum import Enum
from typing import Optional


class BlastProgram(Enum):
    BLASTN = "blastn"


class BlastDatabase(Enum):
    NT = "nt"


def make_blast_call(
    program: Optional[BlastProgram] = None,
    database: Optional[BlastDatabase] = None,
    sequence: str = "",
) -> None:
    
    print("hello")