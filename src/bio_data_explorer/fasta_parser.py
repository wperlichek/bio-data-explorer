import logging, gzip
from typing import List
from .gene_explorer import Gene
from .config import path_config

FASTA_SEQUENCE_CHARS_DNA_ONLY = {"ACGT"}

logger = logging.getLogger(__name__)


class FastaParsingError(Exception):
    pass


def parse_fasta_file(fasta_file_name: str = "") -> List[Gene]:
    genes: List[Gene] = []
    try:
        with gzip.open(
            f"{path_config.DATA_DIR}/{fasta_file_name}",
            "rt",
            encoding="utf-8",
        ) as File:
            identifier = ""
            description = ""
            sequence = ""
            found_format_problem = False
            for line in File:
                stripped_line = line.strip()
                if not stripped_line:
                    pass
                else:
                    if _line_is_formatted_correctly(stripped_line):
                        if stripped_line[0] == ">":
                            if sequence:
                                if not found_format_problem:
                                    genes.append(
                                        Gene(identifier, description, sequence.upper())
                                    )
                                else:
                                    found_format_problem = False
                            identifier_and_description = (
                                stripped_line[1::].strip().split(" ")
                            )
                            identifier = identifier_and_description[0]
                            description = " ".join(identifier_and_description[1::])
                            sequence = ""
                        else:
                            sequence += stripped_line
                    else:
                        logger.warning(
                            f"{line} is not formatted correctly, skipping this entire FASTA entry"
                        )
                        if sequence and not found_format_problem:
                            genes.append(
                                Gene(identifier, description, sequence.upper())
                            )
                        found_format_problem = True
            if not found_format_problem:
                genes.append(Gene(identifier, description, sequence.upper()))
        return genes
    except FileNotFoundError as e:
        logger.error(f"Could not open {fasta_file_name}: {e.strerror}")
        raise FastaParsingError(e)


def _line_is_formatted_correctly(line: str = "") -> bool:
    if line[0] == ">":
        there_are_contents = len(line) > 1
        if there_are_contents:
            return True
        else:
            logger.warning(f"Line {line} is missing identifier")
            return False
    else:
        for ch in line:
            if ch not in FASTA_SEQUENCE_CHARS_DNA_ONLY:
                logger.warning(f"Found non-DNA FASTA character {ch} in {line}")
                return False
        return True
