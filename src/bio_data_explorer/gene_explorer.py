import logging
from typing import Dict, Optional, List

logger = logging.getLogger(__name__)


class Gene:
    def __init__(self, identifier: str = "", description: str = "", sequence: str = ""):
        self.gene_name = identifier  # for simplicity, assume identifier is a gene_name
        self.description = description
        self.sequence = sequence


class GenesExplorer:
    def __init__(self, genes: Optional[List[Gene]] = None):
        self.gene_to_sequence: Dict[str, str] = {}
        self.sequence_to_nucleotide_counts: Dict[str, Dict[str, int]] = {}
        self.gene_name_casing_map: Dict[str, str] = {}
        self.gene_to_description: Dict[str, str] = {}
        if genes:
            for gene in genes:
                self.add_gene_sequence(gene.gene_name, gene.description, gene.sequence)
                self.add_sequence_nucleotide_counts(gene.sequence)

    def add_gene_sequence(
        self, gene_name: str = "", description: str = "", sequence: str = ""
    ) -> None:
        gene_name_lower = gene_name.lower()
        if gene_name_lower not in self.gene_to_sequence:
            self.gene_to_sequence[gene_name_lower] = sequence.upper()
            self.gene_name_casing_map[gene_name_lower] = gene_name
            self.gene_to_description[gene_name_lower] = description
        else:
            logger.warning(f"{gene_name} already exists, not adding it to genes data")

    def add_sequence_nucleotide_counts(self, sequence: str = "") -> None:
        if sequence.upper() not in self.sequence_to_nucleotide_counts:
            self.sequence_to_nucleotide_counts[sequence.upper()] = (
                self.count_nucleotides_in_sequence(sequence.upper())
            )

    def get_gene_sequence(self, gene_name: str = "") -> Optional[str]:
        if gene_name == "":
            logger.info("Must provide gene name to get its sequence")
            return None
        elif gene_name.lower() not in self.gene_to_sequence:
            logger.warning(f"Gene not found: {gene_name}")
            return None
        else:
            return self.gene_to_sequence[gene_name.lower()]

    def get_gene_name_original_casing(
        self, gene_name_case_insensitive: str = ""
    ) -> Optional[str]:
        if gene_name_case_insensitive.lower() not in self.gene_name_casing_map:
            logger.warning(f"Gene not found: {gene_name_case_insensitive}")
            return None
        else:
            return self.gene_name_casing_map[gene_name_case_insensitive.lower()]

    def count_nucleotides_in_sequence(self, sequence: str = "") -> Dict[str, int]:
        counts = {"A": 0, "C": 0, "T": 0, "G": 0}
        if sequence == "":
            logger.warning("Must provide non-empty sequence to count nucleotides")
            return counts
        else:
            for nucleotide in sequence.upper():
                if nucleotide in counts:
                    counts[nucleotide] += 1
                else:
                    logger.warning(
                        f"Unknown nucleotide present in sequence {sequence.upper()}: {nucleotide}"
                    )
            return counts

    def get_count_nucleotides(self, sequence: str = "") -> Optional[Dict[str, int]]:
        if sequence == "":
            logger.warning("Must provide sequence to count nucleotides")
            return None
        elif sequence.upper() not in self.sequence_to_nucleotide_counts:
            logger.warning("Sequence not found")
            return None
        else:
            return self.sequence_to_nucleotide_counts[sequence.upper()]

    def print_all_genes(self) -> None:
        print(f"There are {len(self.gene_name_casing_map)} genes loaded: ")
        number = 1
        for k, v in self.gene_name_casing_map.items():
            print(f"{number}. {v} | {self.gene_to_description[k]}")
            number += 1

    def pretty_print_count_nucleotides(
        self, gene_name: str = "", nucleotide_counts: Optional[Dict[str, int]] = None
    ) -> None:
        if nucleotide_counts is None:
            logger.warning("Must provide nucleotide count map to print")
        else:
            parts = [f"{k}={v}" for k, v in nucleotide_counts.items()]
            pretty_printed = " ".join(parts)
            print(f"{self.get_gene_name_original_casing(gene_name)}: {pretty_printed}")
