import logging, sys, gzip
from typing import Dict, Optional, List

GENES_FILE = "sample_genes.fasta.gz"

FASTA_SEQUENCE_CHARS_DNA = set("ACGT")

logging.basicConfig(
    level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s"
)


class GenesFileParsingError(Exception):
    pass


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
            logging.warning(f"{gene_name} already exists, not adding it to genes data")

    def add_sequence_nucleotide_counts(self, sequence: str = "") -> None:
        if sequence.upper() not in self.sequence_to_nucleotide_counts:
            self.sequence_to_nucleotide_counts[sequence.upper()] = (
                self.count_nucleotides_in_sequence(sequence.upper())
            )

    def get_gene_sequence(self, gene_name: str = "") -> Optional[str]:
        if gene_name == "":
            logging.info("Must provide gene name to get its sequence")
            return None
        elif gene_name.lower() not in self.gene_to_sequence:
            logging.warning(f"Gene not found: {gene_name}")
            return None
        else:
            return self.gene_to_sequence[gene_name.lower()]

    def get_gene_name_original_casing(
        self, gene_name_case_insensitive: str = ""
    ) -> Optional[str]:
        if gene_name_case_insensitive.lower() not in self.gene_name_casing_map:
            logging.warning(f"Gene not found: {gene_name_case_insensitive}")
            return None
        else:
            return self.gene_name_casing_map[gene_name_case_insensitive.lower()]

    def count_nucleotides_in_sequence(self, sequence: str = "") -> Dict[str, int]:
        counts = {"A": 0, "C": 0, "T": 0, "G": 0}
        if sequence == "":
            logging.warning("Must provide non-empty sequence to count nucleotides")
            return counts
        else:
            for nucleotide in sequence.upper():
                if nucleotide in counts:
                    counts[nucleotide] += 1
                else:
                    logging.warning(
                        f"Unknown nucleotide present in sequence {sequence.upper()}: {nucleotide}"
                    )
            return counts

    def get_count_nucleotides(self, sequence: str = "") -> Optional[Dict[str, int]]:
        if sequence == "":
            logging.warning("Must provide sequence to count nucleotides")
            return None
        elif sequence.upper() not in self.sequence_to_nucleotide_counts:
            logging.warning("Sequence not found")
            return None
        else:
            return self.sequence_to_nucleotide_counts[sequence.upper()]

    def print_all_genes(self) -> None:
        # TODO :: print the description too
        print(f"There are {len(self.gene_name_casing_map)} genes loaded: ")
        number = 1
        for k, v in self.gene_name_casing_map.items():
            print(f"{number}. {v} | {self.gene_to_description[k]}")
            number += 1

    def pretty_print_count_nucleotides(
        self, gene_name: str = "", nucleotide_counts: Optional[Dict[str, int]] = None
    ) -> None:
        if nucleotide_counts is None:
            logging.warning("Must provide nucleotide count map to print")
        else:
            parts = [f"{k}={v}" for k, v in nucleotide_counts.items()]
            pretty_printed = " ".join(parts)
            print(f"{self.get_gene_name_original_casing(gene_name)}: {pretty_printed}")


def parse_genes_data(genes_file: str = "") -> List[Gene]:
    genes: List[Gene] = []
    try:
        with gzip.open(genes_file, "rt", encoding="utf-8") as File:
            identifier = ""
            description = ""
            sequence = ""
            found_format_problem = False
            for line in File:
                stripped_line = line.strip()
                if not stripped_line:
                    pass
                else:
                    if line_is_formatted_correctly(stripped_line):
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
                        logging.warning(
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
        logging.error(f"Could not open {genes_file}: {e.strerror}")
        raise GenesFileParsingError(e)


def line_is_formatted_correctly(line: str = "") -> bool:
    if line[0] == ">":
        there_are_contents = len(line) > 1
        if there_are_contents:
            return True
        else:
            logging.warning(f"Line {line} is missing identifier")
            return False
    else:
        for ch in line:
            if ch not in FASTA_SEQUENCE_CHARS_DNA:
                logging.warning(f"Found non-FASTA character {ch} in {line}")
                return False
        return True


def cli_app() -> None:

    logging.info("Starting app")

    try:
        genes = parse_genes_data(GENES_FILE)
    except GenesFileParsingError as e:
        logging.critical(f"Application can't start due to {e}, exiting application")
        sys.exit(1)

    genes_explorer = GenesExplorer(genes)

    while True:
        print("** Gene Sequence Explorer **")
        print("1. List all genes")
        print("2. View sequence of gene")
        print("3. Count nucleotides of gene")
        print("4. Exit application")

        menu_choice = input("Enter choice (1-4): ").strip()

        if menu_choice == "1":
            genes_explorer.print_all_genes()
        elif menu_choice == "2":
            gene_name = input("Enter gene name: ").strip().lower()
            result = genes_explorer.get_gene_sequence(gene_name)
            if result:
                print(
                    f"{genes_explorer.get_gene_name_original_casing(gene_name)}: {result}"
                )
        elif menu_choice == "3":
            gene_name = input("Enter gene name: ").strip().lower()
            sequence = genes_explorer.get_gene_sequence(gene_name)
            if sequence:
                result = genes_explorer.get_count_nucleotides(sequence)
                if result:
                    genes_explorer.pretty_print_count_nucleotides(gene_name, result)
        elif menu_choice == "4":
            logging.info("Exiting app")
            break
        else:
            print("Invalid choice")


if __name__ == "__main__":
    cli_app()
