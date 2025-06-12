import sys
from typing import Dict, Optional, List

GENES_FILE = "genes.txt"


class Gene:
    def __init__(self, gene_name, sequence):
        self.gene_name = gene_name
        self.sequence = sequence


class GenesData:

    def __init__(self):
        self.gene_to_sequence: Dict[str, str] = {}
        self.sequence_to_nucleotide_counts: Dict[str, Dict[str, int]] = {}

    def add_gene_to_sequence(self, gene_name: str = "", sequence: str = "") -> None:
        if gene_name not in self.gene_to_sequence:
            self.gene_to_sequence[gene_name] = sequence
        else:
            print(f"Warning: {gene_name} already exists, not adding it to genes data")

    def add_sequence_to_nucleotide_counts(self, sequence: str = "") -> None:
        if sequence not in self.sequence_to_nucleotide_counts:
            self.sequence_to_nucleotide_counts[sequence] = (
                self.count_nucleotides_in_sequence(sequence)
            )

    def get_gene_sequence(self, gene_name: str = "") -> Optional[str]:
        if gene_name == "":
            print("Must provide gene name to get its sequence")
            return None
        elif gene_name not in self.gene_to_sequence:
            print("Gene not found")
            return None
        else:
            return self.gene_to_sequence[gene_name]

    def count_nucleotides_in_sequence(self, sequence: str = "") -> Dict[str, int]:
        if sequence == "":
            print("Must provide non-empty sequence to count nucleotides")
        else:
            counts = {"A": 0, "C": 0, "T": 0, "G": 0}
            for nucleotide in sequence:
                if nucleotide in counts:
                    counts[nucleotide] += 1
                else:
                    print(
                        "Warning: Unknown nucleotide present in sequence "
                        + sequence
                        + ": "
                        + nucleotide
                    )
                    pass
            return counts

    def get_count_nucleotides(self, sequence: str = "") -> Optional[Dict[str, int]]:
        if sequence == "":
            print("Must provide sequence to count nucleotides")
            return None
        elif sequence not in self.sequence_to_nucleotide_counts:
            print("Sequence not found")
            return None
        else:
            return self.sequence_to_nucleotide_counts[sequence]

    def print_all_genes(self) -> None:
        # TODO - consistent use of f-strings
        print(f"There are {len(self.gene_to_sequence)} genes loaded: ")
        number = 1
        for k, _ in self.gene_to_sequence.items():
            print(str(number) + ": " + k)
            number += 1


def parse_genes_data(genes_file: str = "") -> List[Gene]:
    genes = []
    try:
        with open(genes_file) as File:
            lines = File.readlines()
            for line in lines:
                gene_and_sequence = line.strip().split(":")

                if len(gene_and_sequence) != 2:
                    print(
                        "Line "
                        + line
                        + " in "
                        + genes_file
                        + " is incorrectly formatted, will not parse it"
                    )
                    pass
                genes.append(Gene(gene_and_sequence[0], gene_and_sequence[1]))
        return genes
    except FileNotFoundError as e:
        # TODO :: use python logging module
        print("Could not open " + genes_file + ": " + e.strerror)
        sys.exit()


def pretty_print_count_nucleotides(
    gene_name: str = "", nucleotide_counts: Dict[str, int] = None
) -> None:
    if not nucleotide_counts:
        print("Must provide nucleotide count map to print")
    else:
        pretty_printed = ""
        for k, v in nucleotide_counts.items():
            pretty_printed += f"{k}={v} "
        print(gene_name + ": " + pretty_printed)


def cli_app() -> None:

    genes = parse_genes_data(GENES_FILE)

    gene_data = GenesData()

    for gene in genes:
        gene_data.add_gene_to_sequence(gene.gene_name, gene.sequence)
        gene_data.add_sequence_to_nucleotide_counts(gene.sequence)

    while True:
        print("** Gene Sequence Explorer **")
        print("1. List all genes")
        print("2. View sequence of gene")
        print("3. Count nucleotides of gene")
        print("4. Exit application")

        menu_choice = input("Enter choice (1-4): ").strip()

        if menu_choice == "1":
            gene_data.print_all_genes()
        elif menu_choice == "2":
            gene_name = input("Enter gene name: ").strip()
            result = gene_data.get_gene_sequence(gene_name)
            if result:
                print(result)
        elif menu_choice == "3":
            gene_name = input("Enter gene name: ").strip()
            sequence = gene_data.get_gene_sequence(gene_name)
            if sequence:
                result = gene_data.get_count_nucleotides(sequence)
                if result:
                    pretty_print_count_nucleotides(gene_name, result)
        elif menu_choice == "4":
            break
        else:
            print("Invalid choice")


if __name__ == "__main__":
    cli_app()
