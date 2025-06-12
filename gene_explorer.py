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

    def get_gene_sequence(self) -> Optional[str]:
        if self.gene_name == "":
            print("Must provide gene name to view sequence")
            return None
        elif self.gene_name not in self.gene_to_sequence:
            print("Gene not found")
            return None
        else:
            return (
                self.gene_name + " sequence: " + self.gene_to_sequence[self.gene_name]
            )

    def get_gene_sequence(self) -> Optional[str]:
        if self.gene_name == "":
            print("Must provide gene name to view sequence")
            return None
        elif self.gene_name not in self.gene_to_sequence:
            print("Gene not found")
            return None
        else:
            return (
                self.gene_name + " sequence: " + self.gene_to_sequence[self.gene_name]
            )

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

    def get_count_nucleotides(self) -> Optional[Dict[str, int]]:
        if self.gene_name == "":
            print("Must provide gene name to count nucleotides")
            return None
        elif self.gene_name not in self.gene_to_sequence:
            print("Gene not found")
            return None
        else:
            sequence = self.ene_to_sequence[self.gene_name]
            if sequence not in self.sequence_to_nucleotide_counts:
                print("Sequence not found")
            else:
                return self.sequence_to_nucleotide_counts[sequence]

    def print_all_genes(self) -> None:
        # TODO - consistent use of f-strings
        print(f"There are {len(self.gene_to_sequence)} genes loaded: ")
        number = 1
        for k, v in self.gene_to_sequence.items():
            print(str(number) + ": " + k)
            number += 1

    def pretty_print_count_nucleotides(self) -> None:
        if not self.nucleotide_counts:
            print("Must provide nucleotide count map to print")
        else:
            pretty_printed = ""
            for k, v in self.nucleotide_counts.items():
                pretty_printed += f"{k}={v} "
            print(self.gene_name + ": " + pretty_printed)


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


def cli_app() -> None:

    genes = parse_genes_data()

    gene_data = GenesData()

    for gene in genes:
        gene_data.add_gene_to_sequence(gene.gene_name, gene.sequence)
        gene_data.add_sequence_to_nucleotide_counts(
            gene.sequence,
        )

    while True:
        print("Gene Sequence Explorer")
        print("1. List all genes")
        print("2. View gene sequence of gene")
        print("3. Count nucleotides of gene")
        print("4. Exit application")

        user_input = input("Enter choice (1-4): ").strip()
        result = ""

        if user_input == "1":
            print_all_genes()
        elif user_input == "2":
            user_input = input("Enter gene name: ").strip()
            result = get_gene_sequence(user_input)
            if result:
                print(result)
        elif user_input == "3":
            user_input = input("Enter gene name: ").strip()
            result = get_count_nucleotides(user_input)
            if result:
                pretty_print_count_nucleotides(user_input, result)
        elif user_input == "4":
            break
        else:
            print("Invalid choice")


if __name__ == "__main__":
    cli_app()
