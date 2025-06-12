from sys import exit
from typing import Dict, Optional

GENES_FILE = "genes.txt"

# TODO::refactor into class
# TODO:: should these be caps in the class?
gene_to_sequence: Dict[str, str] = {}
sequence_to_nucleotide_counts: Dict[str, Dict[str,int]] = {}


def count_nucleotides_in_sequence(sequence: str = "") -> None:
    if sequence == "":
        print("Must provide non-empty sequence to count nucleotides")
    else:
        counts = {"A": 0, "C": 0, "T": 0, "G": 0}
        sequence_to_nucleotide_counts[sequence] = counts
        for nucleotide in sequence:
            if nucleotide in counts:
                sequence_to_nucleotide_counts[sequence][nucleotide] += 1
            else:
                print(
                    "Warning: Unknown nucleotide present in sequence "
                    + sequence
                    + ": "
                    + nucleotide
                )
                pass
                


def parse_genes_data() -> None:
    try:
        with open(GENES_FILE) as File:
            lines = File.readlines()
            for line in lines:
                gene_and_sequence = line.strip().split(":")

                if (len(gene_and_sequence) != 2):
                    print("Line " + line + " in " + GENES_FILE + 
                          " is incorrectly formatted, will not parse it")
                    pass

                gene_name = gene_and_sequence[0]
                gene_sequence = gene_and_sequence[1]
                if gene_name not in gene_to_sequence:
                    gene_to_sequence[gene_name] = gene_sequence
                else:
                    print("Gene " + gene_name + " already loaded")

                if gene_sequence not in sequence_to_nucleotide_counts:
                    count_nucleotides_in_sequence(gene_sequence)
                else:
                    print(
                        "Warning: There is already a nucleotide count for sequence "
                        + gene_sequence
                    )
    except FileNotFoundError as e:
        # TODO :: use python logging module
        print("Could not open " + GENES_FILE + ": " + e.strerror)
        sys.exit()


def print_all_genes() -> None:
    # TODO - consistent use of f-strings
    print(f"There are " + {len(gene_to_sequence)} + " genes loaded: ")
    number = 1
    for k, v in gene_to_sequence.items():
        print(str(number) + ": " + k)
        number += 1


def get_gene_sequence(gene_name: str = "") -> Optional[str]:
    if gene_name == "":
        print("Must provide gene name to view sequence")
        return None
    elif gene_name not in gene_to_sequence:
        print("Gene not found")
        return None
    else:
        return gene_name + " sequence: " + gene_to_sequence[gene_name]


def get_count_nucleotides(gene_name: str = "") -> Optional[Dict[str, int]]:
    if gene_name == "":
        print("Must provide gene name to count nucleotides")
        return None
    elif gene_name not in gene_to_sequence:
        print("Gene not found")
        return None
    else:
        sequence = gene_to_sequence[gene_name]
        if sequence not in sequence_to_nucleotide_counts:
            print("Sequence not found")
        else:
            return sequence_to_nucleotide_counts[sequence]


def pretty_print_count_nucleotides(gene_name: str = "", nucleotide_counts: Optional[Dict[str, int]] = None) -> None:
    if not nucleotide_counts:
        print("Must provide nucleotide count map to print")
    else:
        pretty_printed = ""
        for k, v in nucleotide_counts.items():
            pretty_printed += f"{k}={v} "
        print(gene_name + ": " + pretty_printed)


def cli_app() -> None:

    parse_genes_data()

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

# TODO :: wrap in standard python startup idiom
cli_app()
