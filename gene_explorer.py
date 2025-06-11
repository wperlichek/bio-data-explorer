import sys

genes_data = "genes.txt"

gene_to_sequence = {}
sequence_to_nucleotide_counts = {}


def count_nucleotides_in_sequence(sequence="") -> None:
    if sequence == "":
        print("Sequence must be a non-empty string")
    else:
        counts = {"A": 0, "C": 0, "T": 0, "G": 0}
        sequence_to_nucleotide_counts[sequence] = counts
        for idx in range(len(sequence)):
            nucleotide = sequence[idx]
            if (
                nucleotide != "A"
                and nucleotide != "C"
                and nucleotide != "T"
                and nucleotide != "G"
            ):
                print(
                    "Warning: unknown nucleotide present in sequence "
                    + sequence
                    + " at position "
                    + str(idx)
                    + " :"
                    + nucleotide
                )
                pass
            else:
                sequence_to_nucleotide_counts[sequence][sequence[idx]] += 1


def parse_genes_data() -> None:
    try:
        with open(genes_data) as File:
            lines = File.readlines()
            for line in lines:
                gene_and_sequence = line.strip().split(":")
                gene_name = gene_and_sequence[0]
                gene_sequence = gene_and_sequence[1]
                if gene_name not in gene_to_sequence:
                    gene_to_sequence[gene_name] = gene_sequence
                else:
                    print("Gene " + gene_name + ", already loaded")

                if gene_to_sequence[1] not in sequence_to_nucleotide_counts:
                    count_nucleotides_in_sequence(gene_sequence)
                else:
                    print(
                        "There is already a sequence count for sequence "
                        + gene_sequence
                    )
    except FileNotFoundError as e:
        print("Could not open " + genes_data + ": " + e.strerror)
        sys.exit()


def print_all_genes() -> None:
    print("There are " + str(len(gene_to_sequence)) + " genes loaded: ")
    number = 1
    for k, v in gene_to_sequence.items():
        print(str(number) + ": " + k)
        number += 1


def get_gene_sequence(gene_name="") -> str:
    if gene_name == "":
        print("Must provide gene name to view sequence")
        return None
    elif gene_name not in gene_to_sequence:
        print("Gene not found")
        return None
    else:
        return gene_name + " sequence: " + gene_to_sequence[gene_name]


def get_count_nucleotides(gene_name="") -> str:
    if gene_name == "":
        print("Must provide gene name to count nucleotides")
        return None
    elif gene_name not in gene_to_sequence:
        print("Gene not found")
        return None
    else:
        sequence = gene_to_sequence[gene_name]
        return (
            gene_name
            + " nucleotide count(s): "
            + str(sequence_to_nucleotide_counts[sequence])
        )  # TODO :: clean up the output so it's more human readable


def cli_app() -> None:

    parse_genes_data()

    while True:
        print("Gene Sequence Explorer")
        print("1. List all genes")
        print("2. View gene sequence of gene")
        print("3. Count nucleotides of gene")
        print("4. Exit application")

        user_input = input("Enter choice (1-4): ").strip()
        res = ""

        if user_input == "1":
            print_all_genes()
        elif user_input == "2":
            user_input = input("Enter gene name: ")
            res = get_gene_sequence(user_input).strip()
            print(res)
        elif user_input == "3":
            user_input = input("Enter gene name: ").strip()
            res = get_count_nucleotides(user_input)
            print(res)
        elif user_input == "4":
            sys.exit()
        else:
            print("Invalid choice")


cli_app()
