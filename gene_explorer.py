import sys

genes_data = "genes.txt"

genes_sequence = {}  # gene_name:dna_sequence
sequences_counts = {}  # sequence: nucleotide counts


def store_nucleotide_sequence_counts(sequence="") -> None:
    if sequence == "":
        print("No sequence to examine")
    else:
        counts = {"A": 0, "C": 0, "T": 0, "G": 0}  # nucleotide letter: count
        sequences_counts[sequence] = counts
        for idx in range(len(sequence)):
            sequences_counts[sequence][sequence[idx]] += 1


# parse text file containing gene data and store the data in maps
def parse_genes_data() -> None:
    # TODO :: better parsing, e.g. make sure all nucleotides are G, T, C, A
    try:
        with open(genes_data) as File:
            lines = File.readlines()
            for line in lines:
                gene_to_sequence = line.strip().split(":")
                gene_name = gene_to_sequence[0]
                gene_sequence = gene_to_sequence[1]
                if gene_name not in genes_sequence:
                    genes_sequence[gene_name] = gene_sequence
                else:
                    print("Gene " + gene_name + ", already loaded!")

                if gene_to_sequence[1] not in sequences_counts:
                    store_nucleotide_sequence_counts(gene_sequence)
                else:
                    print(
                        "There is already a sequence count for sequence "
                        + gene_sequence
                    )
    except FileNotFoundError as e:
        print("Could not open " + genes_data + ": " + e.strerror)
        sys.exit()


def list_all_genes() -> None:
    print("There are " + str(len(genes_sequence)) + " genes loaded: ")
    number = 1
    for k, v in genes_sequence.items():
        print(str(number) + ": " + k)
        number += 1


def view_gene_sequence(gene_name="") -> str:
    if gene_name == "":
        print("Must provide gene name to view sequence")
        return None
    elif gene_name not in genes_sequence:
        print("Gene not found")
        return None
    else:
        return gene_name + " sequence: " + genes_sequence[gene_name]


def count_nucleotides(gene_name="") -> str:
    if gene_name == "":
        print("Must provide gene name to count nucleotides")
        return None
    elif gene_name not in genes_sequence:
        print("Gene not found")
        return None
    else:
        sequence = genes_sequence[gene_name]
        return (
            gene_name + " nucleotide count(s): " + str(sequences_counts[sequence])
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
            list_all_genes()
        elif user_input == "2":
            user_input = input("Enter gene name: ")
            res = view_gene_sequence(user_input).strip()
            print(res)
        elif user_input == "3":
            user_input = input("Enter gene name: ").strip()
            res = count_nucleotides(user_input)
            print(res)
        elif user_input == "4":
            sys.exit()
        else:
            print("Invalid choice.")


cli_app()
