import sys

genes_data = "genes.txt"

genes_sequence = {}  # gene_name:dna_sequence
sequences_counts = {}  # sequence: nucleotide counts


def store_nucleotide_sequence_counts(sequence=""):
    if sequence == "":
        print("No sequence to examine")
    else:
        counts = {"A": 0, "C": 0, "T": 0, "G": 0}  # nucleotide letter: count
        sequences_counts[sequence] = counts
        for idx in range(len(sequence)):
            sequences_counts[sequence][sequence[idx]] += 1


# parse text file containing gene data and store the data in maps
def parse_genes_data():
    # TODO :: better parsing, e.g. make sure all nucleotides are G, T, C, A
    try:
        with open(genes_data) as File:
            print("Opened " + genes_data)
            lines = File.readlines()
            for line in lines:
                gene_to_sequence = line.strip().split(":")
                print("Gene name: " + gene_to_sequence[0])
                print("Sequence: " + gene_to_sequence[1])

                if gene_to_sequence[0] not in genes_sequence:
                    genes_sequence[genes_data[0]] = gene_to_sequence[1]
                else:
                    print("Gene " + gene_to_sequence[0] + ", already loaded!")

                if gene_to_sequence[1] not in sequences_counts:
                    store_nucleotide_sequence_counts(gene_to_sequence[1])
                else:
                    print(
                        "There is already a sequence count for sequence "
                        + gene_to_sequence[1]
                    )
    except FileNotFoundError as e:
        print("Could not open " + genes_data + ": " + e.strerror)
        sys.exit()


parse_genes_data()
