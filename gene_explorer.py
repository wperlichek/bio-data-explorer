import sys

genes_data = "genes.txt"

genes_sequence = {} # gene_name:dna_sequence
sequences_counts = {} # sequence: nucleotide counts

def store_nucleotide_sequence_counts(sequence=""):
    if sequence == "":
        print("No sequence to examine")
    else:
        counts = {"A": 0, "C": 0, "T": 0, "G": 0} # nucleotide letter: count
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
                gene_data = line.strip().split(":")
                print("Gene name: " + gene_data[0])
                print("Sequence: " + gene_data[1])

                if gene_data[0] not in genes_sequence:
                    genes_sequence[genes_data[0]] = genes_data[1]
                else:
                    print("Gene " + genes_data[0] + ", already loaded!")
                    
                if genes_data[1] not in sequences_counts:
                    store_nucleotide_sequence_counts(genes_data[1])
                else:
                    print("There is already a sequence count for sequence " + genes_data[1])
    except FileNotFoundError as e:
        print("Could not open " + genes_data + ": " + e.strerror)
        sys.exit()

parse_genes_data()