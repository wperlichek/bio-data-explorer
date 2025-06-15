import sys, logging
from .gene_explorer import GenesExplorer
from .fasta_parser import parse_genes_data, GenesFileParsingError
from .blast_client import make_blast_call, BlastDatabase, BlastProgram

logging.basicConfig(
    level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s"
)

DEFAULT_GENES_FILE = "sample_genes.fasta.gz"


def main() -> None:
    logging.info("Starting app")

    try:
        genes_file = sys.argv[1] if len(sys.argv) > 1 else DEFAULT_GENES_FILE
        genes = parse_genes_data(genes_file)
    except GenesFileParsingError as e:
        logging.critical(f"Application can't start due to {e}, exiting application")
        sys.exit(1)

    genes_explorer = GenesExplorer(genes)

    while True:
        print("** Gene Sequence Explorer **")
        print("1. List all genes")
        print("2. View sequence of gene")
        print("3. Count nucleotides of gene")
        print("4. Make BLAST call")
        print("5. Exit application")

        menu_choice = input("Enter choice (1-4): ").strip()

        if menu_choice == "1":
            genes_explorer.print_all_genes()
        elif menu_choice == "2":
            gene_name = input("Enter gene name: ").strip().lower()
            sequence = genes_explorer.get_gene_sequence(gene_name)
            if sequence:
                print(
                    f"{genes_explorer.get_gene_name_original_casing(gene_name)}: {sequence}"
                )
                print(f"Sequence length: {genes_explorer.get_sequence_length(sequence)}")
                print(f"Reverse compliment: {genes_explorer.get_reverse_compliment(sequence)}")
                print(f"DNA to RNA transcription: {genes_explorer.get_dna_to_rna_transcription(sequence)}")
        elif menu_choice == "3":
            gene_name = input("Enter gene name: ").strip().lower()
            sequence = genes_explorer.get_gene_sequence(gene_name)
            if sequence:
                count_nucleotides = genes_explorer.get_count_nucleotides(sequence)
                if count_nucleotides:
                    genes_explorer.pretty_print_count_nucleotides(gene_name, count_nucleotides)
        elif menu_choice == "4":
            sequence = input("Input sequence: ").strip().lower()
            make_blast_call(BlastProgram.BLASTN, BlastDatabase.NT, "ATGGCAGATTAGTGCAATGGGAGCCTTCGGAGCCATGGCCAACCTCCTCCTAGC")
        elif menu_choice == "5":
            logging.info("Exiting app")
            break
        else:
            print("Invalid choice")


if __name__ == "__main__":
    main()
