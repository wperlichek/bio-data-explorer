import sys, logging
from gene_explorer import GenesExplorer
from fasta_parser import parse_genes_data, GenesFileParsingError

logging.basicConfig(
    level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s"
)

GENES_FILE = "sample_genes.fasta.gz"

def main() -> None: 
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
    main()