import sys, logging
from pathlib import Path
from .config import path_config, logging_config
from .util import cli_util
from .gene_explorer import GenesExplorer
from .fasta_parser import FastaParsingError, parse_fasta_file
from .fastq_parser import parse_fastq_file
from .vcf_parser import show_low_confidence_variants
from .blast_client import make_blast_call, BlastDatabase, BlastProgram
from .sam_bam_parser import open_alignment_file, get_read_alignment_stats_summary

logging_config.setup_logging()

DEFAULT_GENES_FILE = "sample_genes.fasta.gz"

logger = logging.getLogger(__name__)


def main() -> None:
    logger.info("Starting app")

    try:
        genes_file = sys.argv[1] if len(sys.argv) > 1 else DEFAULT_GENES_FILE
        if cli_util.validate_file_input(genes_file, [".fasta.gz", "fastq.gz"]):
            if genes_file.endswith(".fasta.gz"):
                genes = parse_fasta_file(f"{path_config.FASTA_PATH}/{genes_file}")
            else:
                genes = parse_fastq_file(f"{path_config.FASTQ_PATH}/{genes_file}")
        else:
            print(
                f"Unaccepted file format in {genes_file}. Accepted file formats: .fasta.gz, .fastq.gz"
            )
            return
    except FastaParsingError as e:
        logger.critical(f"Problem parsing {genes_file}: {e}, exiting application")
        sys.exit(1)

    genes_explorer = GenesExplorer(genes)

    while True:
        print("Bio Data Explorer")
        print("1. List all genes")
        print("2. View sequence info of gene")
        print("3. Count nucleotides of gene")
        print("4. Run BLAST search")
        print("5. Read alignment analysis (SAM/BAM)")
        print("6. Filter variants (VCF)")
        print("7. Exit application")

        menu_choice = input("Enter choice (1-7): ").strip()

        if menu_choice == "1":
            genes_explorer.print_all_genes()
        elif menu_choice == "2":
            gene_name = input("Enter gene name: ").strip().lower()
            sequence = genes_explorer.get_gene_sequence(gene_name)
            if sequence:
                cli_util.print_sequence_info(genes_explorer, gene_name, sequence)
        elif menu_choice == "3":
            gene_name = input("Enter gene name: ").strip().lower()
            sequence = genes_explorer.get_gene_sequence(gene_name)
            if sequence:
                count_nucleotides = genes_explorer.get_count_nucleotides(sequence)
                if count_nucleotides:
                    genes_explorer.pretty_print_count_nucleotides(
                        gene_name, count_nucleotides
                    )
        elif menu_choice == "4":
            sequence = input("Enter sequence: ").strip().lower()
            print("Processing BLAST call, this takes some time...")
            blast_record = make_blast_call(
                BlastProgram.BLASTN, BlastDatabase.NT, sequence
            )
            cli_util.print_blast_record(blast_record)
        elif menu_choice == "5":
            bam_or_same_file_input = input(
                "Enter file name of alignment file in data/sam-bam/, ex: align.sam or align.bam: "
            ).strip()
            if cli_util.validate_file_input(bam_or_same_file_input, [".sam", ".bam"]):
                alignment_file = open_alignment_file(
                    f"{path_config.SAM_BAM_PATH}/{bam_or_same_file_input}"
                )
                alignment_stats = get_read_alignment_stats_summary(alignment_file)
                cli_util.print_alignment_summary(alignment_stats)
                alignment_file = open_alignment_file(
                    f"{path_config.SAM_BAM_PATH}/{bam_or_same_file_input}"
                )

                chrom = None
                start = None
                end = None

                if Path(bam_or_same_file_input).suffix == ".bam":
                    print("Region-specific exploration for .bam file...")
                    chrom = input("Enter chromosome (e.g. chr1):").strip()
                    start = input("Enter start (1-based):").strip()
                    end = input("Enter end:").strip()
                    if cli_util.sequence_input_range_is_valid(start, end):
                        cli_util.print_alignment_core_details(
                            alignment_file, chrom, int(start), int(end)
                        )
        elif menu_choice == "6":
            variant_file = input(
                "Enter file_name of compressed .vcf file in data/vcf, ex: file_name.vcf.gz:"
            ).strip()
            if cli_util.validate_file_input(variant_file, ["vcf.gz"]):
                low_confidence_variants = show_low_confidence_variants(
                    f"{path_config.VCF_PATH}/{variant_file}"
                )
                if len(low_confidence_variants) > 0:
                    cli_util.print_low_confidence_variants(low_confidence_variants)
                else:
                    print(f"{variant_file} had no low confidence variants")
        elif menu_choice == "7":
            logger.info("Exiting app")
            break
        else:
            print("Invalid choice")


if __name__ == "__main__":
    main()
