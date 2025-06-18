import sys, logging
from typing import Any, Dict, Optional
from pathlib import Path
from .gene_explorer import GenesExplorer
from .fasta_parser import FastaParsingError, parse_fasta_file
from .fastq_parser import parse_fastq_file
from .vcf_parser import show_low_confidence_variants
from .blast_client import make_blast_call, BlastDatabase, BlastProgram
from .sam_bam_parser import open_alignment_file, get_read_alignment_stats_summary
from pysam import AlignmentFile

logging.basicConfig(
    level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s"
)

DEFAULT_GENES_FILE = "sample_genes.fasta.gz"


def main() -> None:
    logging.info("Starting app")

    try:
        genes_file = sys.argv[1] if len(sys.argv) > 1 else DEFAULT_GENES_FILE
        if "fasta" in genes_file:  # TODO :: stricter parsing requirements
            genes = parse_fasta_file(genes_file)
        else:
            genes = parse_fastq_file(genes_file)
    except FastaParsingError as e:
        logging.critical(f"Application can't start due to {e}, exiting application")
        sys.exit(1)

    genes_explorer = GenesExplorer(genes)

    while True:
        print("** Bio Data Explorer **")
        print("1. List all genes")
        print("2. View sequence of gene")
        print("3. Count nucleotides of gene")
        print("4. Make BLAST call")
        print("5. Show low confidence variants in VCF file")
        print("6. Sam analysis")
        print("7. Exit application")

        menu_choice = input("Enter choice (1-7): ").strip()

        if menu_choice == "1":
            genes_explorer.print_all_genes()
        elif menu_choice == "2":
            gene_name = input("Enter gene name: ").strip().lower()
            sequence = genes_explorer.get_gene_sequence(gene_name)
            if sequence:
                print(
                    f"{genes_explorer.get_gene_name_original_casing(gene_name)}: {sequence}"
                )
                print(
                    f"Sequence length: {genes_explorer.get_sequence_length(sequence)}"
                )
                print(
                    f"Reverse compliment: {genes_explorer.get_reverse_compliment(sequence)}"
                )
                print(
                    f"DNA to RNA transcription: {genes_explorer.get_dna_to_rna_transcription(sequence)}"
                )
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
            sequence = input("Input sequence: ").strip().lower()
            logging.info("Processing BLAST call, this takes some time...")
            records = make_blast_call(BlastProgram.BLASTN, BlastDatabase.NT, sequence)
            print_blast_record(records)
        elif menu_choice == "5":
            variant_file = input(
                "Enter file_name of compressed .vcf file in data/, ex: file_name.vcf.gz:"
            )
            low_confidence_variants = show_low_confidence_variants(variant_file)
            if len(low_confidence_variants) > 0:
                print(
                    "These are the low confidence variants in format CHROM:POS_REF>ALT(s):"
                )
                for variant in low_confidence_variants:
                    print(variant)
        elif menu_choice == "6":
            bam_or_same_file_input = (
                input("Enter .sam or .bam file name: ").strip().lower()
            )
            alignment_file = open_alignment_file(bam_or_same_file_input)
            alignment_stats = get_read_alignment_stats_summary(alignment_file)
            print_alignment_summary(alignment_stats)
            alignment_file = open_alignment_file(bam_or_same_file_input)

            chrom = None
            start = None
            end = None

            if Path(bam_or_same_file_input).suffix == ".bam":
                # TODO :: this is too brittle, e.g. putting empty for int causes exception
                chrom = input("Enter chromosome (e.g. chr1):").strip().lower()
                start = int(input("Enter start (1-based):").strip())
                end = int(input("Enter end:").strip())

            print_alignment_core_details(alignment_file, chrom, start, end)
        elif menu_choice == "7":
            logging.info("Exiting app")
            break
        else:
            print("Invalid choice")


# TODO :: these "print" functions should probably be in utils
def print_blast_record(blast_records: Any = None) -> None:
    for record in blast_records:
        print("****")
        print("Record alignments: ")
        for alignment in record.alignments:
            print(f"{alignment.title}")
            print("Alignment hsps:")
            for hsp in alignment.hsps:
                # TODO :: percentage identity and coverage calculations
                # TODO :: should printing and calculating for in a util file?
                print(f"Score: {hsp.score}")
                print(f"E-Value: {hsp.expect}")
        print("****")


def print_alignment_summary(alignment_stats: Dict[str, int]):
    print("Alignment stats:")
    for k, v in alignment_stats.items():
        print(f"{k}: {v}")
    return


def print_alignment_core_details(
    alignment_file: AlignmentFile,
    chrom: Optional[str] = None,
    start: Optional[int] = None,
    end: Optional[int] = None,
):
    print("Core aligment details for each read:")
    print(f"Has index: {alignment_file.has_index()}")
    if alignment_file.has_index() and chrom and start and end:
        logging.info(
            f"Using region-specific alignment view at {chrom} in range {start}:{end}"
        )
        for read in alignment_file.fetch(chrom, start, end):
            print("***")
            print(f"QNAME: {read.query_name}")
            print(f"FLAG: {read.flag}")
            print(f"MAPQ: {read.mapping_quality}")
            print(f"CIGAR: {read.cigarstring}")
            print("***")
    else:
        for read in alignment_file.fetch(until_eof=True):
            print("***")
            print(f"QNAME: {read.query_name}")
            print(f"FLAG: {read.flag}")
            print(f"MAPQ: {read.mapping_quality}")
            print(f"CIGAR: {read.cigarstring}")
            print("***")


if __name__ == "__main__":
    main()
