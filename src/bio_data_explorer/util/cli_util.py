import logging
from typing import Dict, Iterator, List, Optional
from Bio.Blast import Record
from pysam import AlignmentFile
from ..gene_explorer import GenesExplorer

logger = logging.getLogger(__name__)


def print_blast_record(blast_record: Optional[Iterator[Record]] = None) -> None:
    if not blast_record:
        logger.warning("Blast record is empty can't print it")
        return
    for hit in blast_record:
        print("****")
        print("Record alignments: ")
        for alignment in hit.alignments:  # type: ignore
            print(f"{alignment.title}")
            print("Alignment hsps:")
            for hsp in alignment.hsps:
                print(f"Score: {hsp.score}")
                print(f"E-Value: {hsp.expect}")
        print("****")


def print_sequence_info(
    genes_explorer: GenesExplorer, gene_name: str = "", sequence: str = ""
):
    if genes_explorer is None or not gene_name or not sequence:
        logger.warning(
            "Must provide valid genes_explorer, gene_name, and sequence to print sequence info"
        )
        return
    else:
        print(f"{genes_explorer.get_gene_name_original_casing(gene_name)}: {sequence}")
        print(f"Sequence length: {genes_explorer.get_sequence_length(sequence)}")
        print(f"Reverse compliment: {genes_explorer.get_reverse_compliment(sequence)}")
        print(
            f"DNA to RNA transcription: {genes_explorer.get_dna_to_rna_transcription(sequence)}"
        )


def print_low_confidence_variants(low_confidence_variants: List[str]):
    print("These are the low confidence variants in format CHROM:POS_REF>ALT(s):")
    for variant in low_confidence_variants:
        print(variant)


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
        logger.info(
            f"Using region-specific alignment view at {chrom} in range {start}:{end}"
        )
        for read in alignment_file.fetch(chrom, start, end):
            print("***")
            print(f"QNAME: {read.query_name}")
            print(f"SEQ: {read.query_sequence}")
            print(f"FLAG: {read.flag}")
            print(f"MAPQ: {read.mapping_quality}")
            print(f"CIGAR: {read.cigarstring}")
            print("***")
    else:
        for read in alignment_file.fetch(until_eof=True):
            print("***")
            print(f"QNAME: {read.query_name}")
            print(f"SEQ: {read.query_sequence}")
            print(f"FLAG: {read.flag}")
            print(f"MAPQ: {read.mapping_quality}")
            print(f"CIGAR: {read.cigarstring}")
            print("***")
