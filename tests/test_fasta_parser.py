import pytest
from bio_data_explorer import fasta_parser as fp
from .config import path_config as test_path_config


def test_parse_fasta_file_missing_file():
    with pytest.raises(fp.FastaParsingError):
        fp.parse_fasta_file("file_doesnt_exist.txt")


def test_parse_fasta_file_returns_all_records_if_formatted_correctly():
    file_name = "sample_genes_good_format.fasta.gz"
    expected_genes_parsed = 11
    genes = fp.parse_fasta_file(f"{test_path_config.FASTA_PATH}/{file_name}")
    assert len(genes) == expected_genes_parsed


def test_parse_fasta_file_exludes_records_missing_identifier():
    file_name = "sample_genes_bad_header.fasta.gz"
    original_record_count = 11
    expected_genes_parsed = (
        original_record_count - 1
    )  # contains 1 record with missing identifier
    genes = fp.parse_fasta_file(f"{test_path_config.FASTA_PATH}/{file_name}")
    assert len(genes) == expected_genes_parsed


def test_parse_fasta_file_exludes_records_with_non_dna_bases():
    file_name = "sample_genes_bad_sequences.fasta.gz"
    original_record_count = 11
    expected_genes_parsed = (
        original_record_count
        - 2  # contains 2 records with bases that are non-DNA bases
    )
    genes = fp.parse_fasta_file(f"{test_path_config.FASTA_PATH}/{file_name}")
    assert len(genes) == expected_genes_parsed
