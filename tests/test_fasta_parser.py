import pytest
from .config import path_config as test_path_config
from bio_data_explorer import fasta_parser as fp


def test_parse_fasta_file_missing_file():
    with pytest.raises(fp.FastaParsingError):
        fp.parse_fasta_file("file_doesnt_exist.txt")


def test_parse_fasta_file_exludes_records_missing_identifier():
    file_name = "sample_genes_bad_header.fasta"
    fp.parse_fasta_file(f"{test_path_config.DATA_DIR}/{file_name}")
