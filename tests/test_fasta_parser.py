import pytest
from bio_data_explorer import fasta_parser as fp


def test_parse_fasta_file_missing_file():
    with pytest.raises(fp.FastaParsingError):
        fp.parse_fasta_file("file_doesnt_exist.txt")


def test_parse_fasta_file_exludes_records_missing_identifier():
    return
