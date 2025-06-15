import pytest
from bio_data_explorer import fasta_parser as fp


def test_line_is_formatted_correctly_missing_identifier():
    with pytest.raises(fp.GenesFileParsingError) as exc_info:
        fp.parse_genes_data("file_doesnt_exist.txt")
