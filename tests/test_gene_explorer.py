import pytest
from bio_data_explorer import gene_explorer as ge


@pytest.fixture(scope="module")
def get_genes_explorer() -> ge.GenesExplorer:
    return ge.GenesExplorer()


def test_get_reverse_compliment_valid_sequence(get_genes_explorer: ge.GenesExplorer):
    positive_strand = "GCTACCCATA"
    expected_reverse_compliment = "TATGGGTAGC"
    result = get_genes_explorer.get_reverse_compliment(positive_strand)
    assert result == expected_reverse_compliment
