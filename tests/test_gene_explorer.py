import pytest
from bio_data_explorer import gene_explorer as ge


@pytest.fixture(scope="module")
def genes_explorer() -> ge.GenesExplorer:
    return ge.GenesExplorer()


def test_get_reverse_compliment_valid_sequence(genes_explorer):
    positive_strand = "GCTACCCATA"
    expected_reverse_compliment = "TATGGGTAGC"
    result = genes_explorer.get_reverse_compliment(positive_strand)
    assert result == expected_reverse_compliment


def test_add_gene_sequence_adds_gene(genes_explorer):
    genes_explorer.add_gene_sequence("GeneA", "desc", "ACTG")
    seq = genes_explorer.get_gene_sequence("GeneA")
    assert seq == "ACTG"


def test_add_sequence_nucleotide_counts_adds_counts(genes_explorer):
    seq = "ACTG"
    genes_explorer.add_sequence_nucleotide_counts(seq)
    counts = genes_explorer.get_count_nucleotides(seq)
    assert counts == {"A": 1, "C": 1, "T": 1, "G": 1}


def test_get_gene_sequence_returns_none_for_missing(genes_explorer):
    result = genes_explorer.get_gene_sequence("NoSuchGene")
    assert result is None


def test_get_gene_name_original_casing_returns_original(genes_explorer):
    genes_explorer.add_gene_sequence("GeneB", sequence="TTTT")
    original = genes_explorer.get_gene_name_original_casing("geneb")
    assert original == "GeneB"


def test_get_sequence_length_returns_correct_length(genes_explorer):
    length = genes_explorer.get_sequence_length("ACTG")
    assert length == 4


def test_get_reverse_compliment_returns_correct_result(genes_explorer):
    result = genes_explorer.get_reverse_compliment("ATCG")
    assert result == "CGAT"


def test_get_dna_to_rna_transcription_transcribes_t_to_u(genes_explorer):
    rna = genes_explorer.get_dna_to_rna_transcription("ATTGC")
    assert rna == "AUUGC"


def test_count_nucleotides_in_sequence_counts_correctly(genes_explorer):
    counts = genes_explorer.count_nucleotides_in_sequence("AACTG")
    assert counts == {"A": 2, "C": 1, "T": 1, "G": 1}


def test_get_count_nucleotides_returns_none_for_unknown_sequence(genes_explorer):
    result = genes_explorer.get_count_nucleotides("NNNN")
    assert result is None


def test_get_gc_percentage_computes_percentage_correctly(genes_explorer):
    counts = {"A": 1, "C": 1, "T": 1, "G": 1}
    gc_pct = genes_explorer.get_gc_percentage(counts)
    assert gc_pct == 50.0
