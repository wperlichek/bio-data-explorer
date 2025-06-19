from bio_data_explorer import fastq_parser as fp
from .config import path_config as test_path_config


def test_parse_fastq_file_keeps_good_reads():
    file_name = "sample_genes_good_format.fastq.gz"
    expected_genes_parsed = 5
    genes = fp.parse_fastq_file(f"{test_path_config.FASTQ_PATH}/{file_name}")
    assert len(genes) == expected_genes_parsed


def test_parse_fastq_file_discards_reads_with_len_too_short():
    file_name = "sample_genes_short_length.fastq.gz"
    expected_genes_parsed = 1
    genes = fp.parse_fastq_file(f"{test_path_config.FASTQ_PATH}/{file_name}")
    assert len(genes) == expected_genes_parsed


def test_parse_fastq_file_discards_reads_with_high_unknown_base_count():
    file_name = "sample_genes_high_unknown_bases.fastq.gz"
    expected_genes_parsed = 1
    genes = fp.parse_fastq_file(f"{test_path_config.FASTQ_PATH}/{file_name}")
    assert len(genes) == expected_genes_parsed
