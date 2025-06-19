from bio_data_explorer import fastq_parser as fqp

# def test_test_parse_fastq_file_keeps_good_reads():


# def test_quality_trim_sequence_end_retains_original_for_good_reads():
#     file_name = "sample_genes_good_format.fasta.gz"
#     expected_genes_parsed = 11
#     genes = fp.parse_fasta_file(f"{test_path_config.DATA_DIR}/{file_name}")
#     assert actual_sequence == fqp._quality_trim_sequence_end(
#         actual_sequence, [20, 20, 20, 40, 20]
#     )


# def test_quality_trim_sequence_end_trims_trailing():
#     actual_sequence = "TCTGA"
#     expected_sequence = "TC"
#     assert expected_sequence == fqp._quality_trim_sequence_end(
#         actual_sequence, [20, 20, 10, 40, 20]
#     )


# def test_quality_trim_sequence_end_trims_from_first_bad_read():
#     actual_sequence = "TCTGA"
#     expected_sequence = "TC"
#     assert expected_sequence == fqp._quality_trim_sequence_end(
#         actual_sequence, [20, 20, 10, 40, 20]
#     )
