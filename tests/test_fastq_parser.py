from bio_data_explorer import fastq_parser as fp


def test_quality_trim_sequence_end_retains_original_for_good_reads():
    actual_sequence = "TCTGA"
    assert actual_sequence == fp.quality_trim_sequence_end(
        actual_sequence, [20, 20, 20, 40, 20]
    )


def test_quality_trim_sequence_end_trims_trailing():
    actual_sequence = "TCTGA"
    expected_sequence = "TC"
    assert expected_sequence == fp.quality_trim_sequence_end(
        actual_sequence, [20, 20, 10, 40, 20]
    )


def test_quality_trim_sequence_end_trims_from_first_bad_read():
    actual_sequence = "TCTGA"
    expected_sequence = "TC"
    assert expected_sequence == fp.quality_trim_sequence_end(
        actual_sequence, [20, 20, 10, 40, 20]
    )
