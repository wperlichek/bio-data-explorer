from unittest.mock import MagicMock
from bio_data_explorer import sam_bam_parser as sbp


def test_get_read_alignment_summary_returns_correct_alignment_summary():
    mock_alignment_file = MagicMock()

    read1 = MagicMock()
    read1.is_mapped = True
    read1.mapping_quality = 60
    read1.is_duplicate = False
    read1.is_supplementary = False
    read1.query_name = "read1"

    read2 = MagicMock()
    read2.is_mapped = True
    read2.mapping_quality = 10
    read2.is_duplicate = True
    read2.is_supplementary = True
    read2.query_name = "read2"

    read3 = MagicMock()
    read3.is_mapped = False
    read3.mapping_quality = 60
    read3.is_duplicate = False
    read3.is_supplementary = False
    read3.query_name = "read3"

    mock_alignment_file = MagicMock()
    mock_alignment_file.__iter__.return_value = [read1, read2, read3]

    result = sbp.get_read_alignment_stats_summary(mock_alignment_file)

    assert result == {
        "reads_count": 3,
        "mapped_reads": 2,
        "unmapped_reads": 1,
        "mapping_rate_percentage": 67.0,
        "low_quality_reads": 1,
        "duplicate_reads": 1,
        "supplementary_reads": 1,
    }
