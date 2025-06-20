from pysam import AlignmentFile
from unittest.mock import patch, MagicMock
from bio_data_explorer import sam_bam_parser as sbp


def test_get_read_alignment_summary_returns_alignment_summary():
    mock_read1 = MagicMock()
    mock_read2 = MagicMock()
    mock_bam_file = [mock_read1, mock_read2]

    with patch("your_module.AlignmentFile", return_value=mock_bam_file):
        summary = sbp.get_read_alignment_stats_summary(mock_bam_file)
