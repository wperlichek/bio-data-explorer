from _pytest.logging import LogCaptureFixture
from unittest.mock import patch
from bio_data_explorer import blast_client as bc
import logging


def test_make_blast_call_empty_sequence(caplog: LogCaptureFixture):
    with caplog.at_level(logging.WARNING):
        result = bc.make_blast_call(bc.BlastProgram.BLASTN, bc.BlastDatabase.NT, "")
    assert result is None
    assert "Must provide sequence to make a BLAST call" in caplog.text


def test_make_blast_uses_qblast():
    with patch("bio_data_explorer.blast_client.NCBIWWW.qblast") as mock_call:
        seq = "AGCTAG"
        mock_call.return_value = "mock_handle"
        result = bc.make_blast_call(bc.BlastProgram.BLASTN, bc.BlastDatabase.NT, seq)
        mock_call.assert_called_once_with(
            bc.BlastProgram.BLASTN, bc.BlastDatabase.NT, seq
        )
