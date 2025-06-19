from _pytest.logging import LogCaptureFixture
from unittest.mock import patch
from bio_data_explorer import blast_client as bc
import logging
from typing import Iterator
from Bio.Blast import Record


def test_make_blast_call_empty_sequence(caplog: LogCaptureFixture):
    with caplog.at_level(logging.WARNING):
        result = bc.make_blast_call(bc.BlastProgram.BLASTN, bc.BlastDatabase.NT, "")
    assert result is None
    assert "Must provide sequence to make a BLAST call" in caplog.text


def test_make_blast_gets_blast_search_results():
    with (
        patch("bio_data_explorer.blast_client.NCBIWWW.qblast") as mock_qblast_call,
        patch("bio_data_explorer.blast_client.NCBIXML.parse") as mock_parse_call,
    ):
        seq = "AGCTAG"
        mock_handle = "mock_handle"
        stub_record_iterator: Iterator[Record] = iter([Record()])
        mock_qblast_call.return_value = mock_handle
        mock_parse_call.return_value = stub_record_iterator
        result = bc.make_blast_call(bc.BlastProgram.BLASTN, bc.BlastDatabase.NT, seq)
        mock_qblast_call.assert_called_once_with(
            bc.BlastProgram.BLASTN.value, bc.BlastDatabase.NT.value, seq
        )
        mock_parse_call.assert_called_once_with(mock_handle)
        assert result is stub_record_iterator
