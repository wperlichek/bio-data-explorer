from _pytest.logging import LogCaptureFixture
from bio_data_explorer import blast_client as bc
import logging


def test_make_blast_call_empty_sequence(caplog: LogCaptureFixture):
    with caplog.at_level(logging.WARNING):
        result = bc.make_blast_call(bc.BlastProgram.BLASTN, bc.BlastDatabase.NT, "")
    assert result is None
    assert "Must provide sequence to make a BLAST call" in caplog.text
