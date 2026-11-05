import numpy as np
import pandas as pd
import pytest

from nanoremix.tables import sequence81


def row(short="AAAACUUUU"):
    return pd.DataFrame([{"center_ref_base": "C", "canonical_context": short}])


def test_explicit_sequence_and_full_context_check():
    frame = row()
    frame["seq"] = "A" * 36 + "AAAACTTTT" + "T" * 36
    seq, status = sequence81(frame, sequence_column="seq")
    assert status[0] == "ok" and len(seq[0]) == 81
    frame["seq"] = "A" * 40 + "C" + "A" * 40
    _, status = sequence81(frame, sequence_column="seq")
    assert status[0] == "81nt_context_mismatch"


@pytest.mark.parametrize("coordinate_base", [0, 1])
def test_fasta_explicit_coordinate_system(coordinate_base):
    frame = row()
    frame["id"] = "ref"
    frame["center"] = 40 + coordinate_base
    seq, status = sequence81(frame, fasta={"ref": "A" * 40 + "C" + "T" * 40},
                             id_column="id", center_column="center", coordinate_base=coordinate_base)
    assert status[0] == "ok"


def test_negative_strand_is_reverse_complemented():
    frame = row()
    frame["id"], frame["center"], frame["strand"] = "ref", 40, "-"
    forward = "A" * 40 + "C" + "T" * 40
    reverse = forward.translate(str.maketrans("ACGT", "TGCA"))[::-1]
    sequences, status = sequence81(frame, fasta={"ref": reverse}, id_column="id",
                                   center_column="center", strand_column="strand", coordinate_base=0)
    assert status[0] == "ok"
    assert sequences[0] == forward


def test_edges_are_not_padded():
    frame = row()
    frame["id"], frame["center"] = "ref", 4
    _, status = sequence81(frame, fasta={"ref": "A" * 100}, id_column="id",
                           center_column="center", coordinate_base=0)
    assert status[0] == "insufficient_81nt_flanks"


@pytest.mark.parametrize("coordinate", ["not-a-number", None, float("inf")])
def test_bad_coordinate_has_a_machine_readable_error_status(coordinate):
    frame = pd.DataFrame([{
        "ref": "transcript", "center": coordinate, "center_ref_base": "A", "canonical_context": "",
    }])
    sequences, statuses = sequence81(frame, fasta={"transcript": "A" * 100},
                                    id_column="ref", center_column="center", coordinate_base=0)
    assert statuses[0].startswith("invalid")
    assert sequences[0] == ""


def test_missing_sequence_metadata_is_explicit():
    frame = pd.DataFrame([{"center_ref_base": "A", "canonical_context": ""}])
    _, statuses = sequence81(frame, fasta={"transcript": "A" * 100},
                            id_column="ref", center_column="center", coordinate_base=0)
    assert statuses[0] == "missing_sequence_metadata"
