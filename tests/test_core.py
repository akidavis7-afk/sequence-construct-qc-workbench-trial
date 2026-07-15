from src.core import motif_hits

def test_iupac_motif():
    hits=motif_hits('AATATATATAAA',[{'name':'t','pattern':'TATAWA'}])
    assert 't' in hits


def test_equal_length_fast_variant_call(tmp_path):
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    from src.core import align_and_call

    reference = SeqRecord(Seq("AACCGGTT"), id="ref")
    sample = SeqRecord(Seq("AACCTGTT"), id="sample")
    variants = align_and_call(reference, sample)
    assert len(variants) == 1
    assert variants[0].kind == "SNV"
    assert variants[0].ref_start == 5
