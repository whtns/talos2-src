import pathlib

"""
Tests for the two behavioural gates this deployment adds to run_small_filtering.

Both default to off, so the first assertion of each group is that stock behaviour is
unchanged. See patches_talos2/0006 and 0010 in the talos deployment repository for the
measurements that motivated them.
"""

from talos2.run_small_filtering import (
    AUTOSOME_OR_PAR,
    HET,
    HOM_ALT,
    HOM_REF,
    MITO,
    Y_NONPAR,
    DeNovoConfig,
    candidate_configuration,
    category_high_impact,
    parent_has_evidence,
)

CRITICAL = {'frameshift', 'stop_gained'}


def csq(consequence='frameshift', biotype='protein_coding', mane_id=''):
    return {'consequence': consequence, 'biotype': biotype, 'mane_id': mane_id}


# --- patch 0006: a parent's hom-ref needs a real call behind it --------------------


def test_parent_evidence_off_accepts_anything():
    """Default config: an uncalled parent still counts, which is upstream behaviour."""
    conf = DeNovoConfig()
    assert conf.require_parent_evidence is False
    assert parent_has_evidence(gq=-1, has_depth=False, conf=conf) is True


def test_parent_evidence_rejects_merge_filler():
    """`bcftools merge -0` filler: GT 0/0, no GQ (cyvcf2 reports -1), no depth."""
    conf = DeNovoConfig(require_parent_evidence=True)
    assert parent_has_evidence(gq=-1, has_depth=False, conf=conf) is False
    # a real GQ with no depth is still not two pieces of evidence
    assert parent_has_evidence(gq=50, has_depth=False, conf=conf) is False
    assert parent_has_evidence(gq=50, has_depth=True, conf=conf) is True


def test_parent_evidence_applies_the_gq_floor():
    conf = DeNovoConfig(require_parent_evidence=True, min_parent_gq=40)
    assert parent_has_evidence(gq=39, has_depth=True, conf=conf) is False
    assert parent_has_evidence(gq=40, has_depth=True, conf=conf) is True


def test_candidate_configuration_requires_both_parents_on_an_autosome():
    args = (AUTOSOME_OR_PAR, None, HET, HOM_REF, HOM_REF)
    assert candidate_configuration(*args) is True
    assert candidate_configuration(*args, dad_called=False) is False
    assert candidate_configuration(*args, mom_called=False) is False


def test_candidate_configuration_only_tests_the_parent_each_region_reads():
    """On Y the mother is not read, so her evidence cannot matter; mito is the mirror."""
    assert candidate_configuration(Y_NONPAR, False, HOM_ALT, HOM_REF, None, mom_called=False) is True
    assert candidate_configuration(Y_NONPAR, False, HOM_ALT, HOM_REF, None, dad_called=False) is False
    assert candidate_configuration(MITO, None, HOM_ALT, None, HOM_REF, dad_called=False) is True
    assert candidate_configuration(MITO, None, HOM_ALT, None, HOM_REF, mom_called=False) is False


# --- patch 0010: High Impact requires a MANE transcript ---------------------------


def test_high_impact_off_fires_on_any_transcript():
    assert category_high_impact([csq()], CRITICAL) == 1


def test_high_impact_with_mane_required():
    assert category_high_impact([csq()], CRITICAL, require_mane=True) == 0
    assert category_high_impact([csq(mane_id='NM_001234.5')], CRITICAL, require_mane=True) == 1


def test_high_impact_exempts_snrna():
    """MANE is protein-coding only, so a strict gate would discard every snRNA report."""
    assert category_high_impact([csq(biotype='snRNA')], CRITICAL, require_mane=True) == 1


def test_high_impact_needs_the_mane_transcript_to_carry_the_consequence():
    """A MANE transcript with a harmless consequence does not license a minor one's."""
    consequences = [
        csq(consequence='frameshift', mane_id=''),
        csq(consequence='synonymous', mane_id='NM_001234.5'),
    ]
    assert category_high_impact(consequences, CRITICAL) == 1
    assert category_high_impact(consequences, CRITICAL, require_mane=True) == 0


# --- patch 0007: a carrier needs a genotype quality worth reporting ----------------
#
# These drive the real create_small_variant against upstream's own labelled fixture,
# whose single carrier is called at GQ 99, rather than re-implementing the predicate
# in the test. config_retrieve is patched because the floor is read from TALOS_CONFIG,
# which conftest points at upstream's test config.


def small_variant_with_floor(monkeypatch, floor):
    """create_small_variant over the labelled fixture, with min_carrier_gq = floor."""
    from cyvcf2 import VCF

    from talos2 import utils

    real = utils.config_retrieve

    def fake(key, default=None, **kwargs):
        if key == ['ValidateMOI', 'min_carrier_gq']:
            return floor
        return real(key, default, **kwargs) if default is not None else real(key)

    monkeypatch.setattr(utils, 'config_retrieve', fake)

    reader = VCF(str(pathlib.Path(__file__).parent / 'input' / '1_labelled_variant.vcf.bgz'))
    variant = utils.create_small_variant(next(reader), reader.samples)
    reader.close()
    return variant


def test_carrier_gq_floor_off_keeps_the_carrier():
    """floor 0 is the default, and the fixture's GQ 99 carrier is reported."""
    import pytest

    mp = pytest.MonkeyPatch()
    try:
        variant = small_variant_with_floor(mp, 0)
    finally:
        mp.undo()
    assert variant.het_samples == {'male'}


def test_carrier_gq_floor_drops_a_carrier_below_it():
    """The same GQ 99 carrier disappears once the floor is above its GQ."""
    import pytest

    mp = pytest.MonkeyPatch()
    try:
        variant = small_variant_with_floor(mp, 100)
    finally:
        mp.undo()
    assert variant.het_samples == set()
    assert variant.hom_samples == set()
    # and it cannot come back through a sample category
    for category in variant.sample_categories:
        assert 'male' not in variant.info[category]


def test_carrier_gq_floor_keeps_the_carrier_just_at_the_floor():
    """The comparison is `gq < floor`, so a carrier exactly at the floor survives."""
    import pytest

    mp = pytest.MonkeyPatch()
    try:
        variant = small_variant_with_floor(mp, 99)
    finally:
        mp.undo()
    assert variant.het_samples == {'male'}
