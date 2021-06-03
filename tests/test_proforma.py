
from ast import parse
from os import path
import unittest
import pyteomics
pyteomics.__path__ = [path.abspath(
    path.join(path.dirname(__file__), path.pardir, 'pyteomics'))]
from pyteomics import proforma
from pyteomics.proforma import (
    ProForma, TaggedInterval, parse, MassModification,
    ModificationRule, StableIsotope, GenericModification, to_proforma,
    obo_cache)


class ProFormaTest(unittest.TestCase):
    maxDiff = None

    def test_complicated_short(self):
        complicated_short = r"<[Carbamidomethyl]@C><13C>[Hydroxylation]?{HexNAc}[Hex]-ST[UNIMOD:Oxidation](EPP)[+18.15]ING"
        tokens, properties = parse(complicated_short)
        assert len(tokens) == 8
        assert len(properties['n_term']) == 1
        assert properties['n_term'][0] == 'Hex'
        assert len(properties['intervals']) == 1
        assert properties['intervals'][0] == TaggedInterval(2, 5, [MassModification(18.15)])
        assert len(properties['isotopes']) == 1
        assert properties['isotopes'][0] == StableIsotope("13C")
        assert properties['fixed_modifications'][0] == ModificationRule(
            GenericModification('Carbamidomethyl', None, None), ['C'])
        assert to_proforma(tokens, **properties) == complicated_short
        self.assertAlmostEqual(
            ProForma(tokens, properties).mass, 1210.5088, 3)


    def test_ranges(self):
        seq = "PRQT(EQC[Carbamidomethyl]FQRMS)[+19.0523]ISK"
        parsed = proforma.ProForma.parse(seq)
        assert str(parsed) == seq

    def test_error_on_nested_range(self):
        self.assertRaises(proforma.ProFormaError, lambda: parse(
            "PRQT(EQ(CFQR)[Carbamidomethyl]MS)[+19.0523]ISK"))

    def test_localization_scores(self):
        seq = "EM[Oxidation]EVT[#g1(0.01)]S[#g1(0.09)]ES[Phospho#g1(0.90)]PEK"
        obj = ProForma.parse(seq)
        tags = obj.find_tags_by_id("#g1")
        solutions = {4: 0.01, 5: 0.09, 7: 0.9}
        for i, tag in tags:
            marker = tag.find_tag_type(proforma.TagTypeEnum.localization_marker)[0]
            expected = solutions[i]
            assert expected == marker.value

    def test_multiple_info(self):
        i = proforma.ProForma.parse(
            "ELVIS[Phospho|INFO:newly discovered|info:really awesome]K")
        tags = i[4][1][0].find_tag_type(proforma.TagTypeEnum.info)
        messages = set(['newly discovered', 'really awesome'])
        assert len(tags) == 2
        for tag in tags:
            messages.remove(tag.value)
        assert len(messages) == 0

    def test_formula(self):
        i = proforma.ProForma.parse("SEQUEN[Formula:[13C2]CH6N]CE")
        mod = i[-3][1][0]
        assert mod.composition == proforma.Composition(
            {'H': 6, 'C[13]': 2, 'C': 1, 'N': 1})

    def test_gnome(self):
        gp = proforma.ProForma.parse("NEEYN[GNO:G59626AS]K")
        self.assertAlmostEqual(gp.mass, 2709.016, 3)

    def test_glycan(self):
        gp = proforma.ProForma.parse("NEEYN[Glycan:Hex5HexNAc4NeuAc1]K")
        self.assertAlmostEqual(gp.mass, 2709.016, 3)


if __name__ == '__main__':
    unittest.main()
