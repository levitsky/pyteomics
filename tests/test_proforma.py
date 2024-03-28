from os import path
import unittest
import pickle
import pyteomics
pyteomics.__path__ = [path.abspath(
    path.join(path.dirname(__file__), path.pardir, 'pyteomics'))]
from pyteomics.proforma import (
    PSIModModification, ProForma, TaggedInterval, parse, MassModification, ProFormaError, TagTypeEnum,
    ModificationRule, StableIsotope, GenericModification, Composition, to_proforma, ModificationMassNotFoundError,
    obo_cache, process_tag_tokens)


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

    def test_range(self):
        seq = "PRQT(EQC[Carbamidomethyl]FQRMS)[+19.0523]ISK"
        parsed = ProForma.parse(seq)
        assert str(parsed) == seq
        chunk = parsed[:6]
        assert chunk.intervals

    def test_ambiguous_range(self):
        seq = "PRQT(?EQC[Carbamidomethyl]FQRMS)ISK"
        parsed = ProForma.parse(seq)
        assert str(parsed) == seq
        self.assertRaises(ValueError, lambda: parsed[:6])

    def test_error_on_nested_range(self):
        self.assertRaises(ProFormaError, lambda: parse(
            "PRQT(EQ(CFQR)[Carbamidomethyl]MS)[+19.0523]ISK"))

    def test_localization_scores(self):
        seq = "EM[Oxidation]EVT[#g1(0.01)]S[#g1(0.09)]ES[Phospho#g1(0.90)]PEK"
        obj = ProForma.parse(seq)
        tags = obj.find_tags_by_id("#g1")
        solutions = {4: 0.01, 5: 0.09, 7: 0.9}
        for i, tag in tags:
            marker = tag.find_tag_type(TagTypeEnum.localization_marker)[0]
            expected = solutions[i]
            assert expected == marker.value

    def test_multiple_info(self):
        i = ProForma.parse(
            "ELVIS[Phospho|INFO:newly discovered|info:really awesome]K")
        tags = i[4][1][0].find_tag_type(TagTypeEnum.info)
        messages = set(['newly discovered', 'really awesome'])
        assert len(tags) == 2
        for tag in tags:
            messages.remove(tag.value)
        assert len(messages) == 0

    def test_formula(self):
        i = ProForma.parse("SEQUEN[Formula:[13C2]CH6N]CE")
        mod = i[-3][1][0]
        assert mod.composition == Composition(
            {'H': 6, 'C[13]': 2, 'C': 1, 'N': 1})

    def test_gnome(self):
        gp = ProForma.parse("NEEYN[GNO:G59626AS]K")
        self.assertAlmostEqual(gp.mass, 2709.016, 3)

    def test_glycan(self):
        gp = ProForma.parse("NEEYN[Glycan:Hex5HexNAc4NeuAc1]K")
        self.assertAlmostEqual(gp.mass, 2709.016, 3)

    def test_c_terminal_modification(self):
        i = ProForma.parse(
            "[iTRAQ4plex]-EM[U:Oxidation]EVNES[Phospho]PEK[iTRAQ4plex]-[Methyl]")
        self.assertEqual(i.c_term[0].name, "Methyl")
        self.assertEqual(i[-1][1][0].name, "iTRAQ4plex")

    def test_fragments(self):
        i = ProForma.parse("PEPTIDE")
        masses = i.fragments('b', 1)

        expected = [98.06004032, 227.1026334, 324.15539725, 425.20307572,
                    538.2871397, 653.31408272]

        for o, e in zip(masses, expected):
            self.assertAlmostEqual(o, e, 3)

        masses = i.fragments('y', 1)
        expected = [148.06043424, 263.08737726, 376.17144124, 477.21911971,
                    574.27188356, 703.31447664]

        for o, e in zip(masses, expected):
            self.assertAlmostEqual(o, e, 3)

    def test_slice(self):
        i = ProForma.parse('[U:1]-MPEP-[UNIMOD:2]/2')
        assert i.n_term is not None
        assert i.c_term is not None

        assert i[:1].n_term is not None
        assert i[:1].c_term is None

        assert i[1:].n_term is None
        assert i[1:].c_term is not None


class TestTagProcessing(unittest.TestCase):
    def test_process_tag_tokens(self):
        tokens = list('UNIMOD:Deamidation')
        tag = process_tag_tokens(tokens)
        assert tag.value == "Deamidation"
        assert tag.type == TagTypeEnum.unimod

    def test_process_tag_tokens_generic(self):
        tokens = list('Deamidation')
        tag = process_tag_tokens(tokens)
        assert tag.value == "Deamidation"
        assert tag.type == TagTypeEnum.generic

    def test_process_tag_tokens_contains_colon(self):
        tokens = list('UNIMOD:Cation:Na')
        tag = process_tag_tokens(tokens)
        assert tag.value == "Cation:Na"
        assert tag.type == TagTypeEnum.unimod

    def test_process_tag_tokens_generic_contains_colon(self):
        for name in ['Cation:Na', 'Cation:Li', 'Unknown:210', 'QAT:2H(3)',
                     'Dimethyl:2H(4)', 'Label:13C(9)', 'Cation:K']:
            tag = process_tag_tokens(list(name))
            assert tag.value == name
            assert tag.type == TagTypeEnum.generic
            state = tag.resolve()
            assert state['name'] == name
            assert state['provider'] == 'unimod'


class GenericModificationResolverTest(unittest.TestCase):
    def test_generic_resolver(self):
        mod = "Oxidation"
        state = GenericModification(mod)
        state.resolve()
        self.assertEqual(state.provider, 'unimod')
        self.assertAlmostEqual(state.mass, 15.994915, 3)

    def test_generic_resolver_with_dangerous_synonyms(self):
        mod = "TMT6plex"
        state = GenericModification(mod)
        state.resolve()
        self.assertEqual(state.provider, 'unimod')
        self.assertAlmostEqual(state.mass, 229.162932, 3)


class PSIModModificationResolverTest(unittest.TestCase):
    def test_unknown_mass(self):
        mod = "MOD:01716"  # 'TMT6plex reporter fragment'
        state = PSIModModification(mod)
        self.assertRaises(ModificationMassNotFoundError, lambda: state.resolve())


class ModificationHashingTest(unittest.TestCase):
    def test_mass_modification(self):
        mod = MassModification(57.08)

        container = set()
        container.add(mod.key)
        self.assertIn(mod.key, container)

        mod2 = MassModification(57.08 + 1e-19)
        self.assertIn(mod2.key, container)
        self.assertIn(mod2, container)


class ModificationPicklingTest(unittest.TestCase):
    def test_pickle(self):
        mod = GenericModification("UNIMOD:1")
        payload = pickle.dumps(mod)
        dup = pickle.loads(payload)
        self.assertEqual(mod, dup)
        assert mod.mass is not None
        payload = pickle.dumps(mod)
        dup = pickle.loads(payload)
        self.assertEqual(mod, dup)


if __name__ == '__main__':
    unittest.main()
