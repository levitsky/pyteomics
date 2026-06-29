import os
from os import path
import unittest
import pickle
import math
import pyteomics
pyteomics.__path__ = [path.abspath(
    path.join(path.dirname(__file__), path.pardir, 'pyteomics'))]
from pyteomics.proforma import (
    PSIModModification, ProForma, TaggedInterval, parse, MassModification, ProFormaError, TagTypeEnum,
    ModificationRule, StableIsotope, GenericModification, Composition, to_proforma, ModificationMassNotFoundError,
    UnimodModification, ModificationTarget,
    AdductParser, ChargeState, proteoforms, _coerce_string_to_modification,
    std_aa_comp, obo_cache, process_tag_tokens, peptidoforms)
from pyteomics import mass


CACHE_PATH = os.environ.get('OBO_CACHE_PATH')
if CACHE_PATH:
    obo_cache.cache_path = CACHE_PATH

obo_cache.enabled = bool(os.environ.get("OBO_CACHE_ENABLED"))


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
        # The newer serializer enforces the ordering of sections
        assert (
            to_proforma(tokens, **properties)
            == r"<13C><[Carbamidomethyl]@C>[Hydroxylation]?{HexNAc}[Hex]-ST[UNIMOD:Oxidation](EPP)[+18.15]ING"
        )
        self.assertAlmostEqual(ProForma(tokens, properties).mass, 1228.6588, 3)

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

        b_expected = [98.06004032, 227.1026334, 324.15539725, 425.20307572,
                    538.2871397, 653.31408272]

        for o, e in zip(masses, b_expected):
            self.assertAlmostEqual(o, e, 3)

        masses = i.fragments('y', 1)
        y_expected = [148.06043424, 263.08737726, 376.17144124, 477.21911971,
                    574.27188356, 703.31447664]

        for o, e in zip(masses, y_expected):
            self.assertAlmostEqual(o, e, 3)

        # Test include labile
        i = ProForma.parse("{+204}PEPTIDE")
        for o, e in zip(i.fragments('y', include_labile=True), masses):
            self.assertAlmostEqual(o, e + 204, 3)
        for o, e in zip(i.fragments("y", include_labile=False), masses):
            self.assertAlmostEqual(o, e, 3)

        # Test include unlocalized
        i = ProForma.parse("[+204]?PEPTIDE")
        for o, e in zip(i.fragments("y", include_unlocalized=True), masses):
            self.assertAlmostEqual(o, e + 204, 3)
        for o, e in zip(i.fragments("y", include_unlocalized=False), masses):
            self.assertAlmostEqual(o, e, 3)

        # Test C-terminal modification
        i = ProForma.parse("PEPTIDE-[+204]")
        for o, e in zip(i.fragments("y"), y_expected):
            self.assertAlmostEqual(o, e + 204, 3)
        for o, e in zip(i.fragments("b"), b_expected):
            self.assertAlmostEqual(o, e, 3)

        # Test N-terminal modification
        i = ProForma.parse("[+204]-PEPTIDE")
        for o, e in zip(i.fragments("b"), b_expected):
            self.assertAlmostEqual(o, e + 204, 3)
        for o, e in zip(i.fragments("y"), y_expected):
            self.assertAlmostEqual(o, e, 3)

        # Test fixed modifications
        i = ProForma.parse("<[+204]@P>PEPTIDE")
        for o, e in zip(i.fragments("b"), b_expected):
            d = 204
            if (o - e) > 205:
                d += 204
            self.assertAlmostEqual(o, e + d, 3)

        # Test regular modifications
        i = ProForma.parse("P[+204]EP[+204]TIDE")
        for o, e in zip(i.fragments("b"), b_expected):
            d = 204
            if (o - e) > 205:
                d += 204
            self.assertAlmostEqual(o, e + d, 3)

        i = ProForma.parse("(PEP)[+204](TIDE)[+204]")
        for o, e in zip(i.fragments("b"), b_expected):
            d = 204
            if (o - e) > 205:
                d += 204
            self.assertAlmostEqual(o, e + d, 3)
        self.assertEqual(d, 408)

    def test_slice(self):
        seq = ProForma.parse('[U:1]-MPEP-[UNIMOD:2]/2')
        assert seq.n_term is not None
        assert seq.c_term is not None

        assert seq[:1].n_term is not None
        assert seq[:1].c_term is None

        assert seq[1:].n_term is None
        assert seq[1:].c_term is not None

        seq = ProForma.parse("MPE[#1]PET[+204#1]ID[#1]E")
        sub = seq[:2]
        assert not sub.find_tags_by_id('1')

        sub = seq[:3]
        of = sub.find_tags_by_id('1')
        assert of
        assert of[0][0] == 2
        self.assertAlmostEqual(of[0][1].mass, 204.0, 3)



    def test_charge_adducts(self):
        sequences = ['PEPTIDE/1[+2Na+,-H+]', 'PEPTIDE/-1[+e-]', 'PEPTIDE/1[+2H+,+e-]']
        charges = [1, -1, 1]
        adducts_list = [[('Na', 1, 2), ('H', 1, -1)], [('e-', -1, 1)], [('H', 1, 2), ('e-', -1, 1)]]
        for seq, charge, adducts in zip(sequences, charges, adducts_list):
            i = ProForma.parse(seq)
            self.assertEqual(i.charge_state.charge, charge)
            self.assertEqual(i.charge_state.adducts, adducts)

    def test_chimeric_parse_requires_opt_in(self):
        with self.assertRaisesRegex(ProFormaError, 'chimeric=True'):
            parse('PEPTIDE+ELVIS')

    def test_chimeric_parse(self):
        seq = 'EMEVEESPEK/2+ELVISLIVER/3'
        parsed = parse(seq, chimeric=True)
        self.assertEqual(len(parsed), 2)
        self.assertEqual(parsed[0][1]['charge_state'].charge, 2)
        self.assertEqual(parsed[1][1]['charge_state'].charge, 3)

        forms = ProForma.parse(seq, chimeric=True)
        self.assertEqual(len(forms), 2)
        self.assertTrue(all(isinstance(form, ProForma) for form in forms))
        self.assertEqual(forms[0].charge_state.charge, 2)
        self.assertEqual(forms[1].charge_state.charge, 3)

    def test_chimeric_single_component_opt_in(self):
        forms = ProForma.parse('PEPTIDE/+2', chimeric=True)
        self.assertEqual(len(forms), 1)
        self.assertEqual(forms[0].charge_state.charge, 2)

    def test_chimeric_empty_component(self):
        with self.assertRaisesRegex(ProFormaError, 'Empty peptidoform'):
            parse('+PEPTIDE', chimeric=True)
        with self.assertRaisesRegex(ProFormaError, 'Empty peptidoform'):
            parse('PEPTIDE+', chimeric=True)

    def test_chimeric_shared_fixed_modifications(self):
        seq = '<[Carbamidomethyl]@C>AC+CC'
        parsed = parse(seq, chimeric=True)
        self.assertEqual(len(parsed), 2)
        self.assertEqual(len(parsed[0][1]['fixed_modifications']), 1)
        self.assertEqual(len(parsed[1][1]['fixed_modifications']), 1)

        forms = ProForma.parse(seq, chimeric=True)
        aa_comp = std_aa_comp.copy()
        aa_comp['cam'] = Composition(formula='H3C2NO')
        self.assertEqual(forms[0].composition(), Composition(sequence='AcamC', aa_comp=aa_comp))
        self.assertEqual(forms[1].composition(), Composition(sequence='camCcamC', aa_comp=aa_comp))

    def test_chimeric_shared_names_and_isotopes(self):
        parsed = parse('(>>>pair)(>sample)<13C>AC+CC', chimeric=True)
        self.assertEqual(parsed[0][1]['names'], {1: 'sample', 3: 'pair'})
        self.assertEqual(parsed[1][1]['names'], {3: 'pair'})
        self.assertEqual(parsed[0][1]['isotopes'], [StableIsotope('13C')])
        self.assertEqual(parsed[1][1]['isotopes'], [StableIsotope('13C')])
        parsed = ProForma.parse("(>>>pair)(>sample)<13C>AC+CC", chimeric=True)
        self.assertEqual(str(parsed[0]), "(>>>pair)<13C>(>sample)AC")
        self.assertEqual(str(parsed[1]), "(>>>pair)<13C>CC")

    def test_empty_name(self):
        for i in range(1, 4):
            p = ProForma.parse("({})PEPTIDE".format(">" * i))
            assert p.names[i] == ''
            assert str(p) == "({})PEPTIDE".format(">" * i)

    def test_chimeric_adduct_separator(self):
        forms = ProForma.parse('PEPTIDE/[Na:z+1,H:z+1]+ELVIS/2', chimeric=True)
        self.assertEqual(len(forms), 2)
        self.assertEqual(forms[0].charge_state.charge, 2)
        self.assertEqual(forms[1].charge_state.charge, 2)

    def test_composition_with_adducts(self):
        sequences = ['PEPTIDE/1[+2Na+,-H+]', 'PEPTIDE/-1[+e-]', 'PEPTIDE/1[+2H+,+e-]', 'PEPTIDE', 'PEPTIDE/1']
        neutral_comp = Composition(sequence='PEPTIDE')
        adducts_list = [Composition({'Na': 2, 'H': -1, 'e-': -1}),
                        Composition({'e-': 1}),
                        Composition({'H': 2, 'e-': -1}),
                        Composition({}),
                        Composition({'H': 1, 'e-': -1})]
        for seq, adducts in zip(sequences, adducts_list):
            with self.subTest(f'{seq} + {adducts}'):
                i = ProForma.parse(seq)
                self.assertEqual(i.composition(), neutral_comp)
                self.assertEqual(i.composition(include_charge=True), neutral_comp + adducts)

    def test_adduct_formatting(self):
        ap = AdductParser()
        ap.buffer.extend('+2Na+')
        ap.bound()
        ap.buffer.extend('-H+')
        c = ChargeState.from_adducts(ap())
        self.assertEqual(str(c), "[Na:z+1^2,H:z+1]")
        ap = AdductParser("Na:z+1^2")
        ap.bound()
        ap.extend("H:z+1")
        c = ChargeState.from_adducts(ap())
        self.assertEqual(str(c), "[Na:z+1^2,H:z+1]")

    def test_default_adduct_formatting(self):
        c = ChargeState(2, None)
        self.assertEqual(str(c), '2')

    def test_composition_fixed(self):
        sequences = ['<[UNIMOD:4]@C>ATPEILTCNSIGCLK']
        aa_comp = std_aa_comp.copy()
        aa_comp['cam'] = Composition(formula='H3C2NO')
        comps = [Composition(sequence='ATPEILTcamCNSIGcamCLK', aa_comp=aa_comp)]

        for seq, comp in zip(sequences, comps):
            i = ProForma.parse(seq)
            self.assertEqual(i.composition(), comp)

    def test_missing_composition(self):
        sequences = ['P[+79.966]EPTIDE']
        comps = [Composition(sequence='PEPTIDE')]
        for seq, comp in zip(sequences, comps):
            i = ProForma.parse(seq)
            self.assertEqual(i.composition(ignore_missing=True), comp)
            with self.assertRaises(ProFormaError):
                ProForma.parse(seq).composition()

    def test_localization_composition(self):
        seq0 = "EMEVT[Phospho]SESPEK"
        test_seq = [
            "[Phospho]?EMEVTSESPEK",
            "EMEVT[#g1]S[#g1]ES[Phospho#g1]PEK",
            "EMEV(TS)[Phospho]ESPEK"
        ]
        base_comp = ProForma.parse(seq0).composition()
        for seq in test_seq:
            with self.subTest(seq=seq):
                i = ProForma.parse(seq)
                self.assertEqual(i.composition(), base_comp)

    def test_proteoform_copiable(self):
        test_seq = [
            "[Phospho]?EMEVTSESPEK",
            "EMEVT[#g1]S[#g1]ES[Phospho#g1]PEK",
            "EMEV(TS)[Phospho]ESPEK"
        ]
        for seq in test_seq:
            with self.subTest(seq=seq):
                i = ProForma.parse(seq)
                icopy = i.copy()
                self.assertEqual(i, icopy)
                self.assertEqual(i.mass, icopy.mass)
                self.assertEqual(i.composition(), icopy.composition())

    def test_from_spec(self):
        positive = [
            "AA",
            "A[+1]",
            "AA[+1]",
            "A(AAAA)[+1][+1]",
            "UWAKJDNLASNOIJPojkjjdakjn[U:Oxidation]",
            "[+1]-A[+1]-[+1]",
            # "AA+AA",
            "EMK[XLMOD:02000#XL1]EVTKSE[XLMOD:02010#XL2]SK[#XL1]PEK[#XL2]AR",
            # "SEK[XLMOD:02001#XL1]UENCE//EMEVTK[XLMOD:02001#XL1]SESPEK",
            "EM[Oxidation]EVEES[Phospho]PEK",
            # "EM[R: Methionine sulfone]EVEES[O-phospho-L-serine]PEK",
            "EMEVTK[X:DSS#XL1]SESPEK",
            "EM[U:Oxidation]EVEES[U:Phospho]PEK",
            "EM[+15.9949]EVEES[+79.9663]PEK",
            "EM[U:+15.995]EVEES[U:+79.966]PEK",
            "EM[U:+15.995]EVEES[Obs:+79.978]PEK",
            "RTAAX[+367.0537]WT",
            "{Glycan:Hex}EM[Oxidation]EVNES[Phospho]PEK[iTRAQ4plex]",
            "[iTRAQ4plex]-EM[Oxidation]EVNES[Phospho]PEK",
            "[iTRAQ4plex]-EM[Oxidation]EVNES[Phospho]PEK[iTRAQ4plex]-[Methyl]",
            "<[S-carboxamidomethyl-L-cysteine]@C>ATPEILTCNSIGCLK",
            "<[MOD:01090]@C>ATPEILTCNSIGCLK",
            "[Phospho]?EM[Oxidation]EVTSESPEK",
            "[Phospho][Phospho]?[Acetyl]-EM[Oxidation]EVTSESPEK",
            "EM[Oxidation]EVT[#g1]S[#g1]ES[Phospho#g1]PEK",
            "EM[Oxidation]EVT[#g1(0.01)]S[#g1(0.09)]ES[Phospho#g1(0.90)]PEK",
            "[Phospho#s1]?EM[Oxidation]EVT[#s1(0.01)]S[#s1(0.90)]ES[#s1(0.90)]PEK",
            "PROT(EOSFORMS)[+19.0523]ISK",
            "PROT(EOC[Carbamidomethyl]FORMS)[+19.0523]ISK",
            "SEQUEN[Formula:C12H20O2]CE",
            "SEQUEN[Formula:HN-1O2]CE",
            "SEQUEN[Formula:[13C2][12C-2]H2N]CE",
            "SEQUEN[Glycan:HexNAc]CE",
            "EMEVTK[XLMOD:02001#XL1]SESPEK[#XL1]",
            "EMEVTK[XLMOD:02001#XL1]SESPEK",
            # "SEK[XLMOD:02001#XL1]UENCE//EMEVTK[XLMOD:02001#XL1]SESPEK",
            # "ETFGD[MOD:00093#BRANCH]//R[#BRANCH]ATER",
            "(?DQ)NGTWEM[Oxidation]ESNENFEGYM[Oxidation]K",
            "ELVIS[Phospho|+79.966331]K",
            "ELVIS[Phospho|Obs:+79.978]K",
            "ELV[INFO:xxxxx]IS",
            "ELVIS[Phospho|INFO:newly discovered|INFO:really awesome]K",
            "ELVIS[Phospho|INFO:newly discovered|INFO:Created on 2021-06]K",
            "ELVIS[Phospho|INFO:newly discovered|INFO:Created by software Tool1]K",
            "<13C>ATPEILTVNSIGQLK",
            "EMEVEESPEK/2",
            "EMEVEESPEK+ELVISLIVER",
            "EMEVEESPEK/2+ELVISLIVER/3",
            # "A[X:DSS#XL1]//B[#XL1]+C[X:DSS#XL1]//D[#XL1]",
            "<[Carbamidomethyl]@C>ATPEILTCNSIGCLK",
            "<[Oxidation]@C,M>MTPEILTCNSIGCLK",
            "<[TMT6plex]@K,N-term>ATPEILTCNSIGCLK",
            "<[TMT6plex]@K,N-term:A>ATPEILTCNSIGCLK",
            "<[TMT6plex]@K,N-term:A,N-term:B>ATPEILTCNSIGCLK",
            "EM[Oxidation]EVEES[Phospho]PEK",
            "EM[L-methionine sulfoxide]EVEES[O-phospho-L-serine]PEK",
            # "EM[R: L-methionine sulfone]EVEES[O-phospho-L-serine]PEK", # don't support RESID
            "EMEVTK[X:DSS#XL1]SESPEK",
            "NEEYN[GNO:G59626AS]K",
            "NEEYN[G:G59626AS]K",
            "EM[U:Oxidation]EVEES[U:Phospho]PEK",
            "EM[M:L-methionine sulfoxide]EVEES[M:O-phospho-L-serine]PEK",
            "EM[U:Oxidation]EVEES[M:O-phospho-L-serine]PEK",
            "EM[Oxidation]EVEES[O-phospho-L-serine]PEK",
            "EM[Oxidation]EVE[Cation:Mg[II]]ES[Phospho]PEK",
            "EM[MOD:00719]EVEES[MOD:00046]PEK",
            "EM[UNIMOD:35]EVEES[UNIMOD:56]PEK",
            # "EM[RESID:AA0581]EVEES[RESID:AA0037]PEK", # don't support RESID
            "EMEVTK[XLMOD:02001#XL1]SESPEK[#XL1]",
            "EMK[XLMOD:02000#XL1]EVTKSE[XLMOD:02010#XL2]SK[#XL1]PEK[#XL2]AR",
            "EMEVTK[XLMOD:02001#XL1]SESPEK",
            "EMEVTK[XLMOD:02001]SESPEK",
            # "SEK[XLMOD:02001#XL1]UENCE//EMEVTK[XLMOD:02001#XL1]SESPEK",
            # "SEK[XLMOD:02001#XL1]UENCE//EMEVTK[#XL1]SESPEK",
            "EVTSEKC[MOD:00034#XL1]LEMSC[#XL1]EFD",
            "EVTSEKC[L-cystine (cross-link)#XL1]LEMSC[#XL1]EFD",
            "FVNQHLC[MOD:00034#XL1]GSHLVEALYLVC[MOD:00034#XL2]GERGFFYTPK",
            # "A//GIVEQC[MOD:00034#XL3]C[#XL1]TSIC[#XL3]SLYQLENYC[#XL2]N",
            "EVTSEKC[XLMOD:02009#XL1]LEMSC[#XL1]EFD",
            "EVTSEKC[X:Disulfide#XL1]LEMSC[#XL1]EFD",
            "EVTSEKC[half cystine]LEMSC[half cystine]EFD",
            "EVTSEKC[MOD:00798]LEMSC[MOD:00798]EFDEVTSEKC[MOD:00798]LEMSC[MOD:00798]EFD",
            "EVTSEKC[UNIMOD:374#XL1]LEMSC[#XL1]EFD",
            "EVTSEKC[Dehydro#XL1]LEMSC[#XL1]EFD",
            # "ETFGD[MOD:00093#BRANCH]//R[#BRANCH]ATER",
            # "AVTKYTSSK[MOD:00134#BRANCH]//AGKQLEDGRTLSDYNIQKESTLHLVLRLRG-[#BRANCH]",
            "NEEYN[GNO:G59626AS]K",
            "YPVLN[GNO:G62765YT]VTMPN[GNO:G02815KT]NSNGKFDK",
            "EM[+15.9949]EVEES[+79.9663]PEK",
            "EM[+15.995]EVEES[-18.01]PEK",
            "EM[U:+15.9949]EVEES[U:+79.9663]PEK",
            "EM[U:+15.995]EVEES[U:+79.966]PEK",
            "EM[U:+15.995]EVEES[Obs:+79.978]PEK",
            "EM[U:+15.995]EVEES[Obs:+79.978]PEK",
            "RTAAX[+367.0537]WT",
            "SEQUEN[Formula:C12H20O2]CE",
            "SEQUEN[Formula:[13C2]CH6N]CE",
            "SEQUEN[Formula:[13C2][12C-2]H2N]CE",
            "SEQUEN[Glycan:HexNAc1Hex2]CE",
            "[iTRAQ4plex]-EM[Oxidation]EVNES[Phospho]PEK",
            "[iTRAQ4plex]-EM[U:Oxidation]EVNES[Phospho]PEK[iTRAQ4plex]-[Methyl]",
            "{Glycan:Hex}EM[U:Oxidation]EVNES[Phospho]PEK[iTRAQ4plex]",
            "{Glycan:Hex}[iTRAQ4plex]-EM[Oxidation]EVNES[Phospho]PEK[iTRAQ4plex]",
            "{Glycan:Hex}[iTRAQ4plex]-EM[Oxidation]EVNES[Phospho]PEK[iTRAQ4plex]-[Methyl]",
            "{Glycan:Hex}{Glycan:NeuAc}EMEVNESPEK",
            "[Phospho]?EM[Oxidation]EVTSESPEK",
            "[Phospho][Phospho]?[Acetyl]-EM[Oxidation]EVTSESPEK",
            "[Phospho]^2?[Acetyl]-EM[Oxidation]EVTSESPEK",
            "[Phospho]^2?[Acetyl]-EM[Oxidation]EVTSESPEK",
            "EM[Oxidation]EVT[#g1]S[#g1]ES[Phospho#g1]PEK",
            "PRT(ESFRMS)[+19.0523]ISK",
            "PRT(EC[Carbamidomethyl]FRMS)[+19.0523]ISK",
            "EM[Oxidation]EVT[#g1(0.01)]S[#g1(0.09)]ES[Phospho#g1(0.90)]PEK",
            "[Phospho#s1]?EM[Oxidation]EVT[#s1(0.01)]S[#s1(0.09)]ES[#s1(0.90)]PEK",
            "MPGLVDSNPAPPESQEKKPLK(PCCACPETKKARDACIIEKGEEHCGHLIEAHKECMRALGFKI)[Oxidation][Oxidation][half cystine][half cystine]",
            "<13C>ATPEILTVNSIGQLK",
            "<15N>ATPEILTVNSIGQLK",
            "<D>ATPEILTVNSIGQLK",
            "<13C><15N>ATPEILTVNSIGQLK",
            "<[S-carboxamidomethyl-L-cysteine]@C>ATPEILTCNSIGCLK",
            "<[MOD:01090]@C>ATPEILTCNSIGCLK",
            "<[Oxidation]@C,M>MTPEILTCNSIGCLK",
            "<[MOD:01090]@C>[Phospho]?EM[Oxidation]EVTSECSPEK",
            "<[MOD:01090]@C>[Acetyl]-EM[Oxidation]EVTSECSPEK",
            "(?DQ)NGTWEM[Oxidation]ESNENFEGYM[Oxidation]K",
            "(?N)NGTWEM[Oxidation]ESNENFEGYM[Oxidation]K",
            "ELV[INFO:AnyString]IS",
            "ELV[info:AnyString]IS",
            "ELVIS[Phospho|INFO:newly discovered]K",
            "ELVIS[Phospho|INFO:newly discovered|INFO:really awesome]K",
            "ELVIS[Phospho|INFO:newly discovered|INFO:Created on 2021-06]K",
            "ELVIS[Phospho|INFO:newly discovered|INFO:Created by software Tool1]K",
            "ELVIS[U:Phospho|+79.966331]K",
            "ELVIS[U:Phospho|Obs:+79.978]K",
            "ELVIS[Phospho|O-phospho-L-serine]K",
            "ELVIS[UNIMOD:21|MOD:00046]K",
            "ELVIS[UNIMOD:21|Phospho]K",
            "ELVIS[Phospho|O-phospho-L-serine|Obs:+79.966]K",
            "ELVIS[Obs:+79.966|Phospho|Sulfo]K",
            "EMEVEESPEK/2",
            "EM[U:Oxidation]EVEES[U:Phospho]PEK/3",
            "[U:iTRAQ4plex]-EM[U:Oxidation]EVNES[U:Phospho]PEK[U:iTRAQ4plex]-[U:Methyl]/3",
            "EMEVEESPEK/2+ELVISLIVER/3",
            "AA(?AA)",
            "AA(?AA)AA",
            "[dehydro]^3?[gln->pyro-glu]-QSC",
            "[deamidated#1]-FEEAQ[#1]A",
            "[#1]-FEEAQ[deamidated#1]A",
            "AHAM[oxidation#1]TEG-[#1]",
            "AHAM[#1]TEG-[oxidation#1]",
            "SEQUEN[Formula:Zn1:z+2]CE",
            "<[TMT6plex]@K,N-term>ATPEILTCNSIGCLK",
            "<[Oxidation]@W,C-term:G>QATPEILTWCNSIGCLKG",
            "<[Gln->pyro-Glu]@N-term:Q><[Oxidation]@W,C-term:G>QATPEILTWCNSIGCLKG",
            "<[Amidated]@C-term>QATPEILTWCNSIGCLKG",
            "PEPTID-[a-type-ion]",
            "PEPTID[Formula:H-1C-1O-2|Info:d-ion]-[a-type-ion]",
            "PEPTIDE/[Na:z+1]",
            "PEPTIDE/[Na:z+1,H:z+1]",
            "PEPTIDE/[Na:z+1^2]",
            "PEPT[Formula:Zn:z+2]IDE/[Na:z+1^2]",
            "PE[Cation:Al[III]]PTIDE/2",
            "PE[Formula:Al H-3:z+1]PTIDE/1",
            "PE[Formula:Al H-3:z+1]PTIDE/[H:z+1]",
            "[Cation:Al[III]]?PEPTIDE/2",
            "PEPTIDE/[Al H-3:z+1,H:z+1]",
            "PEPTIDEG-[Methyl][Amidated]",
            "[Acetyl][Carbamyl]-QPEPTIDE",
            "[Formula:Zn:z+2|Position:N-term,C-term]^5[Carbamidomethyl|Position:C]^5?MDPETCPCPSGGSCTCADSCKCEGCKCTSCKKSCCSCCPAECEKCAKDCVCKGGEAAEAEAEKCSCCQ",
            "PEPTI(MERMERMERM)[Oxidation|Position:M][Oxidation|Position:M]DE",
            "PEPTI(MERMERMERM)[+32|Position:E]PEPTIDE",
            "PETIEM[Dioxidation#1][Oxidation#2]REM[#1][#2]REM[#2]RM[#1]PEPTIDE",
            "[Oxidation|CoMKP]?PEPT[Phospho]IDE",
            "(>Trypsin)AANSIPYQVSLNS",
            # "(>P07225 Vitamin K-dependent protein S OS=Homo sapiens OX=9606 GN=PROS1 PE=1 (SV=1) RANGE=12..42)GGK[xlink:dss[138]#XLDSS]IEVQLK//(>P07225 Vitamin K-dependent protein S OS=Homo sapiens OX=9606 GN=PROS1 PE=1 SV=1)KVESELIK[#XLDSS]PINPR/4",
            # "(>>>Trastuzumab Fab and coeluting Fc)(>>Fab)(>Heavy chain)EVQLVESGGGLVQPGGSLRLSC[M:l-cystine (cross-link)#XL1]AASGFNIKDTYIHWVRQAPGKGLEWVARIYPTNGYTRYADSVKGRFTISADTSKNTAYLQMNSLRAEDTAVYYC[#XL1]SRWGGDGFYAMDYWGQGTLVTVSSASTKGPSVFPLAPSSKSTSGGTAALGC[M:l-cystine (cross-link)#XL2]LVKDYFPEPVTVSWNSGALTSGVHTFPAVLQSSGLYSLSSVVTVPSSSLGTQTYIC[#XL2]NVNHKPSNTKVDKKVEPKSC[M:l-cystine (cross-link)#XL3]DKT//(>Light chain)DIQMTQSPSSLSASVGDRVTITC[M:l-cystine (cross-link)#XL4]RASQDVNTAVAWYQQKPGKAPKLLIYSASFLYSGVPSRFSGSRSGTDFTLTISSLQPEDFATYYC[#XL4]QQHYTTPPTFGQGTKVEIKRTVAAPSVFIFPPSDEQLKSGTASVVC[M:l-cystine (cross-link)#XL5]LLNNFYPREAKVQWKVDNALQSGNSQESVTEQDSKDSTYSLSSTLTLSKADYEKHKVYAC[#XL5]EVTHQGLSSPVTKSFNRGEC[#XL3]+(>Fc)HTCPPCPAPELLGGPSVFLFPPKPKDTLMISRTPEVTCVVVDVSHEDPEVKFNWYVDGVEVHNAKTKPREEQYNSTYRVVSVLTVLHQDWLNGKEYKCKVSNKALPAPIEKTISKAKGQPREPQVYTLPPSREEMTKNQVSLTCLVKGFYPSDIAVEWESNGQPENNYKTTPPVLDSDGSFFLYSKLTVDKSRWQQGNVFSCSVMHEALHNHYTQKSLSLSPGK",
        ]
        for seq in positive:
            with self.subTest(seq=seq):
                parsed = ProForma.parse(seq, chimeric='+' in seq)
                assert parsed is not None

    def test_nonstandard_amino_acid(self):
        seq = ProForma.parse("PEPTX[MOD:01001]IDE")
        bad_seq = ProForma.parse("PEPTXIDE")
        assert seq.mass != bad_seq.mass
        self.assertAlmostEqual(seq.mass, 884.4127280267099, 4)

    def test_charged_tags(self):
        seq = ProForma.parse("SEQUEN[Formula:Zn1:z+2]CE")
        assert seq.tags[0].charge == 2
        assert seq._local_charges() == (2, 1)
        assert seq.charge_state.charge == 2

        # While charge state is 2+, when serializing,
        # do not count that charge at the end of the sequence
        assert str(seq) == "SEQUEN[Formula:Zn1:z+2]CE"

        # Adapted from the spec
        template = "SEQUEN[Formula:Zn:z+2]CE/[Na:z+1^2]"
        mixed = ProForma.parse(template)
        self.assertAlmostEqual(mixed.mass, seq.mass, 6)
        # Check the physical charge state is the sum
        assert mixed.charge_state.charge == 4
        # Check that it collapses back down correctly rather
        # than double counting
        assert str(mixed) == template

        self.assertWarns(UserWarning, lambda: mixed.mz(charge=2))

        template = "<[Formula:Zn:z+2]@E>SEQUENCE"
        seq = ProForma.parse(template)
        assert seq.charge_state == 6

    def test_mass(self):
        sequences = ["PEPTIDE", "PEPTIDE/2"]
        for seq in sequences:
            with self.subTest(seq=seq):
                parsed = ProForma.parse(seq)
                self.assertAlmostEqual(parsed.mass, mass.fast_mass(sequences[0]))
        parsed = ProForma.parse("<[+204]@P>PEPTIDE")
        self.assertAlmostEqual(parsed.mass, mass.fast_mass(sequences[0]) + 408)

        seq = "<[+{}]@X>XEXTIDE".format(mass.std_aa_mass['P'])
        with self.subTest(seq):
            parsed = ProForma.parse(seq)
            self.assertAlmostEqual(parsed.mass, mass.fast_mass(sequences[0]))

    def test_terminal_mass(self):
        parsed = ProForma.parse("[+22]-PEPTIDE-[+26]")
        ref = ProForma.parse("PEPTIDE")
        self.assertAlmostEqual(parsed.mass, mass.fast_mass('PEPTIDE') + 48)
        for a, b in zip(parsed.fragments('b'), ref.fragments('b')):
            self.assertAlmostEqual(a, b + 22)
        for a, b in zip(parsed.fragments("y"), ref.fragments("y")):
            self.assertAlmostEqual(a, b + 26)

    def test_charge_settable(self):
        t = "PEPTIDE"

        parsed = ProForma.parse(t)
        self.assertAlmostEqual(parsed.mass, mass.fast_mass(t))
        self.assertEqual(str(parsed), t)

        parsed.charge_state = 1
        self.assertAlmostEqual(parsed.mass, mass.fast_mass(t))
        self.assertAlmostEqual(parsed.mz(), mass.fast_mass(t, charge=1))
        self.assertEqual(str(parsed), t + '/1')

        parsed.charge_state = 2
        self.assertAlmostEqual(parsed.mass, mass.fast_mass(t))
        self.assertAlmostEqual(parsed.mz(), mass.fast_mass(t, charge=2))
        self.assertEqual(str(parsed), t + '/2')

        parsed.charge_state = None
        self.assertAlmostEqual(parsed.mass, mass.fast_mass(t))
        self.assertEqual(str(parsed), t)

        parsed.charge_state = ChargeState(2)
        self.assertAlmostEqual(parsed.mass, mass.fast_mass(t))
        self.assertAlmostEqual(parsed.mz(), mass.fast_mass(t, charge=2))
        self.assertEqual(str(parsed), t + "/2")

    def test_position_labels(self):
        t = "PETIEM[Dioxidation#1][Oxidation#2]REM[#1][#2]REM[#2]RM[#1]PEPTIDE"
        seq = ProForma.parse(t)
        tags = seq.find_tags_by_id('2')
        self.assertEqual(len(tags), 3)
        self.assertEqual(str(seq), t)
        t = "[Dioxidation#1]?PETIE(MREMREMRM)[#1][Oxidation#2]PEPTIDE"
        seq = ProForma.parse(t)
        tags = seq.find_tags_by_id('2')
        # Matches the interval
        self.assertEqual(len(tags), 1)
        tags = seq.find_tags_by_id("1")
        # Matches the unlocalized tag and the interval
        self.assertEqual(len(tags), 2)


    def test_glycan_composition_resolution(self):
        seqs = [
            ("NEEYN[Glycan:Hex5HexNAc5NeuAc1]K", 2912.0957972884694),
            ("NEEYN[Glycan:{C6H12N4O2S1}5HexNAc4NeuAc1]K", 2919.0929805042692),
            ("NEEYN[Glycan:{+204.068}5HexNAc4NeuAc1]K", 2919.0924972884695),
            ("NEEYN[Glycan:HexHexHex3HexNAc5NeuAc1]K", 2912.0957972884694),
        ]
        for seq, mass_of in seqs:
            parsed = ProForma.parse(seq)
            self.assertAlmostEqual(parsed.mass, mass_of, 2)

        self.assertRaises(
            ValueError, lambda: ProForma.parse("NEEYN[Glycan:HexHexHex3HexNAc5NeuAc1Kxo]K").mass
        )

    def test_post_interval_tag(self):
        seqs = [
            "PEPTI(DE)[INFO:foo]",
            "PEPTI(DE)[INFO:foo]-[INFO:bar]",
            "PEPTI(DE)[INFO:foo]/2",
            "PEPTI(DE)[INFO:foo]+PEPTI(DE)[INFO:foo]",
        ]
        for i, s in enumerate(seqs):
            with self.subTest("seq={s}".format(s=s)):
                [p, *rest] = ProForma.parse(s, chimeric=True)
                if i == 1:
                    assert p.c_term
                elif i == 2:
                    assert p.charge_state == 2
                elif i == 3:
                    assert rest

    def test_mz(self):
        self.assertAlmostEqual(ProForma.parse("PEPTIDE/2").mz(), mass.fast_mass("PEPTIDE", charge=2), 5)

        seq = ProForma.parse("SEQUEN[Formula:Zn1:z+2]CE")
        adducted = ProForma.parse("SEQUENCE/[Zn1:z+2]")

        # Charged formulae contribute to m/z just like adducts do
        self.assertAlmostEqual(seq.mz(), adducted.mz(), 4)
        # But they also contribute to the total mass, unlike adducts
        self.assertAlmostEqual(seq.mass, adducted.mass + adducted.charge_state.for_mz_calculation()[0], 4)

        electron = Composition({"e-": 1}).mass()
        # Verify that the direct mass calculation is correct
        self.assertAlmostEqual(seq.mass, mass.fast_mass("SEQUENCE") + Composition("Zn").mass() - electron * 2)

        proton = Composition("").mass(charge=1)
        salted = ProForma.parse("SEQUEN[Cation:Zn[II]]CE/2")
        # Per discussion, a charged modification will be a greater mass than
        # its neutral salt cousin by as many protons as it has positive charges.
        self.assertAlmostEqual(seq.mass, salted.mass + proton * 2, 4)

        self.assertRaises(ProFormaError, lambda: ProForma.parse("PEPTIDE").mz())

    def test_parse_unresolved(self):
        p = ProForma.parse("PEPT[Cmm]IDE")
        assert p
        assert str(p) == "PEPT[Cmm]IDE"

    def test_to_proforma_with_incomplete_signature(self):
        seq = to_proforma([("I", []), ("P", [])], charge_state=ChargeState(2))
        assert seq == "IP/2"
        seq = to_proforma([("I", []), ("P", [])], charge_state=2)
        assert seq == "IP/2"


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

    def test_tag_limit(self):
        tokens = list('Phospho|Position:S|Position:T|comup|Limit:2|Limit:3')
        tag = process_tag_tokens(tokens)
        with self.assertWarns(UserWarning):
            effective_limit = tag.limit
        self.assertEqual(effective_limit, 2)


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


class ModificationTest(unittest.TestCase):
    def test_mass_modification_hashable(self):
        mod = MassModification(57.08)

        container = set()
        container.add(mod.key)
        self.assertIn(mod.key, container)

        mod2 = MassModification(57.08 + 1e-19)
        self.assertIn(mod2.key, container)
        self.assertIn(mod2, container)

    def test_generic_modifications_copiable(self):
        mod = GenericModification("Phospho")
        modcopy = mod.copy()
        self.assertEqual(mod, modcopy)

    def test_mass_modifications_copiable(self):
        mod = MassModification(57.08)
        modcopy = mod.copy()
        self.assertEqual(mod, modcopy)

    def test_resolve_unimod_by_alias(self):
        mod = UnimodModification("U:Acetylation").resolve()
        self.assertEqual(mod['name'], 'Acetyl')
        mod = GenericModification("Acetylation").resolve()
        self.assertEqual(mod["name"], "Acetyl")


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


class ProteoformCombinatorTest(unittest.TestCase):
    def test_range(self):
        seq = "EMEV(TS)[Phospho]ESPEK"
        pf = ProForma.parse(seq)
        for include_unmodified in [False, True]:
            for deepcopy in [False, True]:
                with self.subTest(include_unmodified=include_unmodified, deepcopy=deepcopy):
                    proteoforms = list(pf.proteoforms(include_unmodified=include_unmodified, deepcopy=deepcopy))
                    self.assertEqual(len(proteoforms), 2 + include_unmodified)   # Phospho on T or S (+ no phospho if include_unmodified)

    def test_unlocalized_position_list_and_count(self):
        k = 2
        seq = f"[Phospho|Position:S|Position:T]^{k}?EMEVTSESPEK"
        nsites = seq.partition('?')[2].count('S') + seq.partition('?')[2].count('T')
        pf = ProForma.parse(seq)
        for include_unmodified in [False, True]:
            with self.subTest(include_unmodified=include_unmodified):
                proteoforms = list(pf.proteoforms(include_unmodified=include_unmodified))
                if not include_unmodified:
                    self.assertEqual(len(proteoforms), math.comb(nsites, k))   # Phospho on T or S, exactly `k` times
                else:
                    self.assertEqual(
                        len(proteoforms),
                        sum([math.comb(nsites, i) for i in range(k + 1)]),  # Phospho on T or S, anywhere from 0 to `k` times
                    )

    def test_localization_tag(self):
        seq = "EMEVT[#g1]S[#g1]ES[Phospho#g1]PEK"
        nsites = seq.count('#g1')
        pf = ProForma.parse(seq)
        for include_unmodified in [False, True]:
            with self.subTest(include_unmodified=include_unmodified):
                proteoforms = list(pf.proteoforms(include_unmodified=include_unmodified))
                self.assertEqual(len(proteoforms), nsites + include_unmodified)

    def test_unlocalized_modification(self):
        seq = "[Phospho]?EMEVTSESPEK"
        pf = ProForma.parse(seq)
        for include_unmodified in [False, True]:
            with self.subTest(include_unmodified=include_unmodified):
                proteoforms = list(pf.proteoforms(include_unmodified=include_unmodified))
                self.assertEqual(len(proteoforms), len(pf) + include_unmodified)

    def test_comup_stacking(self):
        k = 2  # number of modifications to combine
        limit = 2  # stack limit
        seq = f"[Phospho|Position:S|Position:T|comup|Limit:{limit}]^{k}?EMEVTESPEK"
        nsites = seq.partition('?')[2].count('S') + seq.partition('?')[2].count('T')
        self.assertGreaterEqual(nsites * limit, k)  # otherwise we can't place `k` mods even with stacking
        effective_limit = min(limit, k)  # if limit >= k, then we can just treat it as a normal combinatorial expansion
        pf = ProForma.parse(seq)
        proteoforms = list(pf.proteoforms())
        self.assertEqual(len(proteoforms), math.comb(nsites + effective_limit - 1, k))  # number of ways to place `k` indistinguishable mods on `nsites` distinguishable sites with a stack limit of `effective_limit`
        proteoforms = list(pf.proteoforms(True))
        self.assertEqual(len(proteoforms), sum([math.comb(nsites + min(limit, i) - 1, i) for i in range(k + 1)]))  # number of ways to place anywhere from 0 to `k` indistinguishable mods on `nsites` distinguishable sites with a stack limit of `effective_limit`

    def test_labile(self):
        seq = "{Phospho}EMEVTESPEK"
        pf = ProForma.parse(seq)
        proteoforms = list(pf.proteoforms(False, True))
        self.assertEqual(len(proteoforms), len(pf) + 1)  # all possible sites and the form where phospho is kept as labile


class ProteoformsFunctionTest(unittest.TestCase):
    def test_proteoforms(self):
        seq = "EMEV(TS)[Phospho]ESPEK"
        nsites = 2  # length of the range
        pf = ProForma.parse(seq)
        for include_unmodified in [False, True]:
            for deepcopy in [False, True]:
                with self.subTest(include_unmodified=include_unmodified, deepcopy=deepcopy):
                    forms = list(proteoforms(pf, include_unmodified=include_unmodified, deepcopy=deepcopy))
                    self.assertEqual(len(forms), nsites + include_unmodified)   # Phospho on T or S (+ no phospho if include_unmodified)

    def test_coerce_modification(self):
        for s, m in [("Phospho", GenericModification("Phospho")),
                     ("UNIMOD:21", UnimodModification("21")),
                     ("MOD:00046", PSIModModification("00046"))]:
            with self.subTest(s=s):
                self.assertEqual(_coerce_string_to_modification(s), m)

    def test_modification_target_from_str(self):
        for s, t in [("S", ModificationTarget('S')),
                     ("T", ModificationTarget('T')),
                     ("N-term", ModificationTarget(None, True, False)),
                     ("C-term", ModificationTarget(None, False, True)),
                     ("N-term:K", ModificationTarget('K', True, False)),
                     ("C-term:Y", ModificationTarget('Y', False, True))]:
            with self.subTest(s=s):
                self.assertEqual(ModificationTarget.from_str(s), t)

    def test_from_simple_dict(self):
        seq = "EMEVTSESPEK"
        variable_mods = {"Phospho": ["S", "T"]}
        nsites = seq.count("S") + seq.count("T")
        pf = ProForma.parse(seq)
        for include_unmodified in [False, True]:
            for deepcopy in [False, True]:
                with self.subTest(include_unmodified=include_unmodified, deepcopy=deepcopy):
                    forms = list(proteoforms(pf, variable_modifications=variable_mods, include_unmodified=include_unmodified, deepcopy=deepcopy))
                    if include_unmodified:
                        self.assertEqual(len(forms), nsites + 1)   # Phospho on T or S + no phospho
                    else:
                        self.assertEqual(len(forms), nsites)  # Phospho on T or S

        forms = list(proteoforms(pf, variable_modifications=variable_mods, expand_rules=True))
        self.assertEqual(len(forms), 2 ** nsites)  # all combinations of phospho / no phospho on each S or T

    def test_expand(self):
        seq = "EMEVTSESPEK"
        total_sites = seq.count("S") + seq.count("T") + seq.count("M")
        oxidation_sites = seq.count("M")
        variable_mods = {"Phospho": ["S", "T"], "Oxidation": ["M"]}
        pf = ProForma.parse(seq)
        combos = peptidoforms(
            pf,
            variable_modifications=variable_mods,
            expand_rules=True,
        )
        variants = list(combos)
        self.assertEqual(len(variants), 2 ** total_sites)  # all combinations of phospho on S/T and oxidation on M
        self.assertAlmostEqual(
            (2 ** oxidation_sites) / (2 ** oxidation_sites - 1),
            len(variants) / sum(['Oxidation' in str(p) for p in variants])
        )

    def test_from_str(self):
        seq = "EMEVTSESPEK"
        variable_mods = ["Phospho|Position:S|Position:T"]
        nsites = seq.count("S") + seq.count("T")
        pf = ProForma.parse(seq)
        for include_unmodified in [False, True]:
            with self.subTest(include_unmodified=include_unmodified):
                forms = list(proteoforms(pf, variable_modifications=variable_mods, include_unmodified=include_unmodified))
                if include_unmodified:
                    self.assertEqual(len(forms), nsites + 1)   # Phospho on T or S + no phospho
                else:
                    self.assertEqual(len(forms), nsites)  # Phospho on T or S
        forms = list(proteoforms(pf, variable_modifications=variable_mods, expand_rules=True))
        self.assertEqual(len(forms), 2 ** nsites)  # all combinations of phospho / no phospho on each S or T

    def test_expand_mods_from_list(self):
        seq = "EMEVTSESPEK"
        variable_mods = ["Phospho|Position:S", "Phospho|Position:T"]
        nsites = seq.count("S") + seq.count("T")
        pf = ProForma.parse(seq)
        forms = list(proteoforms(pf, variable_modifications=variable_mods, expand_rules=True))
        self.assertEqual(len(forms), 2 ** nsites)  # all combinations of phospho on 0, 1, or 2 of the S or T

    def test_expand_mods_from_dict(self):
        seq = "EMEVTSESPEK"
        variable_mods = {"Phospho": ["S", "T"]}
        nsites = seq.count("S") + seq.count("T")
        pf = ProForma.parse(seq)
        forms = list(proteoforms(pf, variable_modifications=variable_mods, expand_rules=True))
        self.assertEqual(len(forms), 2 ** nsites)  # all combinations of phospho on 0, 1, or 2 of the S or T

    def test_expand_from_dict(self):
        seq = "EMEVTSESPEK"
        variable_mods = {"Phospho": ["S", "T"]}
        nsites = seq.count("S") + seq.count("T")
        pf = ProForma.parse(seq)
        forms = list(proteoforms(pf, variable_modifications=variable_mods, expand_rules=True))
        self.assertEqual(len(forms), 2 ** nsites)  # all combinations of phospho on 0, 1, or 2 of the S or T

    def test_expand_mods_comup(self):
        seq = "EMEVTSESPEK"
        limit = 2
        variable_mods = [f"Phospho|Position:S|Position:T|comup|Limit:{limit}"]
        nsites = seq.count("S") + seq.count("T")
        pf = ProForma.parse(seq)
        forms = list(proteoforms(pf, variable_modifications=variable_mods, expand_rules=True))
        self.assertEqual(len(forms), (limit + 1) ** nsites)  # all combinations of 0 to `limit` phosphos on each S or T


if __name__ == '__main__':
    unittest.main()
