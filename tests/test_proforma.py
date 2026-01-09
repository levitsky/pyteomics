from os import path
import unittest
import pickle
import pyteomics
pyteomics.__path__ = [path.abspath(
    path.join(path.dirname(__file__), path.pardir, 'pyteomics'))]
from pyteomics.proforma import (
    PSIModModification, ProForma, TaggedInterval, parse, MassModification, ProFormaError, TagTypeEnum,
    ModificationRule, StableIsotope, GenericModification, Composition, to_proforma, ModificationMassNotFoundError,
    AdductParser, ChargeState,
    std_aa_comp, obo_cache, process_tag_tokens)


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

    def test_charge_adducts(self):
        sequences = ['PEPTIDE/1[+2Na+,-H+]', 'PEPTIDE/-1[+e-]', 'PEPTIDE/1[+2H+,+e-]']
        charges = [1, -1, 1]
        adducts_list = [[('Na', 1, 2), ('H', 1, -1)], [('e-', -1, 1)], [('H', 1, 2), ('e-', -1, 1)]]
        for seq, charge, adducts in zip(sequences, charges, adducts_list):
            i = ProForma.parse(seq)
            self.assertEqual(i.charge_state.charge, charge)
            self.assertEqual(i.charge_state.adducts, adducts)

    def test_composition_with_adducts(self):
        sequences = ['PEPTIDE/1[+2Na+,-H+]', 'PEPTIDE/-1[+e-]', 'PEPTIDE/1[+2H+,+e-]', 'PEPTIDE', 'PEPTIDE/1']
        neutral_comp = Composition(sequence='PEPTIDE')
        adducts_list = [Composition({'Na': 2, 'H': -1}),
                        Composition({'e-': 1}),
                        Composition({'H': 2, 'e-': 1}),
                        Composition({}),
                        Composition({'H': 1})]
        for seq, adducts in zip(sequences, adducts_list):
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
            "EM[R: Methionine sulfone]EVEES[O-phospho-L-serine]PEK",
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
            # "EMEVEESPEK+ELVISLIVER",
            # "EMEVEESPEK/2+ELVISLIVER/3",
            # "A[X:DSS#XL1]//B[#XL1]+C[X:DSS#XL1]//D[#XL1]",
            "<[Carbamidomethyl]@C>ATPEILTCNSIGCLK",
            "<[Oxidation]@C,M>MTPEILTCNSIGCLK",
            "<[TMT6plex]@K,N-term>ATPEILTCNSIGCLK",
            "<[TMT6plex]@K,N-term:A>ATPEILTCNSIGCLK",
            "<[TMT6plex]@K,N-term:A,N-term:B>ATPEILTCNSIGCLK",
            "EM[Oxidation]EVEES[Phospho]PEK",
            "EM[L-methionine sulfoxide]EVEES[O-phospho-L-serine]PEK",
            "EM[R: L-methionine sulfone]EVEES[O-phospho-L-serine]PEK",
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
            "EM[RESID:AA0581]EVEES[RESID:AA0037]PEK",
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
            # "EMEVEESPEK/2+ELVISLIVER/3",
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
            # "(>Trypsin)AANSIPYQVSLNS+(>Keratin)AKEQFERQTA",
            # "(>P07225 Vitamin K-dependent protein S OS=Homo sapiens OX=9606 GN=PROS1 PE=1 (SV=1) RANGE=12..42)GGK[xlink:dss[138]#XLDSS]IEVQLK//(>P07225 Vitamin K-dependent protein S OS=Homo sapiens OX=9606 GN=PROS1 PE=1 SV=1)KVESELIK[#XLDSS]PINPR/4",
            # "(>>>Trastuzumab Fab and coeluting Fc)(>>Fab)(>Heavy chain)EVQLVESGGGLVQPGGSLRLSC[M:l-cystine (cross-link)#XL1]AASGFNIKDTYIHWVRQAPGKGLEWVARIYPTNGYTRYADSVKGRFTISADTSKNTAYLQMNSLRAEDTAVYYC[#XL1]SRWGGDGFYAMDYWGQGTLVTVSSASTKGPSVFPLAPSSKSTSGGTAALGC[M:l-cystine (cross-link)#XL2]LVKDYFPEPVTVSWNSGALTSGVHTFPAVLQSSGLYSLSSVVTVPSSSLGTQTYIC[#XL2]NVNHKPSNTKVDKKVEPKSC[M:l-cystine (cross-link)#XL3]DKT//(>Light chain)DIQMTQSPSSLSASVGDRVTITC[M:l-cystine (cross-link)#XL4]RASQDVNTAVAWYQQKPGKAPKLLIYSASFLYSGVPSRFSGSRSGTDFTLTISSLQPEDFATYYC[#XL4]QQHYTTPPTFGQGTKVEIKRTVAAPSVFIFPPSDEQLKSGTASVVC[M:l-cystine (cross-link)#XL5]LLNNFYPREAKVQWKVDNALQSGNSQESVTEQDSKDSTYSLSSTLTLSKADYEKHKVYAC[#XL5]EVTHQGLSSPVTKSFNRGEC[#XL3]+(>Fc)HTCPPCPAPELLGGPSVFLFPPKPKDTLMISRTPEVTCVVVDVSHEDPEVKFNWYVDGVEVHNAKTKPREEQYNSTYRVVSVLTVLHQDWLNGKEYKCKVSNKALPAPIEKTISKAKGQPREPQVYTLPPSREEMTKNQVSLTCLVKGFYPSDIAVEWESNGQPENNYKTTPPVLDSDGSFFLYSKLTVDKSRWQQGNVFSCSVMHEALHNHYTQKSLSLSPGK",
        ]
        for seq in positive:
            parsed = ProForma.parse(seq)
            assert parsed is not None

    def test_nonstandard_amino_acid(self):
        seq = ProForma.parse("PEPTX[MOD:01001]IDE")
        bad_seq = ProForma.parse("PEPTXIDE")
        assert seq.mass != bad_seq.mass
        self.assertAlmostEqual(seq.mass, 884.4127280267099, 4)


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
