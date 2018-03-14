from os import path
import pyteomics
pyteomics.__path__ = [path.abspath(path.join(path.dirname(__file__), path.pardir, 'pyteomics'))]
import unittest
from pyteomics.mzid import *
from pyteomics import auxiliary as aux
from data import mzid_spectra
from itertools import product
from io import BytesIO

class MzidTest(unittest.TestCase):
    maxDiff = None
    def testReadPSM(self):
        for rec, refs, rs, it, ui in product((True, False), repeat=5):
            for func in [MzIdentML, read, chain,
                    lambda x, **kw: chain.from_iterable([x], **kw)]:
                with func('test.mzid', recursive=rec, retrieve_refs=refs,
                        read_schema=rs, iterative=it, use_index=ui) as reader:
                    psms = list(reader)
                    self.assertEqual(psms, mzid_spectra[(rec, refs)])

    def test_unit_info(self):
        with MzIdentML('test.mzid') as handle:
            for protocol in handle.iterfind("SpectrumIdentificationProtocol"):
                fragment_tolerance = protocol['FragmentTolerance']
                self.assertEqual(fragment_tolerance['search tolerance minus value'].unit_info, 'dalton')
                parent_tolerance = protocol['ParentTolerance']
                self.assertEqual(parent_tolerance['search tolerance plus value'].unit_info, 'parts per million')

    def test_structure_normalization(self):
        text = '''
        <Inputs>
          <SearchDatabase id="UniProt:Human" location="UniProt:Human:2014_02.fasta" version="2014_02">
            <DatabaseName><cvParam cvRef="MS" accession="MS:1001013" name="database name" value="UniProt"/></DatabaseName>
            <cvParam cvRef="MS" accession="MS:1001468" name="taxonomy: common name" value="Human"/>
            <cvParam cvRef="MS" accession="MS:1001254" name="DB source UniProt" value=""/>
            <cvParam cvRef="MS" accession="MS:1001015" name="database original uri" value="http://www.uniprot.org/uniprot/?query=taxonomy%3a9606+AND+keyword%3a1185&amp;force=yes&amp;format=fasta&amp;include=yes"/>
          </SearchDatabase>
          <SearchDatabase id="RefSeq:Human" location="RefSeq:Human:r63.fasta" version="r63">
            <DatabaseName><cvParam cvRef="MS" accession="MS:1001013" name="database name" value="RefSeq"/></DatabaseName>
            <cvParam cvRef="MS" accession="MS:1001468" name="taxonomy: common name" value="Human"/>
            <cvParam cvRef="MS" accession="MS:1001253" name="DB source NCBI" value=""/>
            <cvParam cvRef="MS" accession="MS:1001015" name="database original uri" value="ftp://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/mRNA_Prot/human.protein.faa.gz"/>
          </SearchDatabase>
          <SearchDatabase id="RefSeq:Mouse" location="RefSeq:Mouse:r63.fasta" version="r63">
            <DatabaseName><cvParam cvRef="MS" accession="MS:1001013" name="database name" value="RefSeq"/></DatabaseName>
            <cvParam cvRef="MS" accession="MS:1001468" name="taxonomy: common name" value="Mouse"/>
            <cvParam cvRef="MS" accession="MS:1001253" name="DB source NCBI" value=""/>
            <cvParam cvRef="MS" accession="MS:1001015" name="database original uri" value="ftp://ftp.ncbi.nlm.nih.gov/refseq/M_musculus/mRNA_Prot/mouse.protein.faa.gz"/>
          </SearchDatabase>
          <SearchDatabase id="UniProt:Mouse" location="UniProt:Mouse:2014_02.fasta" version="2014_02">
            <DatabaseName><cvParam cvRef="MS" accession="MS:1001013" name="database name" value="UniProt"/></DatabaseName>
            <cvParam cvRef="MS" accession="MS:1001468" name="taxonomy: common name" value="Mouse"/>
            <cvParam cvRef="MS" accession="MS:1001254" name="DB source UniProt" value=""/>
            <cvParam cvRef="MS" accession="MS:1001015" name="database original uri" value="http://www.uniprot.org/uniprot/?query=taxonomy%3a10090+AND+keyword%3a1185&amp;force=yes&amp;format=fasta&amp;include=yes"/>
          </SearchDatabase>
          <SearchDatabase id="RefSeq" location=".">
            <DatabaseName><cvParam cvRef="MS" accession="MS:1001013" name="database name" value="RefSeq"/></DatabaseName>
            <cvParam cvRef="MS" accession="MS:1001253" name="DB source NCBI" value=""/>
          </SearchDatabase>
          <SpectraData id="CPTAC_CompRef_00_iTRAQ_01_2Feb12_Cougar_11-10-09" location="CPTAC_CompRef_00_iTRAQ_01_2Feb12_Cougar_11-10-09.mzML">
            <FileFormat><cvParam cvRef="MS" accession="MS:1000584" name="mzML file" value=""/></FileFormat>
            <SpectrumIDFormat><cvParam cvRef="MS" accession="MS:1000768" name="Thermo nativeID format" value=""/></SpectrumIDFormat>
          </SpectraData>
        </Inputs>
        '''
        buff = BytesIO(bytes(text, 'utf8'))
        datum = next(read(buff).iterfind("SpectraData"))
        index = aux.cvquery(datum)
        assert index['MS:1000768'] == 'Thermo nativeID format'
        text = '''
        <Inputs>
            <SourceFile id="Scaffold_Spectrum_SOURCEFILE_OT_141126_16_GFP (F004691)_Mascot" location="D:\inetpub\mascot\data\20141217\F004691.dat">
                <FileFormat>
                    <cvParam accession="MS:1001199" name="Mascot DAT file" cvRef="PSI-MS"/>
                </FileFormat>
            </SourceFile>
            <SearchDatabase id="unreferenced database" location="" name="ArabidopsisTAIR10_20101214.fasta" numDatabaseSequences="35545" version="Unknown">
                <FileFormat>
                    <cvParam accession="MS:1001348" name="FASTA format" cvRef="PSI-MS"/>
                </FileFormat>
                <DatabaseName>
                    <userParam name="name" value="ArabidopsisTAIR10_20101214.fasta"/>
                </DatabaseName>
                <cvParam accession="MS:1001073" name="database type amino acid" cvRef="PSI-MS"/>
            </SearchDatabase>
            <SpectraData id="OT_141126_16_GFP (F004691)" location="Y:\A-J\AJ\will data 141015\will scaff\141216_wh_gfp mzIdent export 05-Mar-2015 08-58-45-AM\gfp\OT_141126_16_GFP (F004691).mzid_OT_141126_16_GFP_(F004691).MGF">
                <FileFormat>
                    <cvParam accession="MS:1001062" name="Mascot MGF file" cvRef="PSI-MS"/>
                </FileFormat>
                <SpectrumIDFormat>
                    <cvParam accession="MS:1000774" name="multiple peak list nativeID format" cvRef="PSI-MS"/>
                </SpectrumIDFormat>
            </SpectraData>
        </Inputs>
        '''
        buff = BytesIO(bytes(text, 'utf8'))
        datum = next(read(buff).iterfind("SpectraData"))
        index = aux.cvquery(datum)
        assert index['MS:1000774'] == 'multiple peak list nativeID format'


if __name__ == '__main__':
    unittest.main()
