import os
import shutil
import tempfile
import unittest
import json

from pathfinderlib.main import main, get_config, Aligner, Assigner

class PathfinderTest(unittest.TestCase):
    def setUp(self):
        self.config = get_config(None) # get user specific config
        self.output_dir = tempfile.mkdtemp()
        self.summary_fp = os.path.join(self.output_dir, "summary.txt")

    def tearDown(self):
        shutil.rmtree(self.output_dir)

    def ko_fp_from_aln(self, r1_name):
        return os.path.join(self.output_dir, os.path.basename(os.path.splitext(r1_name)[0]+'.ko'))

    def run_pipeline(self, R1, R2, config_file):
        r1 = tempfile.NamedTemporaryFile(suffix=".fastq")
        r1.write(R1)
        r1.seek(0)
        r2 = tempfile.NamedTemporaryFile(suffix=".fastq")
        r2.write(R2)
        r2.seek(0)
        
        args = [
            "--forward-reads", r1.name,
            "--reverse-reads", r2.name,
            "--summary-file", self.summary_fp,
            "--output-dir", self.output_dir,
            "--config-file", config_file.name
            ]

        main(args)
        return r1.name

    def check_results(self, results_fp, summary_fp, expected_results, expected_summary):
        #check if the output is correct
        observed = open(results_fp).read()
        self.assertEqual(observed.strip(), expected_results)
        
        # check if the summary file is correct
        with open(self.summary_fp) as f:
            observed = json.load(f)
        self.assertEqual(observed.get('data', {}), expected_summary)
        
    def test_main_RAP_bestHit(self):
        #generate the correct config file
        self.config['mapping_method'] = 'best_hit'
        config_file = tempfile.NamedTemporaryFile(suffix=".json")
        json.dump(self.config, config_file)
        config_file.seek(0)
        
        # run the pipeline and check results
        r1_name = self.run_pipeline(MOCK_R1, MOCK_R2, config_file)
        self.check_results(self.ko_fp_from_aln(r1_name), self.summary_fp, EXPECTED_OUTPUT, EXPECTED_SUMMARY)

    def test_empty_alignment(self):
        config_file = tempfile.NamedTemporaryFile(suffix=".json")
        json.dump(self.config, config_file)
        config_file.seek(0)

        # run the pipeline and check results
        r1_name = self.run_pipeline(MOCK_NO_MATCH_R1, MOCK_NO_MATCH_R2, config_file)
        self.check_results(self.ko_fp_from_aln(r1_name), self.summary_fp, EXPECTED_NO_MATCH_OUTPUT, EXPECTED_NO_MATCH_SUMMARY)

    def test_make_index_bestHit(self):
        # create the mock kegg file
        kegg = tempfile.NamedTemporaryFile()
        kegg.write(MOCK_KEGG)
        kegg.seek(0)
        
        kegg2ko = tempfile.NamedTemporaryFile()
        self.config['kegg_fp'] = kegg.name
        self.config['kegg_to_ko_fp'] = kegg2ko.name

        # make index
        assignerApp = Assigner(self.config)
        assignerApp.make_index()

        # check equal
        observed = open(self.config["kegg_to_ko_fp"]).read()
        self.assertEqual(observed.strip(), KEGG_TO_KO.strip())


# get known sequences from KEGG to test the cases:
# a sequence having no hits in the kegg database (s0)
# Multiple proteins mapping to the same KO (s3-6)
# Single sequence having multiple KO values (s7)
# Multiple sequences sharing multiple KO values (s1, s2)
# R1 and R2 having exclusive KOs (s8 and s9)
MOCK = """\
>s0_nomatch
GCTATAAGCAGTAACGGCGGTGATGTTTTAACAAACAACGGAGTCAATAAAGTCATCAGGTAGTACTGCATCATCAAGAAGAGGATATTTAACT
+
EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
>s1_ptr_452965
cgctggccagtccctaagcacgtgggttgggttgtcctgcttggctgcggagggagtggaacctcgatattggtggtgtccatcgtgggcagcggactaataaaggccatggcgccagcagaaat
+
EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
>s2_pap_PSPA7_3507
caacgcgtgaccgagcgccgccagcaaggcctgcgcgttcccggcctggcggtgatcctggtgggcaccgatccggcctctcaggtctatgtggcgcacaagcgcaaggactgcgaggaagtcgg
+
EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
>s3_mcc_693781
ttcaagcgcaaagttgggggcctgggattcctggtgaaggagcgggtcagtaagccacccgtgatcatctccgacctgattcgtgggggcgccgcggaacagagtggcctcatccaggccggaga
+
EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
>s4_ecb_791246
caggatgaccccaagggtcacagcctcggcaagcacaggaatgagtccctgcagcccgtcaccggaatggcaaagaagtctccagaatccctggtcaagctggatgtgcccccctcggcctgccc
+
EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
>s5_ecb_791246
caggatgaccccaagggtcacagcctcggcaagcacaggaatgagtccctgcagcccgtcaccggaatggcaaagaagtctccagaatccctggtcaagctggatgtgcccccctcggcctgccc
+
EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
>s6_cfa_403784
gcaccggccctggcgcccccaccctcgcccccaccagcaccagaccacagcagccccccactcacccggcctccagatgggcccaagttccctcgtgtgaagaactgggaggtggggagcatcac
+
EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
>s7_mmu_14874
atggagtacctggaagagactcggcctatcccacggctcctgcctcaggacccacagaaaagagccatcgtgcgcatgatttctgacctcatcgctagtggcatccagccccttcagaacctgtc
+
EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
"""

MOCK_R1 = MOCK + """\
>s8_bta_617935
ctactgactctgtctgcacttacagacagacctaaattgcctgataactacactcaggacacctggcagaagcttcacgaagcggtgagagccatacagagcagcacctccatccgctacaacct
+
EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
"""

MOCK_R2 = MOCK + """\
>s9_bta_505481
ctgcaagttcccctctttataatattcattctcatttacttggtcaatgtggttggaaacgtgggcatcatcctgctggtgctcctggactcgcatctccacacgcccatgtactttttcctcag
taacctgtctctggtggactttggttactccacagctgtcattcccacagtcctggctggattcctgaaaggccgcggggtcatctcctataatgtgtgtgctgctcagatgttc
+
EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
"""

EXPECTED_OUTPUT = """\
KO	ko_abundance
K00288 methylenetetrahydrofolate dehydrogenase (NADP+) [EC:1.5.1.5]	1.5
K00295 methylenetetrahydrofolate dehydrogenase (NAD+) [EC:1.5.1.15]	0.5
K00491 nitric-oxide synthase [EC:1.14.13.39]	8.0
K00799 glutathione S-transferase [EC:2.5.1.18]	1.0
K01491 methenyltetrahydrofolate cyclohydrolase [EC:3.5.4.9]	1.5
K01800 maleylacetoacetate isomerase [EC:5.2.1.2]	1.0
K01938 formate--tetrahydrofolate ligase [EC:6.3.4.3]	0.5
K10609 cullin 4	1.0
K04257 olfactory receptor	1.0"""

EXPECTED_SUMMARY = {"unique_prot_hits": 14, "mapped_sequences": 16, "ko_hits": 14, "unique_ko_hits": 16, "mapped_sequences_evalue": 16}


MOCK_NO_MATCH_R1 = """\
@HWI-D00727:9:C6JHHANXX:2:1101:4623:2240 1:N:0:CGTACTAGCTAAGCCT
TTACTACCTTATATAACGAAAACACAAATGTAAAACT
+
:@BB@FGGGGGGDFGGGGGGGGGGGGGGGGGCGGGGG
@HWI-D00727:9:C6JHHANXX:2:1101:6139:2116 1:N:0:CGTACTAGCTAAGCCT
CGATCGCTTCCTTCTCAGTCACGATCTCTGTTTCCCAGTCACGATCTATTATTCAGTCACGATCTCTCGCATAAGACCCCGATCG
+
A3@@BCGGGGGGGGGGGGGGGGGGGGGGGGGEGGGGGGGGGGGGFGGGGGGGGGGGGGGGGGGGGGGGGGGGGGEGGGGGGGGGG
@HWI-D00727:9:C6JHHANXX:2:1101:7152:2193 1:N:0:CGTACTAGCTAAGCCT
GTTTCCCACTTAGCAATATTTAGGGACCTTAGCTGGCGGTCTGGGTTGTTTCCCTCTTGACACCGGACGTTAGCACCCGATGTCTGTCTCCCGTGATTGCACTCTTCGGTATTCGGAGTTTGCTAT
+
3@B=@EFGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG@DGGGGGG<EGGGGCGGGGGGGGGFGG@AGGGGEGGGGEGGCDGGGGEGGGGGGGGGGDF=GGGGGGEGGGGEGBEGED:DDG8EGGG
@HWI-D00727:9:C6JHHANXX:2:1101:8281:2188 1:N:0:CGTACTAGCTAAGCCT
CTTCCCACTTCGTTTCCCACTTAGCAATATTTAGGGACCTTAGCTGGCGGTCTGGGTTGTTTCCCTCTTGACACCGGACGTTAGCACCCGATGTCTGTCTCCCGTGATTGCACTCTTCGGTATTCG
+
BBB@BFBFGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG>FGGGDG/EDGGBFGGGGGG.9FFGGGGGGGGFADB.:DGG//@@9DGGDC.D=EG
@HWI-D00727:9:C6JHHANXX:2:1101:12300:2173 1:N:0:CGTACTAGCTAAGCCT
CTTCTGTACGCAGTCACGATCTTCTGTTTCCCAGTCACGATCCGTAAATTAACCCCCCGATCGTGACTGCGTGTTGAGATCGTGAC
+
BCBBCEFGGGGGGGGGGGGGGGGGGGGFGGFGGGGGGGGCGGGGGGFGGGGGGGGGGGGGGGGGGGFGGGGGGGGGGGGGGGGGGG
"""

MOCK_NO_MATCH_R2 = """\
@HWI-D00727:9:C6JHHANXX:2:1101:4623:2240 2:N:0:CGTACTAGCTAAGCCT
AGTTTTACATTTGTGTTTTCGTTATATAAGGTAGTAA
+
BCCCBGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG
@HWI-D00727:9:C6JHHANXX:2:1101:6139:2116 2:N:0:CGTACTAGCTAAGCCT
CGATCGGGGTCTTATGCGAGAGATCGTGACTGAATAATAGATCGTGACTGGGAAACAGAGATCGTGACTGAGAAGGAAGCGATCG
+
BCCCBGGGGEGGGGGGGGGGGGGGGGGGGGGGGGGFGGGGGGGGGGGGGGGGGGGGGGGGGGGGGEGGGGEGGEGGGGGGGGGGG
@HWI-D00727:9:C6JHHANXX:2:1101:7152:2193 2:N:0:CGTACTAGCTAAGCCT
GTCATGGTTGGGGGGTCTATTGCTGATTACCCCGCCATAGCAAACTCCGAATACCGAAGAGTGCAATCACGGGAGACAGACATCGGGTGCTAACGTCCGGTGTCAAGAGGGAAACAACCCAGACCG
+
ABB@BGGGGGGGFG@BFGGGGFGGGGG=DFGGFGGGGGGGEG=GGGGGBGGGGGGGGGDGGGGGGGGGGGGDGGCGGGGGGGGGGDGGGGGGGGBDGGDDCGBG=DD@DDCGBGGGGGDBG.BGG:
@HWI-D00727:9:C6JHHANXX:2:1101:8281:2188 2:N:0:CGTACTAGCTAAGCCT
GGGTAGAGCACTGTCATGGTTGGGGGGTCTATTGCTGATTACCCCGCCATAGCAAACTCCGAATACCGAAGAGTGCAATCACGGGAGACAGACATCGGGTGCTAACGTCCGGTGTCAAGAGGGAAA
+
BBBCCGCGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGDGGGGGCDGGGGGGGGGGGGGGGGGGGGGGGGEGGGGGGGGDGGGEGGGGGGGGB
@HWI-D00727:9:C6JHHANXX:2:1101:12300:2173 2:N:0:CGTACTAGCTAAGCCT
GTCACGATCTCAACACGCAGTCACGATCGGGGGGTTAATTTACGGATCGTGACTGGGAAACAGAAGATCGTGACTGCGTACAGAAG
+
BBCBBGGGEGGGGGGGGBGGFCGGGDGGGFGGGG:/9CEGG@GCGD<EGBD@EDGGG<DB0;@CDGGEGGGGG<@DDDGCADED@@
"""

EXPECTED_NO_MATCH_OUTPUT = "KO	ko_abundance"

EXPECTED_NO_MATCH_SUMMARY = {"unique_prot_hits": 0, "mapped_sequences": 0, "ko_hits": 0, "unique_ko_hits": 0, "mapped_sequences_evalue": 0}

# data to test the kegg2ko index code. Tests for:
# Protein to a single K0 prefix
# Protein to a single K1 prefix
# Protein to a single K0 with  no description
# Protein to multiple KOs
# Protein with no KO associations
MOCK_KEGG = """\
>hsa:1145  CHRNE; cholinergic receptor, nicotinic, epsilon ; K04817 cholinergic receptor, nicotinic, epsilon
MARAPLGVLLLLGLLGRGVGKNEELRLYHHLFNNYDPGSRPVREPEDTVTISLKVTLTNL
>hsa:27284  SULT1B1; sulfotransferase family, cytosolic, 1B, member 1 (EC:2.8.2.-); K01025
MLSPKDILRKDLKLVHGYPMTCAFASNWEKIEQFHSRPDDIVIATYPKSGTTWVSEIIDM
>hsa:3797  KIF3C; kinesin family member 3C ; K10394 kinesin family member 3/17
MASKTKASEALKVVARCRPLSRKEEAAGHEQILTMDVKLGQVTLRNPRAAPGELPKTFTF
>hsa:223  ALDH9A1; aldehyde dehydrogenase 9 family, member A1 (EC:1.2.1.3 1.2.1.19 1.2.1.47); K00128 aldehyde dehydrogenase (NAD+) [EC:1.2.1.3]; K00137 aminobutyraldehyde dehydrogenase [EC:1.2.1.19]; K00149 4-trimethylammoniobutyraldehyde dehydrogenase [EC:1.2.1.47]
MFLRAGLAALSPLLRSLRPSPVAAMSTGTFVVSQPLNYRGGARVEPADASGTEKAFEPAT
>foo:001; FOO; spam spam bacon K0 spam
SPAMSPAMSPAMSPAMSPAMEGGSSPAMBACON
"""

KEGG_TO_KO = """\
Subject	KO
hsa:1145	K04817 cholinergic receptor, nicotinic, epsilon
hsa:27284	K01025
hsa:3797	K10394 kinesin family member 3/17
hsa:223	K00128 aldehyde dehydrogenase (NAD+) [EC:1.2.1.3]
hsa:223	K00137 aminobutyraldehyde dehydrogenase [EC:1.2.1.19]
hsa:223	K00149 4-trimethylammoniobutyraldehyde dehydrogenase [EC:1.2.1.47]"""
