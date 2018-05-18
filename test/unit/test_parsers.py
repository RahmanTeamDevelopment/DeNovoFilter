"""Unit tests for the parsers module"""

from unittest import TestCase
from denovo_ import parsers



class TestVCFParser(TestCase):


    def test_parse_vcf_record_single_alt_single_transcript(self):

        vcf_record = '11      1642299 .       T       A       200     PASS    ABPV=1.00e+00;FPV=5.41e-02;FR=0.5000;' \
                     'HP=1;MMLQ=34;NF=14;NR=13;PP=200;RMP=57.26;RPV=9.05e-01;SC=CCCACCTTGTTGCAGGTGGGA;TC=60;TCF=39;' \
                     'TCR=21;TR=27;TYPE=Substitution;TRANSCRIPT=ENST00000399682;GENE=KRTAP5-4;GENEID=ENSG00000241598;' \
                     'TRINFO=-/1.2kb/1/1.2kb/228;LOC=3UTR;CSN=c.+338A>T;PROTPOS=.;PROTREF=.;PROTALT=.;CLASS=3PU;' \
                     'SO=3_prime_UTR_variant;IMPACT=3;ALTANN=.;ALTCLASS=.;ALTSO=.'

        expected = {}
        expected[('11', '1642299', 'T', 'A')] = [
            {
                'quality': 'high',
                'TR': 27,
                'TC': 60,
                'NF': 14,
                'NR': 13,
                'gene': 'KRTAP5-4',
                'csn': 'c.+338A>T',
                'class_': '3PU',
                'altann': '.',
                'altclass': '.'
            }
        ]

        parsed_record = parsers.parse_vcf_record(vcf_record)
        self.assertDictEqual(parsed_record, expected)


    def test_parse_vcf_record_multiple_alt_single_transcript(self):

        vcf_record = 'chrX       3523053 .       GACACACACAC     GAC,G   200     PASS    ABPV=9.76e-02;FPV=3.81e-06;' \
                     'FR=0.5000,0.5000;HP=3;NF=0,3;NR=5,3;PP=200,200;RMP=66.35;RPV=7.17e-02;SC=CCTTCACCAAGACACACACAC;' \
                     'TC=35;TCF=18;TCR=17;TR=5,6;TYPE=Deletion,Deletion;TRANSCRIPT=ENST00000262848,ENST00000262848;' \
                     'GENE=PRKX,PRKX;GENEID=ENSG00000183943,ENSG00000183943;TRINFO=-/109.2kb/9/6.0kb/358,-/109.2kb/9/' \
                     '6.0kb/358;LOC=3UTR,3UTR;CSN=c.+3949_+3956del8,c.+3947_+3956del10;PROTPOS=.,.;PROTREF=.,.;' \
                     'PROTALT=.,.;CLASS=3PU,3PU;SO=3_prime_UTR_variant,3_prime_UTR_variant;IMPACT=3,3;ALTANN=c.+3909_+' \
                     '3916del8,c.+3909_+3918del10;ALTCLASS=.,.;ALTSO=.,.        GT:GL:GQ:NR     1/2:-1,-1,-1:100:10922'

        expected = {}
        expected[('X', '3523053', 'GACACACACAC', 'GAC')] = [
            {
                'quality': 'low',
                'TR': 5,
                'TC': 35,
                'NF': 0,
                'NR': 5,
                'gene': 'PRKX',
                'csn': 'c.+3949_+3956del8',
                'class_': '3PU',
                'altann': 'c.+3909_+3916del8',
                'altclass': '.'
            }
        ]

        expected[('X', '3523053', 'GACACACACAC', 'G')] = [
            {
                'quality': 'low',
                'TR': 6,
                'TC': 35,
                'NF': 3,
                'NR': 3,
                'gene': 'PRKX',
                'csn': 'c.+3947_+3956del10',
                'class_': '3PU',
                'altann': 'c.+3909_+3918del10',
                'altclass': '.'
            }
        ]

        parsed_record = parsers.parse_vcf_record(vcf_record)
        self.assertDictEqual(parsed_record, expected)


    def test_parse_vcf_record_single_alt_multiple_transcript(self):

        vcf_record = 'X       152728155       .       T       C       200     PASS    ABPV=1.00e+00;FPV=2.03e-03;' \
                     'FR=1.0000;HP=2;MMLQ=32;NF=5;NR=5;PP=200;RMP=60.61;RPV=1.00e+00;SC=GAGCCCATCATGGCAGGCACC;TC=12;' \
                     'TCF=7;TCR=5;TR=10;TYPE=Substitution;TRANSCRIPT=ENST00000330912:ENST00000370219;GENE=TREX2:HAUS7;' \
                     'GENEID=ENSG00000183479:ENSG00000213397;TRINFO=-/25.9kb/13/2.3kb/236:-/23.5kb/10/1.9kb/368;' \
                     'LOC=In3/4:In3/4;CSN=c.-1229-26A>G:c.323-26A>G;PROTPOS=.:.;PROTREF=.:.;PROTALT=.:.;CLASS=5PU:INT;' \
                     'SO=5_prime_UTR_variant:intron_variant;IMPACT=3:3;ALTANN=.:.;ALTCLASS=.:.;ALTSO=.:.    GT:GL:GQ:NR' \
                     '     1/1:-89.91,-15.44,-8.51:57:12'

        expected = {}
        expected[('X', '152728155', 'T', 'C')] = [
            {
                'quality': 'high',
                'TR': 10,
                'TC': 12,
                'NF': 5,
                'NR': 5,
                'gene': 'TREX2',
                'csn': 'c.-1229-26A>G',
                'class_': '5PU',
                'altann': '.',
                'altclass': '.'
            },
            {
                'quality': 'high',
                'TR': 10,
                'TC': 12,
                'NF': 5,
                'NR': 5,
                'gene': 'HAUS7',
                'csn': 'c.323-26A>G',
                'class_': 'INT',
                'altann': '.',
                'altclass': '.'
            }
        ]

        parsed_record = parsers.parse_vcf_record(vcf_record)
        self.assertDictEqual(parsed_record, expected)


    def test_parse_vcf_record_multiple_alt_multiple_transcript(self):

        vcf_record = 'chrX       106184601       .       TGAGAGAGAGAGA   TGAGA,T 200     PASS    ABPV=4.76e-01;FPV=3.02e-07;' \
                     'FR=0.5000,0.5000;HP=4;NF=5,4;NR=5,2;PP=200,200;RMP=61.78;RPV=2.07e-02;SC=AAGGTGAGGGTGAGAGAGAGA;' \
                     'TC=39;TCF=19;TCR=20;TR=10,6;TYPE=Deletion,Deletion;TRANSCRIPT=ENST00000355610:ENST00000358740,' \
                     'ENST00000355610:ENST00000358740;GENE=MORC4:AL158821.1,MORC4:AL158821.1;GENEID=ENSG00000133131:' \
                     'ENSG00000196333,ENSG00000133131:ENSG00000196333;TRINFO=-/59.5kb/17/3.8kb/937:' \
                     '-/179.5kb/2/2.0kb/246,-/59.5kb/17/3.8kb/937:-/179.5kb/2/2.0kb/246;LOC=3UTR:In1/2,3UTR:In1/2;' \
                     'CSN=c.+100_+107del8:c.-139+51851_-139+51858del8,c.+96_+107del12:c.-139+51847_-139+51858del12;' \
                     'PROTPOS=.:.,.:.;PROTREF=.:.,.:.;PROTALT=.:.,.:.;CLASS=3PU:5PU,3PU:5PU;SO=3_prime_UTR_variant:' \
                     '5_prime_UTR_variant,3_prime_UTR_variant:5_prime_UTR_variant;IMPACT=3:3,3:3;ALTANN=c.+68_+75del8:' \
                     'c.-139+51819_-139+51826del8,c.+68_+79del12:c.-139+51819_-139+51830del12;ALTCLASS=.:.,.:.;ALTSO=.' \
                     ':.,.:. GT:GL:GQ:NR     1/2:-1,-1,-1:100:10922'

        expected = {}
        expected[('X', '106184601', 'TGAGAGAGAGAGA', 'TGAGA')] = [
            {
                'quality': 'high',
                'TR': 10,
                'TC': 39,
                'NF': 5,
                'NR': 5,
                'gene': 'MORC4',
                'csn': 'c.+100_+107del8',
                'class_': '3PU',
                'altann': 'c.+68_+75del8',
                'altclass': '.'
            },
            {
                'quality': 'high',
                'TR': 10,
                'TC': 39,
                'NF': 5,
                'NR': 5,
                'gene': 'AL158821.1',
                'csn': 'c.-139+51851_-139+51858del8',
                'class_': '5PU',
                'altann': 'c.-139+51819_-139+51826del8',
                'altclass': '.'
            }
        ]
        expected[('X', '106184601', 'TGAGAGAGAGAGA', 'T')] = [
            {
                'quality': 'low',
                'TR': 6,
                'TC': 39,
                'NF': 4,
                'NR': 2,
                'gene': 'MORC4',
                'csn': 'c.+96_+107del12',
                'class_': '3PU',
                'altann': 'c.+68_+79del12',
                'altclass': '.'
            },
            {
                'quality': 'low',
                'TR': 6,
                'TC': 39,
                'NF': 4,
                'NR': 2,
                'gene': 'AL158821.1',
                'csn': 'c.-139+51847_-139+51858del12',
                'class_': '5PU',
                'altann': 'c.-139+51819_-139+51830del12',
                'altclass': '.'
            }
        ]

        parsed_record = parsers.parse_vcf_record(vcf_record)
        self.assertDictEqual(parsed_record, expected)


    def test_parse_vcf_record_tc_is_zero(self):

        vcf_record = '11      1642299 .       T       A       200     PASS    ABPV=1.00e+00;FPV=5.41e-02;FR=0.5000;' \
                     'HP=1;MMLQ=34;NF=14;NR=13;PP=200;RMP=57.26;RPV=9.05e-01;SC=CCCACCTTGTTGCAGGTGGGA;TC=0;TCF=39;' \
                     'TCR=21;TR=27;TYPE=Substitution;TRANSCRIPT=ENST00000399682;GENE=KRTAP5-4;GENEID=ENSG00000241598;' \
                     'TRINFO=-/1.2kb/1/1.2kb/228;LOC=3UTR;CSN=c.+338A>T;PROTPOS=.;PROTREF=.;PROTALT=.;CLASS=3PU;' \
                     'SO=3_prime_UTR_variant;IMPACT=3;ALTANN=.;ALTCLASS=.;ALTSO=.'

        self.assertIsNone(parsers.parse_vcf_record(vcf_record))

