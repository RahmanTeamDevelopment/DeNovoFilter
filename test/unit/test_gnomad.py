from unittest import TestCase
from main import gnomad
from mock import patch


class TestGnomadDBReader(TestCase):


    @patch('main.gnomad.GnomadDBReader._extract_pops')
    @patch('main.gnomad.GnomadDBReader._read_variants_in_vicinity')
    @patch('main.gnomad.GnomadDBReader.__init__')
    def test_get_max_frequency_autosome(self, mocked_init, mocked_read_variants_in_vicinity, mocked_extract_pops):

        mocked_init.return_value = None
        mocked_extract_pops.return_value = ['AFR', 'AMR', 'ASJ', 'EAS', 'FIN', 'NFE', 'OTH', 'SAS']

        flags = {
            'pos': 20699518,
            'ref': 'C',
            'alts': ['T'],
            'GENE': 'NELL1',
            'CSN': 'c.96C>T_p.=',
            'GC_AFR': '7650,0,0',
            'GC_AMR': '16786,1,0',
            'GC_ASJ': '4924,0,0',
            'GC_EAS': '8624,0,0',
            'GC_FIN': '11149,0,0',
            'GC_NFE': '55800,3,0',
            'GC_OTH': '2742,0,0',
            'GC_SAS': '15376,15,0',
        }
        mocked_read_variants_in_vicinity.return_value = [(1, flags)]

        reader = gnomad.GnomadDBReader()

        result = reader.get_max_frequency(
            ('11', '20699518', 'C', 'T'),
            ('NELL1', 'c.96C>T_p.=')
        )
        self.assertEquals(round(result[0], 11), 0.09745955428)
        self.assertEquals(result[1], 'SAS')


    @patch('main.gnomad.GnomadDBReader._extract_pops')
    @patch('main.gnomad.GnomadDBReader._read_variants_in_vicinity')
    @patch('main.gnomad.GnomadDBReader.__init__')
    def test_get_max_frequency_chrX(self, mocked_init, mocked_read_variants_in_vicinity, mocked_extract_pops):

        mocked_init.return_value = None
        mocked_extract_pops.return_value = ['AFR', 'AMR', 'ASJ', 'EAS', 'FIN', 'NFE', 'OTH', 'SAS']

        flags = {
            'pos': 66765242,
            'ref': 'GGCA',
            'alts': ['G', 'GGCAGCA'],
            'GENE': 'AR',
            'CSN': 'c.271_273delCAG_p.Gln91del',
            'GC_AFR_Male': '1836,0,0,0,0,0',
            'GC_AMR_Male': '5064,0,0,0,0,0',
            'GC_ASJ_Male': '2089,0,0,0,0,0',
            'GC_EAS_Male': '2891,0,0,0,0,0',
            'GC_FIN_Male': '3727,0,0,0,0,0',
            'GC_NFE_Male': '17980,0,0,0,0,1',
            'GC_OTH_Male': '1031,0,0,0,0,0',
            'GC_SAS_Male': '9069,0,2,0,0,1',
            'GC_AFR_Female': '3710,0,0,2,0,0',
            'GC_AMR_Female': '8517,0,0,0,0,0',
            'GC_ASJ_Female': '2132,0,0,0,0,0',
            'GC_EAS_Female': '3549,0,0,0,0,0',
            'GC_FIN_Female': '4570,0,0,0,0,0',
            'GC_NFE_Female': '19787,2,0,3,0,0',
            'GC_OTH_Female': '1145,0,0,0,0,0',
            'GC_SAS_Female': '3418,0,0,1,0,0'
        }
        mocked_read_variants_in_vicinity.return_value = [(2, flags)]

        reader = gnomad.GnomadDBReader()

        result = reader.get_max_frequency(
            ('X', '66765242', 'GGCA', 'G'),
            ('AR', 'c.271_273delCAG_p.Gln91del')
        )
        self.assertEquals(round(result[0], 10), 0.0160115283)
        self.assertEquals(result[1], 'SAS')


    @patch('main.gnomad.GnomadDBReader._extract_pops')
    @patch('main.gnomad.GnomadDBReader._read_variants_in_vicinity')
    @patch('main.gnomad.GnomadDBReader.__init__')
    def test_get_max_frequency_chrY(self, mocked_init, mocked_read_variants_in_vicinity, mocked_extract_pops):

        mocked_init.return_value = None
        mocked_extract_pops.return_value = ['AFR', 'AMR', 'ASJ', 'EAS', 'FIN', 'NFE', 'OTH', 'SAS']

        flags = {
            'pos': 27190050,
            'ref': 'A',
            'alts': ['C'],
            'GENE': 'BPY2C',
            'CSN': 'c.78+43T>G',
            'AF_AFR': '1.45985e-01',
            'AF_AMR': '3.16706e-03',
            'AF_ASJ': '0.00000e+00',
            'AF_EAS': '0.00000e+00',
            'AF_FIN': '0.00000e+00',
            'AF_NFE': '0.00000e+00',
            'AF_OTH': '3.86100e-03',
            'AF_SAS': '0.00000e+00'
        }
        mocked_read_variants_in_vicinity.return_value = [(1, flags)]

        reader = gnomad.GnomadDBReader()

        result = reader.get_max_frequency(
            ('Y', '27190050', 'A', 'C'),
            ('BPY2C', 'c.78+43T>G')
        )
        self.assertEquals(result[0], 1.45985e-01)
        self.assertEquals(result[1], 'AFR')

