"""Unit tests for the filters module"""

from unittest import TestCase
from main import filters



class TestFilters(TestCase):


    def setUp(self):

        self.filt = filters.Filters(False)


    def assert_exception(self, errmsg='', fun=lambda: None, argv=None):

        with self.assertRaises(ValueError) as cm:
            fun(*argv)
        err = cm.exception
        self.assertEqual(str(err), errmsg)


    def assert_no_exception(self, fun=lambda: None, argv=None):

        try:
            fun(*argv)
        except ValueError:
            self.fail("{}() raised ValueError unexpectedly".format(fun.func_name))




    def test_check_if_multiallelic(self):

        multi_calls = [('7', 41231134, 'GGG'), ('12', 12345678, 'AC')]

        self.assert_exception(
            errmsg='multi_allele_call',
            fun=self.filt.check_if_multiallelic,
            argv=[
                {'REMOVE_MULTI_ALLELE_CALLS': True},
                {'multiallelic_calls': multi_calls},
                ('12', 12345678, 'AC', 'T')
            ]
        )

        self.assert_no_exception(
            fun=self.filt.check_if_multiallelic,
            argv=[
                {'REMOVE_MULTI_ALLELE_CALLS': True},
                {'multiallelic_calls': multi_calls},
                ('2', 12345678, 'AC', 'T')
            ]
        )

        self.assert_no_exception(
            fun=self.filt.check_if_multiallelic,
            argv=[
                {'REMOVE_MULTI_ALLELE_CALLS': True},
                {'multiallelic_calls': []},
                ('12', 12345678, 'AC', 'T')
            ]
        )


    def test_check_if_called_in_parent(self):

        data = {
            'mother_var': {
                ('12', 12345678, 'AC', 'T'): {},
                ('7', 41231134, 'TTT', 'A'): {},
            },
            'father_var': {
                ('5', 4779121, 'G', 'C'): {}
            }
        }

        self.assert_exception(
            errmsg='called_in_parent',
            fun=self.filt.check_if_called_in_parent,
            argv=[
                data,
                ('12', 12345678, 'AC', 'T')
            ]
        )

        self.assert_exception(
            errmsg='called_in_parent',
            fun=self.filt.check_if_called_in_parent,
            argv=[
                data,
                ('5', 4779121, 'G', 'C')
            ]
        )

        self.assert_no_exception(
            fun=self.filt.check_if_called_in_parent,
            argv=[
                data,
                ('9', 24681012, 'A', 'T')
            ]
        )


    def test_check_if_low_quality(self):

        self.assert_exception(
            errmsg='low_quality',
            fun=self.filt.check_if_low_quality,
            argv=[{'quality': 'low'}]
        )

        self.assert_no_exception(
            fun=self.filt.check_if_low_quality,
            argv=[{'quality': 'else'}]
        )


    def test_check_if_outside_splice_side_boundary(self):

        # More asserts

        self.assert_no_exception(
            fun=self.filt.check_if_outside_splice_side_boundary,
            argv=[
                {'SPLICE_SITE_BOUNDARY': 10},
                {'csn': 'c.111A>C_p.='}
            ]
        )

