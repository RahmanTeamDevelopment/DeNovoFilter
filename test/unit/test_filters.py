"""Unit tests for the filters module"""

from unittest import TestCase
from main import filters



class TestFilters(TestCase):


    def setUp(self):

        self.filt = filters.Filters(False)


    def assert_fails_check(self, errmsg='', fun=lambda: None, argv=None):

        with self.assertRaises(ValueError) as cm:
            fun(*argv)
        err = cm.exception
        self.assertEqual(str(err), errmsg)


    def assert_passes_check(self, fun=lambda: None, argv=None):

        try:
            fun(*argv)
        except ValueError:
            self.fail("{}() raised ValueError unexpectedly".format(fun.func_name))


    def test_check_if_multiallelic(self):

        multi_calls = [('7', 41231134, 'GGG'), ('12', 12345678, 'AC')]

        self.assert_fails_check(
            errmsg='multi_allele_call',
            fun=self.filt.check_if_multiallelic,
            argv=[
                {'REMOVE_MULTI_ALLELE_CALLS': True},
                {'multiallelic_calls': multi_calls},
                ('12', 12345678, 'AC', 'T')
            ]
        )

        self.assert_passes_check(
            fun=self.filt.check_if_multiallelic,
            argv=[
                {'REMOVE_MULTI_ALLELE_CALLS': True},
                {'multiallelic_calls': multi_calls},
                ('2', 12345678, 'AC', 'T')
            ]
        )

        self.assert_passes_check(
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

        self.assert_fails_check(
            errmsg='called_in_parent',
            fun=self.filt.check_if_called_in_parent,
            argv=[
                data,
                ('12', 12345678, 'AC', 'T')
            ]
        )

        self.assert_fails_check(
            errmsg='called_in_parent',
            fun=self.filt.check_if_called_in_parent,
            argv=[
                data,
                ('5', 4779121, 'G', 'C')
            ]
        )

        self.assert_passes_check(
            fun=self.filt.check_if_called_in_parent,
            argv=[
                data,
                ('9', 24681012, 'A', 'T')
            ]
        )


    def test_check_if_low_quality(self):

        self.assert_fails_check(
            errmsg='low_quality',
            fun=self.filt.check_if_low_quality,
            argv=[{'quality': 'low'}]
        )

        self.assert_passes_check(
            fun=self.filt.check_if_low_quality,
            argv=[{'quality': 'else'}]
        )


    def test_check_if_outside_splice_side_boundary(self):

        self.assert_passes_check(
            fun=self.filt.check_if_outside_splice_side_boundary,
            argv=[
                {'SPLICE_SITE_BOUNDARY': 10},
                {'csn': 'c.111A>C_p.='}
            ]
        )

        # TODO: Add more asserts here!
        self.fail('Need to add more asserts here')


    def test_check_tr_in_child(self):

        self.assert_fails_check(
            errmsg='low_child_tr ({})'.format(1),
            fun=self.filt.check_tr_in_child,
            argv=[
                {'CHILD_MIN_TR': 3},
                {'TR': 1}
            ]
        )

        self.assert_passes_check(
            fun=self.filt.check_tr_in_child,
            argv=[
                {'CHILD_MIN_TR': 3},
                {'TR': 4}
            ]
        )

        self.assert_passes_check(
            fun=self.filt.check_tr_in_child,
            argv=[
                {'CHILD_MIN_TR': 3},
                {'TR': 3}
            ]
        )


    def test_check_tc_in_child(self):

        self.assert_fails_check(
            errmsg='low_child_tc ({})'.format(11),
            fun=self.filt.check_tc_in_child,
            argv=[
                {'CHILD_MIN_TC': 15},
                {'TC': 11}
            ]
        )

        self.assert_passes_check(
            fun=self.filt.check_tc_in_child,
            argv=[
                {'CHILD_MIN_TC': 15},
                {'TC': 18}
            ]
        )

        self.assert_passes_check(
            fun=self.filt.check_tc_in_child,
            argv=[
                {'CHILD_MIN_TC': 15},
                {'TC': 15}
            ]
        )


    def test_check_tr_per_tc_in_child(self):

        self.assert_fails_check(
            errmsg='low_child_tr_per_tc ({})'.format(0.1),
            fun=self.filt.check_tr_per_tc_in_child,
            argv=[
                {'CHILD_MIN_TR_PER_TC': 0.2},
                {'TR': 1, 'TC': 10}
            ]
        )

        self.assert_passes_check(
            fun=self.filt.check_tr_per_tc_in_child,
            argv=[
                {'CHILD_MIN_TR_PER_TC': 0.2},
                {'TR': 3, 'TC': 10}
            ]
        )

        self.assert_passes_check(
            fun=self.filt.check_tr_per_tc_in_child,
            argv=[
                {'CHILD_MIN_TR_PER_TC': 0.2},
                {'TR': 2, 'TC': 10}
            ]
        )


    def test_check_control_frequency(self):

        self.assert_fails_check(
            errmsg='high_control_frequency ({})'.format(0.15),
            fun=self.filt.check_control_frequency,
            argv=[
                {'CONTROL_MAX_FREQUENCY': 0.1},
                0.15
            ]
        )

        self.assert_passes_check(
            fun=self.filt.check_control_frequency,
            argv=[
                {'CONTROL_MAX_FREQUENCY': 0.1},
                0.09
            ]
        )

        self.assert_passes_check(
            fun=self.filt.check_control_frequency,
            argv=[
                {'CONTROL_MAX_FREQUENCY': 0.1},
                0.1
            ]
        )


    def test_check_gnomad_exomes_frequency(self):

        self.assert_fails_check(
            errmsg='high_gnomad_exomes_frequency ({})'.format(0.15),
            fun=self.filt.check_gnomad_exomes_frequency,
            argv=[
                {'GNOMAD_MAX_FREQUENCY': 0.1},
                0.15
            ]
        )

        self.assert_passes_check(
            fun=self.filt.check_gnomad_exomes_frequency,
            argv=[
                {'GNOMAD_MAX_FREQUENCY': 0.1},
                0.09
            ]
        )

        self.assert_passes_check(
            fun=self.filt.check_gnomad_exomes_frequency,
            argv=[
                {'GNOMAD_MAX_FREQUENCY': 0.1},
                0.1
            ]
        )


    def test_check_gnomad_genomes_frequency(self):

        self.assert_fails_check(
            errmsg='high_gnomad_genomes_frequency ({})'.format(0.15),
            fun=self.filt.check_gnomad_genomes_frequency,
            argv=[
                {'GNOMAD_MAX_FREQUENCY': 0.1},
                0.15
            ]
        )

        self.assert_passes_check(
            fun=self.filt.check_gnomad_genomes_frequency,
            argv=[
                {'GNOMAD_MAX_FREQUENCY': 0.1},
                0.09
            ]
        )

        self.assert_passes_check(
            fun=self.filt.check_gnomad_genomes_frequency,
            argv=[
                {'GNOMAD_MAX_FREQUENCY': 0.1},
                0.1
            ]
        )


    def test_check_tc_and_tr_in_mother(self):

        self.assert_fails_check(
            errmsg='low_mother_tc ({})'.format(4),
            fun=self.filt.check_tc_and_tr_in_mother,
            argv=[
                {'PARENT_MIN_COVERAGE': 6, 'PARENT_MAX_ALT_ALLELE_COUNT': 1},
                {'mother_tc': 4, 'mother_tr': 0}
            ]
        )

        self.assert_passes_check(
            fun=self.filt.check_tc_and_tr_in_mother,
            argv=[
                {'PARENT_MIN_COVERAGE': 6, 'PARENT_MAX_ALT_ALLELE_COUNT': 1},
                {'mother_tc': 7, 'mother_tr': 0}
            ]
        )

        self.assert_passes_check(
            fun=self.filt.check_tc_and_tr_in_mother,
            argv=[
                {'PARENT_MIN_COVERAGE': 6, 'PARENT_MAX_ALT_ALLELE_COUNT': 1},
                {'mother_tc': 6, 'mother_tr': 0}
            ]
        )

        self.assert_fails_check(
            errmsg='high_mother_tr ({})'.format(2),
            fun=self.filt.check_tc_and_tr_in_mother,
            argv=[
                {'PARENT_MIN_COVERAGE': 6, 'PARENT_MAX_ALT_ALLELE_COUNT': 1},
                {'mother_tc': 7, 'mother_tr': 2}
            ]
        )

        self.assert_passes_check(
            fun=self.filt.check_tc_and_tr_in_mother,
            argv=[
                {'PARENT_MIN_COVERAGE': 6, 'PARENT_MAX_ALT_ALLELE_COUNT': 1},
                {'mother_tc': 7, 'mother_tr': 0}
            ]
        )

        self.assert_passes_check(
            fun=self.filt.check_tc_and_tr_in_mother,
            argv=[
                {'PARENT_MIN_COVERAGE': 6, 'PARENT_MAX_ALT_ALLELE_COUNT': 1},
                {'mother_tc': 7, 'mother_tr': 1}
            ]
        )


    def test_check_tc_and_tr_in_father(self):

        self.assert_fails_check(
            errmsg='low_father_tc ({})'.format(4),
            fun=self.filt.check_tc_and_tr_in_father,
            argv=[
                {'PARENT_MIN_COVERAGE': 6, 'PARENT_MAX_ALT_ALLELE_COUNT': 1},
                {'father_tc': 4, 'father_tr': 0}
            ]
        )

        self.assert_passes_check(
            fun=self.filt.check_tc_and_tr_in_father,
            argv=[
                {'PARENT_MIN_COVERAGE': 6, 'PARENT_MAX_ALT_ALLELE_COUNT': 1},
                {'father_tc': 7, 'father_tr': 0}
            ]
        )

        self.assert_passes_check(
            fun=self.filt.check_tc_and_tr_in_father,
            argv=[
                {'PARENT_MIN_COVERAGE': 6, 'PARENT_MAX_ALT_ALLELE_COUNT': 1},
                {'father_tc': 6, 'father_tr': 0}
            ]
        )

        self.assert_fails_check(
            errmsg='high_father_tr ({})'.format(2),
            fun=self.filt.check_tc_and_tr_in_father,
            argv=[
                {'PARENT_MIN_COVERAGE': 6, 'PARENT_MAX_ALT_ALLELE_COUNT': 1},
                {'father_tc': 7, 'father_tr': 2}
            ]
        )

        self.assert_passes_check(
            fun=self.filt.check_tc_and_tr_in_father,
            argv=[
                {'PARENT_MIN_COVERAGE': 6, 'PARENT_MAX_ALT_ALLELE_COUNT': 1},
                {'father_tc': 7, 'father_tr': 0}
            ]
        )

        self.assert_passes_check(
            fun=self.filt.check_tc_and_tr_in_father,
            argv=[
                {'PARENT_MIN_COVERAGE': 6, 'PARENT_MAX_ALT_ALLELE_COUNT': 1},
                {'father_tc': 7, 'father_tr': 1}
            ]
        )




    def test_apply(self):

        self.filt.detailed = False
        self.filt.filter = []

        with self.assertRaises(ValueError) as cm:
            self.filt._apply(True, 'some message')
        err = cm.exception
        self.assertEqual(str(err), 'some message')
        self.assertEquals(len(self.filt.filter), 0)

        try:
            self.filt._apply(False, 'some message')
        except ValueError:
            self.fail("_apply() raised ValueError unexpectedly")
        self.assertEquals(len(self.filt.filter), 0)

        self.filt.detailed = True

        self.filt.filter = []
        try:
            self.filt._apply(True, 'some message')
        except ValueError:
            self.fail("_apply() raised ValueError unexpectedly")
        self.assertEquals(self.filt.filter[0], 'some message')

        self.filt.filter = []
        try:
            self.filt._apply(False, 'some message')
        except ValueError:
            self.fail("_apply() raised ValueError unexpectedly")
        self.assertEquals(len(self.filt.filter), 0)


    def test_apply_filters(self):

        self.fail('Test not yet implemented')