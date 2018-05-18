"""Unit tests for the filters module"""

from unittest import TestCase
import pytest
from main import filters



class TestFilters(TestCase):


    def setUp(self):

        # SELF.CONFIG TO BE DOWNLOADED LATER
        self.config = {
            'REMOVE_MULTI_ALLELE_CALLS': True,
            'SPLICE_SITE_BOUNDARY': 10,
            'CHILD_MIN_TR': 3,
            'CHILD_MIN_TC': 15,
            'CHILD_MIN_TR_PER_TC': 0.2,
            'CONTROL_MAX_FREQUENCY': 0.1,
            'GNOMAD_MAX_FREQUENCY': 0.1,
            'PARENT_MIN_COVERAGE': 6,
            'PARENT_MAX_ALT_ALLELE_COUNT': 1,
        }

        self.filt = filters.Filters(True)


    def _test_check_if_multiallelic(self):

        multi_calls = [('7', 41231134, 'GGG'), ('12', 12345678, 'AC')]

        assert self.filt.check_if_multiallelic(
            {'REMOVE_MULTI_ALLELE_CALLS': True},
            {'multiallelic_calls': multi_calls},
            ('12', 12345678, 'AC', 'T')
        ) is True

        assert self.filt.check_if_multiallelic(
            {'REMOVE_MULTI_ALLELE_CALLS': True},
            {'multiallelic_calls': multi_calls},
            ('2', 12345678, 'AC', 'T')
        ) is False

        assert self.filt.check_if_multiallelic(
            {'REMOVE_MULTI_ALLELE_CALLS': True},
            {'multiallelic_calls': []},
            ('12', 12345678, 'AC', 'T')
        ) is False


    def _test_check_if_called_in_parent(self):

        data = {
            'mother_var': {
                ('12', 12345678, 'AC', 'T'): {},
                ('7', 41231134, 'TTT', 'A'): {},
            },
            'father_var': {
                ('5', 4779121, 'G', 'C'): {}
            }
        }

        assert self.filt.check_if_called_in_parent(
            data,
            ('12', 12345678, 'AC', 'T')
        ) is True

        assert self.filt.check_if_called_in_parent(
            data,
            ('5', 4779121, 'G', 'C')
        ) is True

        assert self.filt.check_if_called_in_parent(
            data,
            ('9', 24681012, 'A', 'T')
        ) is False


    def _test_check_if_low_quality(self):

        assert self.filt.check_if_low_quality({'quality': 'low'}) is True
        assert self.filt.check_if_low_quality({'quality': 'else'}) is False


    def _test_check_if_outside_splice_side_boundary(self):

        assert self.filt.check_if_outside_splice_side_boundary(
            {'SPLICE_SITE_BOUNDARY': 10},
            {'csn': '....'}
        ) is True

        # More asserts