"""Unit tests for the helper module"""

from unittest import TestCase
from main import helper



class TestHelper(TestCase):




    def test_find_multiallelic_calls(self):

        var_data = {

            (None, None, None, None): None,
            (None, None, None, None): None,
            (None, None, None, None): None,
            (None, None, None, None): None,
            (None, None, None, None): None,

        }
