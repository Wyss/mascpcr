
import json
import os
import unittest

import bitarray as bt

from mascpcr import genbankfeatures
from ._common import RECODED_GB, TEST_OUTPUT_DIR

class TestGenbankfeatures(unittest.TestCase):

    def test_filterFeatures(self):
        feat_list = genbankfeatures.filterFeatures(
            RECODED_GB, 
            ['synth_fragment'], 
            {'label': r'seg23.*'})
        feat_list = [str(ft) for ft in feat_list]
        check_fp = os.path.join(TEST_OUTPUT_DIR, 'test_filterfeatures.out')
        with open(check_fp) as fd:
            check_list = json.loads(fd.read())
        self.assertEqual(feat_list, check_list)

    def test_buildBorderLUT(self):
        # Strand agnostic
        border_lut = genbankfeatures.buildBorderLUT(
            RECODED_GB,
            ['synth_fragment'], 
            {'label': r'seg23.*'})
        check_fp = os.path.join(TEST_OUTPUT_DIR, 'test_buildborderlut.out')
        with open(check_fp) as fd:
            check_lut = bt.bitarray()
            check_lut.fromfile(fd)
        self.assertEqual(border_lut, check_lut)
        # Strand specific
        fwd_border_lut, rev_border_lut = genbankfeatures.buildBorderLUT(
            RECODED_GB,
            ['synth_fragment'], 
            {'label': r'seg23.*'},
            strand_specific=True)
        check_fp_fwd = os.path.join(TEST_OUTPUT_DIR, 
                                    'test_buildborderlut_fwd.out')
        check_fp_rev= os.path.join(TEST_OUTPUT_DIR, 
                                    'test_buildborderlut_rev.out')
        with open(check_fp_fwd) as fd:
            check_lut_fwd = bt.bitarray()
            check_lut_fwd.fromfile(fd)
        self.assertEqual(fwd_border_lut, check_lut_fwd)  
        with open(check_fp_rev) as fd:
            check_lut_rev = bt.bitarray()
            check_lut_rev.fromfile(fd)
        self.assertEqual(rev_border_lut, check_lut_rev)        

    def test_findAggregateBoundaries(self):
        agg_bounds = genbankfeatures.findAggregateBoundaries(
            RECODED_GB,
            ['synth_fragment'], 
            {'label': r'seg23.*'})
        self.assertEqual(agg_bounds, (1039203, 1084477))
