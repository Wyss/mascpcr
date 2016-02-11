
import pickle
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


        EXPECTED_OUTPUT_FP = os.path.join(
            TEST_OUTPUT_DIR, 'test_filterfeatures.out')

        # NOTE: Used to create expected output. If code breaks, you should
        # figure out whether there is really a bug before uncommenting the
        # following and changing the expected output.
        # with open(EXPECTED_OUTPUT_FP, 'w') as fd:
        #     pickle.dump(feat_list, fd)

        # Load expected output.
        with open(EXPECTED_OUTPUT_FP) as fd:
            check_feature_list = pickle.loads(fd.read())

        # Check all features preserved.
        for i in range(len(feat_list)):
            self.assertEqual(str(check_feature_list[i]), str(feat_list[i]))

    def test_buildBorderLUT(self):
        # Strand agnostic
        border_lut = genbankfeatures.buildBorderLUT(
            RECODED_GB,
            ['synth_fragment'], 
            {'label': r'seg23.*'})
        check_fp = os.path.join(TEST_OUTPUT_DIR, 'test_buildborderlut.out')
        # with open(check_fp, 'wb') as fd:
        #     border_lut.tofile(fd)
        with open(check_fp, 'rb') as fd:
            check_lut = bt.bitarray()
            check_lut.fromfile(fd)
        self.assertEqual(border_lut, check_lut[:len(border_lut)])
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
        # with open(check_fp_fwd, 'wb') as fd:
        #     fwd_border_lut.tofile(fd)
        # with open(check_fp_rev, 'wb') as fd:
        #     rev_border_lut.tofile(fd)
        with open(check_fp_fwd, 'rb') as fd:
            check_lut_fwd = bt.bitarray()
            check_lut_fwd.fromfile(fd)
        self.assertEqual(fwd_border_lut, check_lut_fwd[:len(fwd_border_lut)])  
        with open(check_fp_rev, 'rb') as fd:
            check_lut_rev = bt.bitarray()
            check_lut_rev.fromfile(fd)
        self.assertEqual(rev_border_lut, check_lut_rev[:len(rev_border_lut)])        

    def test_findAggregateBoundaries(self):
        agg_bounds = genbankfeatures.findAggregateBoundaries(
            RECODED_GB,
            ['synth_fragment'], 
            {'label': r'seg23.*'})
        self.assertEqual(agg_bounds, (1127000, 1176000))
