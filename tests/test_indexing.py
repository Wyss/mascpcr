
import os
import unittest

import numpy as np

from mascpcr import indexing
from ._common import RECODED_GENOME_FP, REFERENCE_GENOME_FP, TEST_OUTPUT_DIR, \
                     TEST_CACHE_DIR, REFERENCE_GB_STR, RECODED_GB_STR


class TestIndexing(unittest.TestCase):

    def test_buildIdxLUT(self):
        idx_lut = indexing.buildIdxLUT(
            RECODED_GENOME_FP,
            REFERENCE_GENOME_FP,
            RECODED_GB_STR,
            REFERENCE_GB_STR,
            cache_dir=TEST_CACHE_DIR)
        check_idx_lut = np.save(os.path.join(TEST_OUTPUT_DIR,
                                             'test_buildidxlut2.out'), idx_lut)
        check_idx_lut = np.load(os.path.join(TEST_OUTPUT_DIR,
                                             'test_buildidxlut.out'))
        self.assertTrue((idx_lut == check_idx_lut).all())

    def test_buildEdgeLUT(self):
        idx_lut = indexing.buildIdxLUT(
            RECODED_GENOME_FP,
            REFERENCE_GENOME_FP,
            RECODED_GB_STR,
            REFERENCE_GB_STR,
            cache_dir=TEST_CACHE_DIR)
        edge_lut = indexing.buildEdgeLUT(
            idx_lut,
            TEST_CACHE_DIR)
        with open(os.path.join(TEST_OUTPUT_DIR,
                              'test_buildedgelut.out'), 'wb') as fd:
            edge_lut.tofile(fd)
        with open(os.path.join(TEST_OUTPUT_DIR,
                              'test_buildedgelut.out'), 'rb') as fd:
            check_edge_lut = np.fromfile(fd, dtype=np.bool)
        self.assertTrue((edge_lut == check_edge_lut).all())


    def test_buildMismatchLUT(self):
        idx_lut = indexing.buildIdxLUT(
            RECODED_GENOME_FP,
            REFERENCE_GENOME_FP,
            RECODED_GB_STR,
            REFERENCE_GB_STR,
            cache_dir=TEST_CACHE_DIR)
        mismatch_lut = indexing.buildMismatchLUT(
            idx_lut,
            RECODED_GB_STR,
            REFERENCE_GB_STR,
            TEST_CACHE_DIR)
        with open(os.path.join(TEST_OUTPUT_DIR,
                              'test_buildmismatchlut.out'), 'wb') as fd:
            mismatch_lut.tofile(fd)
        with open(os.path.join(TEST_OUTPUT_DIR,
                              'test_buildmismatchlut.out'), 'rb') as fd:
            check_mismatch_lut = np.fromfile(fd, dtype=np.bool)
        self.assertTrue((mismatch_lut == check_mismatch_lut).all())

