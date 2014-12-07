import filecmp
import os
import unittest

from mascpcr import pipeline, genbankfeatures
from ._common import RECODED_GENOME_FP, REFERENCE_GENOME_FP, TEST_OUTPUT_DIR, \
                     TEST_CACHE_DIR, REFERENCE_GB_STR, RECODED_GB_STR, \
                     REFERENCE_GB, RECODED_GB


class TestPipeline(unittest.TestCase):

    def test_findMascPrimers(self):

        start_idx, end_idx = genbankfeatures.findAggregateBoundaries(
            # `sr_obj` expects a SeqRecord object
            sr_obj=RECODED_GB,
            # These are the features types that we want to look for
            feature_types=['synth_fragment'],
            # We will look up the qualifier 'label' in each Feature object and check
            # to make sure that the regular expression "seg23.*" matches its contents
            # (regex refresher: seg23.* will match "seg23" followed by any characters,
            #  e.g., seg23_001)
            qualifier_regexs={'label':'seg23.*'}
        )

        genome_str, ref_genome_str, idx_lut, edge_lut, mismatch_lut, \
            border_lut = pipeline.generateLUTs(
                genome_fp=RECODED_GENOME_FP,
                ref_genome_fp=REFERENCE_GENOME_FP,
                start_idx=start_idx,       
                end_idx=end_idx,          
                border_feature_types=['synth_fragment'],
                cache_luts=True,
                cache_dir=TEST_CACHE_DIR
        )
        # We have to prevent the output file from containing the parameters
        # as it will dump the absolute filepaths, which makes file comparisons
        # more difficult
        params = {
            'dump_params': False,
            'output_fp': TEST_CACHE_DIR,
            'output_basename': 'seg23'  
        }
        pipeline.findMascPrimers(
            idx_lut=idx_lut,
            genome_str=RECODED_GB_STR,
            ref_genome_str=REFERENCE_GB_STR,
            start_idx=start_idx,
            end_idx=end_idx,
            edge_lut=edge_lut,
            mismatch_lut=mismatch_lut,
            border_lut=border_lut,
            params=params
        )
        # Now compare the output files to the expected output files
        output_report_fp = os.path.join(TEST_CACHE_DIR, 
                                        'seg23_masc_report.csv')
        check_output_report_fp = os.path.join(TEST_OUTPUT_DIR, 
                                              'seg23_masc_report.csv')
        self.assertTrue(filecmp.cmp(output_report_fp, check_output_report_fp))
