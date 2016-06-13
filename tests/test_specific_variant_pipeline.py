import os
import tempfile
import unittest

import pandas as pd

from ._common import TEST_INPUT_DIR
from mascpcr.specific_variant_pipeline import _create_modified_genome_str
from mascpcr.specific_variant_pipeline import find_masc_pcr_primers_given_variants
from mascpcr.specific_variant_pipeline import OUTPUT_COLUMNS_ORDERED


TEST_MASC_PCR_INPUT = os.path.join(
        TEST_INPUT_DIR, 'fix_recoli_exp5_variants.csv')

TEST_REF_GENOME = os.path.join(TEST_INPUT_DIR, 'mg1655.NC_000913.2.fa')


class TestSpecificVariantPipeline(unittest.TestCase):

    def test_find_masc_pcr_primers_given_variants(self):
        """Basic test for finding primers given variants.
        """
        with tempfile.NamedTemporaryFile() as output_file_csv_fh:
            output_file_csv = output_file_csv_fh.name
        find_masc_pcr_primers_given_variants(
                TEST_MASC_PCR_INPUT, TEST_REF_GENOME, output_file_csv)

        # Check the results.
        results_df = pd.read_csv(output_file_csv)
        input_df = pd.read_csv(TEST_MASC_PCR_INPUT)

        # Expected columns present.
        self.assertEqual(0,
                len(set(OUTPUT_COLUMNS_ORDERED) - set(results_df.columns)))

        # All variants accounted for.
        self.assertEqual(set(input_df.POSITION), set(results_df.POSITION))

    def test_create_modified_genome_str__single_nucleotide_mutation(self):
        ref_genome_str = 'AAAAAA'
        variant_list = [{
            'POSITION': 3,
            'REF': 'A',
            'ALT': 'G'
        }]
        variants_df = pd.DataFrame(variant_list)
        modified_genome_str = _create_modified_genome_str(
                ref_genome_str, variants_df)
        self.assertEqual('AAGAAA', modified_genome_str)

    def test_create_modified_genome_str__insertion(self):
        ref_genome_str = 'AAAAAA'
        variant_list = [{
            'POSITION': 3,
            'REF': 'A',
            'ALT': 'AG'
        }]
        variants_df = pd.DataFrame(variant_list)
        modified_genome_str = _create_modified_genome_str(
                ref_genome_str, variants_df)
        self.assertEqual('AAAGAAA', modified_genome_str)

    def test_create_modified_genome_str__multiple(self):
        ref_genome_str = 'AATAAATAAA'
        variant_list = [
            {
                'POSITION': 3,
                'REF': 'T',
                'ALT': 'G'
            },
            {
                'POSITION': 7,
                'REF': 'T',
                'ALT': 'TG'
            }
        ]
        variants_df = pd.DataFrame(variant_list)
        modified_genome_str = _create_modified_genome_str(
                ref_genome_str, variants_df)
        self.assertEqual('AAGAAATGAAA', modified_genome_str)
