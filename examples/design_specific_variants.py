"""
design_specific_variants | examples/design_specific_variants.py
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This is an example script that demonstrates designing MASC-PCR primers for
specific variants.

The supporting code is still work in progress and further testing of the
specific_variant_pipeline may reveal additional bugs.
"""

import os
import sys

PACKAGE_DIR = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
MODULE_DIR = os.path.join(PACKAGE_DIR, 'mascpcr')

# This fanciness is only necessary to run the script in place without
# installing the pipeline. If you install the pipeline you can just import
# "mascpcr" and call it a day.

try:
    import mascpcr
except:
    if not 'mascpcr' in os.listdir(PACKAGE_DIR):
        raise IOError('`mascpcr` must be installed in your PYTHONPATH or '
                      'script must be run from within the package directory')
    sys.path.append(MODULE_DIR)
    sys.path.append(PACKAGE_DIR)
    import mascpcr

from mascpcr.specific_variant_pipeline import find_masc_pcr_primers_given_variants


TEST_INPUT_DIR = os.path.join(PACKAGE_DIR, 'tests', 'test_input')

TEST_OUTPUT_DIR = os.path.join(PACKAGE_DIR, 'tests', 'test_output')

TEST_MASC_PCR_INPUT = os.path.join(
        TEST_INPUT_DIR, 'fix_recoli_exp5_variants.csv')

TEST_REF_GENOME = os.path.join(TEST_INPUT_DIR, 'mg1655.NC_000913.2.fa')

TEST_OUTPUT_CSV = os.path.join(
        TEST_OUTPUT_DIR, 'fix_recoli_exp5_mascpcr_primers.csv')

find_masc_pcr_primers_given_variants(
        TEST_MASC_PCR_INPUT, TEST_REF_GENOME, TEST_OUTPUT_CSV)
