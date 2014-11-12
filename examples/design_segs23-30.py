"""
design_seg23 | examples/design_seg23-30.py
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This is an example script that uses the MASC PCR design pipeline to design
primers for multiple inserted segments ("seg23"-"seg30").

You can use the command line interface in a similar manner (using a
shell script for calls) but this will be faster as the genbank files
only need to be read in and indexed once. 

"""
from __future__ import print_function  # Python 2/3 compatibility

import os
import sys

# This fanciness is only necessary to run the script in place without
# installing the pipeline. If you install the pipeline you can just import
# "mascpcr" and call it a day.
PACKAGE_DIR = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
MODULE_DIR = os.path.join(PACKAGE_DIR, 'mascpcr')

try: 
    import mascpcr
except:
    if not 'mascpcr' in os.listdir(PACKAGE_DIR):
        raise IOError('`mascpcr` must be installed in your PYTHONPATH or '
                      'script must be run from within the package directory')
    sys.path.append(MODULE_DIR)
    sys.path.append(PACKAGE_DIR)
    import mascpcr


# Let's collect + build the filepaths to our genbank files. You can get away
# with using relative file paths or hard coded filepaths in your scripts,
# but we will build them dynamically here so that this script will work
# on any system.
RECODED_GENOME_FP = os.path.join(PACKAGE_DIR, 'tests', 'test_input', 
                                 'recoded_genome.gb')
REFERENCE_GENOME_FP = os.path.join(PACKAGE_DIR, 'tests', 'test_input', 
                                   'reference_genome.gb')

# Let's pretend that we don't know the start and end indices of the segments
# for which we are designing primers. You'd probably want to run this 
# separately and check the output, but for the purposes of this example I'll
# just do it inline. 

from Bio import SeqIO  # Need SeqIO to build SeqRecord object

RECODED_GB_DATA = SeqIO.read(RECODED_GENOME_FP, 'gb')
REFERENCE_GB_DATA = SeqIO.read(REFERENCE_GENOME_FP, 'gb')

seg_indices = {}

for seg_num in range(23, 31):

    start_idx, end_idx = mascpcr.genbankfeatures.findAggregateBoundaries(
        # `sr_obj` expects a SeqRecord object
        sr_obj=RECODED_GB_DATA,
        # These are the features types that we want to look for
        feature_types=['synth_fragment'],
        # We will look up the qualifier 'label' in each Feature object and 
        # checkto make sure that the regular expression "seg23.*" matches its 
        # contents(regex refresher: seg23.* will match "seg23" followed by any 
        # characters, e.g., seg23_001)
        qualifier_regexs={'label':'seg%d.*' % seg_num}
    )

    seg_indices[seg_num] = (start_idx, end_idx)

# Now let's build all of the necessary data structures for the design pipline
#   idx_lut - lookup table that maps recoded indices to ref indices
#   edge_lut - lookup table that indicates "edges" (discontinuities) in the 
#              index mapping
#   mismatch_lut - lookup table that indicates single bp mismatches between
#                  the recoded and reference genomes
#   border_lut - lookup table that indicates feature borders for the purpose
#                of designing primers that span assembly junctions
genome_str, ref_genome_str, idx_lut, edge_lut, mismatch_lut, \
    border_lut = mascpcr.generateLUTs(
        genome_fp=RECODED_GENOME_FP,
        ref_genome_fp=REFERENCE_GENOME_FP,
        start_idx=start_idx,
        end_idx=end_idx,          
        # We want our primers to span "synth_fragment" junctions, if possible 
        # to insure that our subassemblies all worked properly 
        border_feature_types=['synth_fragment'],
        # We won't cache the lookup tables for this demonstration, but if you
        # are planning on running the pipeline against a single recoded 
        # genome multiple times it's worth the ~20 mb or so to cache the 
        # lookup tables and save a couple of minutes per call
        cache_luts=False
)

# Now let's call the actual pipeline for each of the respective segments. We 
# will modify the `output_basename` argument of `params` so that each output
# file has an appropriate, respective name

for seg_num in range(23, 31):

    params = {'output_basename': 'example_seg%d' % seg_num}

    mascpcr.findMascPrimers(
        idx_lut=idx_lut,
        genome_str=str(RECODED_GB_DATA.seq),
        ref_genome_str=str(REFERENCE_GB_DATA.seq),
        start_idx=seg_indices[seg_num][0],
        end_idx=seg_indices[seg_num][1],
        edge_lut=edge_lut,
        mismatch_lut=mismatch_lut,
        border_lut=border_lut,
        params=params
    )
