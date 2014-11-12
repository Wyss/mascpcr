"""
design_seg23 | examples/design_seg23.py
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This is an example script that uses the MASC PCR design pipeline to design
primers for a single inserted segment ("seg23").

For one-off design tasks you might be better off using the `mascpcrcli` 
script, but you could also modify this script to suit your needs.

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


# We will set up our params here, note that we only define parameters that 
# deviate from the defaults (see mascpcr/pipeline.py for defaults)
params = {
    # Our output files and sequence names will be prefixed with "seg23"
    'output_basename': 'example_seg23'      
}

# Let's collect + build the filepaths to our genbank files. You can get away
# with using relative file paths or hard coded filepaths in your scripts,
# but we will build them dynamically here so that this script will work
# on any system.
RECODED_GENOME_FP = os.path.join(PACKAGE_DIR, 'tests', 'test_input', 
                                 'recoded_genome.gb')
REFERENCE_GENOME_FP = os.path.join(PACKAGE_DIR, 'tests', 'test_input', 
                                   'reference_genome.gb')

# Let's pretend that we don't know the start and end indices of the segment
# for which we are designing primers. You'd probably want to run this 
# separately and check the output, but for the purposes of this example I'll
# just do it inline. 

from Bio import SeqIO  # Need SeqIO to build SeqRecord object

RECODED_GB_DATA = SeqIO.read(RECODED_GENOME_FP, 'gb')
REFERENCE_GB_DATA = SeqIO.read(REFERENCE_GENOME_FP, 'gb')

start_idx, end_idx = mascpcr.genbankfeatures.findAggregateBoundaries(
    # `sr_obj` expects a SeqRecord object
    sr_obj=RECODED_GB_DATA,
    # These are the features types that we want to look for
    feature_types=['synth_fragment'],
    # We will look up the qualifier 'label' in each Feature object and check
    # to make sure that the regular expression "seg23.*" matches its contents
    # (regex refresher: seg23.* will match "seg23" followed by any characters,
    #  e.g., seg23_001)
    qualifier_regexs={'label':'seg23.*'}
)

# As I already know the expected start and end indicies of seg23, we will 
# just assert that they are correct here as proof of principle

assert start_idx == 1039203
assert end_idx == 1084477

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

# This call runs the actual pipeline with the data structures that we just
# generated. Note that it requires the genome string and reference genome
# string, which we get from the SeqRecord object (need to cast it to a str)
mascpcr.findMascPrimers(
    idx_lut=idx_lut,
    genome_str=str(RECODED_GB_DATA.seq),
    ref_genome_str=str(REFERENCE_GB_DATA.seq),
    start_idx=start_idx,
    end_idx=end_idx,
    edge_lut=edge_lut,
    mismatch_lut=mismatch_lut,
    border_lut=border_lut,
    params=params
)
