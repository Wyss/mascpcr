import sys
import os
import re

LOCAL_DIR = os.path.dirname(os.path.realpath(__file__))

from libnano import util, seqstr, seqfilter
from libnano.fileio.fasta import parseFasta

from genbankfeatures import filterFeatures
from Bio import SeqIO

def findSegBoundaries(gb_rec, seg_num):
    g9_frags = filterFeatures(gb_rec, 'gen9_seg')
    min_idx = 1000000000000
    max_idx = -1
    for frag in g9_frags:
        frag_num = int(re.search('seg(\d+)_.*', str(frag.qualifiers['label'])).group(1))
        if frag_num == seg_num:
            min_idx = min(int(frag.location.start), min_idx)
            max_idx = max(int(frag.location.end), max_idx)
    return min_idx, max_idx

if __name__ == '__main__':
    # GENOME = os.path.join(LOCAL_DIR, 'test_data', '2014_02_18_gen9_54_failed_seg_fixed_with_gc_fixes.gb')
    # genome_rec = rec = SeqIO.read(GENOME, 'gb')
    # frag_idxs = {}
    # for frag in range(5,49):
    #     frag_idxs[frag] = findSegBoundaries(genome_rec, frag)
    # print frag_idxs
    GENOME = os.path.join(LOCAL_DIR, 'test_data', '2014_06_25_sgi_seg59_to_seg68_v2.gb')
    genome_rec = rec = SeqIO.read(GENOME, 'gb')
    frag_idxs = {}
    for frag in range(59,69):
        frag_idxs[frag] = findSegBoundaries(genome_rec, frag)
    print frag_idxs
    GENOME = os.path.join(LOCAL_DIR, 'test_data', '2014_07_30_gen9_seg0_to_seg4.gb')
    genome_rec = rec = SeqIO.read(GENOME, 'gb')
    frag_idxs = {}
    for frag in range(0,5):
        frag_idxs[frag] = findSegBoundaries(genome_rec, frag)
    print frag_idxs