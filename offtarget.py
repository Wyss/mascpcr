import os
import multiprocessing

import numpy as np
import primer3

from libnano.fileio.fasta import parseFasta
from libnano import seqstr


def checkOffTarget(primer, genome, primer_idx, genome_rc=None,
                   hamming_percentile=0.05):
    ''' Check for the strongest off-target hybridization of `primer`
    on both the fwd and rev strands of `genome`.

    The 5' index of the primer on the fwd strand must be provided for masking
    purposes. The reverse complement of the genome, `genome_rc`, may be
    provided as a performance optimization.

    `hamming_percentile` specifies the percentile cutoff for hamming distance
    matches that will be screened with a thermodynamic alignment.
    '''
    genome_rc = genome_rc or seqstr.reverseComplement(genome)
    primer_length = len(primer)

    strand_results = multiprocessing.Queue()

    def _fwdStrand():
        fwd_hamming_distances = seqstr.rollingHammingDistance(primer, genome_rc)
        fwd_hd_thresh = np.percentile(fwd_hamming_distances, hamming_percentile)
        fwd_primer_footprint = (-(primer_idx+primer_length), (-primer_idx))
        fwd_hamming_distances[fwd_primer_footprint[0]: \
                              fwd_primer_footprint[1]] = primer_length
        fwd_hotspots, = np.where((fwd_hamming_distances < fwd_hd_thresh))

        highest_tm_idx = None
        highest_tm = -100

        for idx in fwd_hotspots:
            tm = primer3.calcHeterodimerTm(
                primer, genome[-(idx+primer_length):-idx])
            if tm > highest_tm:
                highest_tm_idx = idx
                highest_tm = tm

        strand_results.put((highest_tm, highest_tm_idx, 1))

    def _revStrand():
        rev_hamming_distances = seqstr.rollingHammingDistance(primer, genome)
        rev_hd_thresh = np.percentile(rev_hamming_distances, hamming_percentile)
        rev_primer_footprint = ((primer_idx)-primer_length,
                                (primer_idx+primer_length))
        rev_hamming_distances[rev_primer_footprint[0]: \
                              rev_primer_footprint[1]] = primer_length
        rev_hotspots, = np.where((rev_hamming_distances < rev_hd_thresh))

        highest_tm_idx = None
        highest_tm = -100

        for idx in rev_hotspots:
            tm = primer3.calcHeterodimerTm(
                primer, genome_rc[idx:idx+primer_length])
            if tm > highest_tm:
                highest_tm_idx = idx
                highest_tm = tm

        strand_results.put((highest_tm, highest_tm_idx, 0))

    fwd_proc = multiprocessing.Process(target=_fwdStrand)
    rev_proc = multiprocessing.Process(target=_revStrand)
    fwd_proc.start()
    rev_proc.start()
    fwd_proc.join()
    rev_proc.join()

    res1 = strand_results.get()
    res2 = strand_results.get()

    return res1 if res1[1] > res2[1] else res2
# end def

def checkOffTarget3p(primer, genome, primer_idx, 
                    num_bases=10, genome_rc=None, thresholds=None):
    ''' Check a primer for 3 prime homology using a threshold for mismatches 
    in the  last 5 bases based on hamming distance of the last num_bases bases 
    from the 3 prime end
    '''
    genome_rc = genome_rc or seqstr.reverseComplement(genome)
    primer_length = len(primer)
    subprimer = primer[-num_bases:]

    if thresholds == None:
        if num_bases > 7:
            thresholds = [0 for x in range(num_bases)]
            thresholds[0] = -1
            thresholds[1] = -1
            thresholds[2] = 2
            thresholds[3] = 2
            thresholds[4] = 2
            thresholds[5] = 1
            thresholds[6] = 1
        else:
            raise ValueError("num_bases too low for defualt thresholds")

    strand_results = multiprocessing.Queue()

    def _fwdStrand():
        a = 1
        if num_bases < primer_idx:
            a = seqstr.can3pMisprime(subprimer, 
                                        genome_rc[:primer_idx],
                                        thresholds,
                                        1)
        b = seqstr.can3pMisprime(subprimer, 
                                    genome_rc[primer_idx+primer_length:],
                                    thresholds,
                                    1)
        strand_results.put(a and b)

    def _revStrand():
        a = 1
        if num_bases < primer_idx:
            a = seqstr.can3pMisprime(subprimer, 
                                        genome[:primer_idx],
                                        thresholds,
                                        1)
        b = seqstr.can3pMisprime(subprimer, 
                                    genome[primer_idx+primer_length:],
                                    thresholds,
                                    1)
        strand_results.put(a and b)

    fwd_proc = multiprocessing.Process(target=_fwdStrand)
    rev_proc = multiprocessing.Process(target=_revStrand)
    fwd_proc.start()
    rev_proc.start()
    fwd_proc.join()
    rev_proc.join()

    res1 = strand_results.get()
    res2 = strand_results.get()

    return res1 and res2
# end def

if __name__ == '__main__':
    import random
    import timeit
    import time
    import numpy as np


    GENOME = os.path.join(LOCAL_DIR,
        '2014_02_18_gen9_54_failed_seg_fixed_with_gc_fixes.fa')

    genome_seq = parseFasta(GENOME)[0][1]
    genome_rc = seqstr.reverseComplement(genome_seq)
    # primer_seq = ''.join(random.choice('ATGC') for x in range(20))
    primer_seq = seqstr.reverseComplement(genome_seq[10:31])
    # t = seqstr.rollingHammingDistance(primer_seq, genome_seq)

    # idx, res = checkOffTarget(primer_seq, genome_seq, 10, genome_rc)

    # checkOffTarget(primer_seq, genome_seq, 10)

    # print(seqstr.rollingHammingDistance(primer_seq, primer_seq))

    # print(primer_seq)
    # print(seqstr.reverseComplement(genome_seq[10:31]))
    ts = time.time()
    # This is an option to avoid memory copies but it's slow to init
    genome_seq_mp = multiprocessing.Array('c', genome_seq, lock=False)
    genome_rc_mp = multiprocessing.Array('c', genome_rc, lock=False)
    print(time.time() - ts)
    # print(checkOffTarget(primer_seq, genome_seq_mp, 10, genome_rc=genome_rc_mp))

    import cProfile as profile
    # profile.run('print(checkOffTarget(primer_seq, genome_seq, 10))')
    profile.run('print(checkOffTarget(primer_seq, genome_seq_mp, 10, genome_rc=genome_rc_mp))')

