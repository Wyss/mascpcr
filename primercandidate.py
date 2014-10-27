import os

import numpy as np

from collections import namedtuple

from primer3 import bindings as p3b

LOCAL_DIR = os.path.dirname(os.path.realpath(__file__))

from libnano import seqstr


CandidatePrimer = namedtuple('CandidatePrimer',
        ['idx',                 # Index of the 5' most nt on the fwd strand
                                #    of the genome
         'seq',                 # 5' to 3' sequence of the primer
         'strand',
         'length',              # Length in bp
         'mismatch_idxs',       # Indices of mismatches from the 3' end
         'tm',                  # Melting temperature of the primer
         'tm_homo',             # Homodimer melting temperature
         'tm_hairpin',          # Hairpin melting temperature
         'score'                # Weighted score based on discriminatory
         ])                     #   power, homo/hairpin tm, etc.

'''
Primer indexing information

The primer index provided in the CandidatePrimer tuple is the 5'-most index
of the primer footprint on the fwd (+) strand of the genome / target:

Fwd Primer Idx         |
Fwd Primer             >>>>>>>>>>>>>>>>>>
Genome          ATTACCGATACCAATTGACCAGTTGGGACCCAGTTGACCAGTTGGACCCAGTTAGC
Rev Primer                                       <<<<<<<<<<<<<<<<<<<
Rev Primer Idx                                   |

'''


def generateWtPrimer(d_idx, d_length, d_strand, index_lut, ref_genome, params):
    wt_seq = ref_genome[index_lut[d_idx]:index_lut[d_idx + d_length]]
    return CandidatePrimer(
        idx=index_lut[d_idx],
        seq=wt_seq,
        strand=d_strand,
        length=d_length,
        mismatch_idxs=[0] * d_length,
        tm=p3b.calcTm(wt_seq, **params['thermo_params']),
        tm_homo=p3b.calcHomodimer(wt_seq, **params['thermo_params']).tm,
        tm_hairpin=p3b.calcHairpin(wt_seq, **params['thermo_params']).tm,
        score=0
    )


def findCandidatePrimer(index, index_lut, genome, ref_genome, edges_lut,
                        mismatches_lut, params, fwd_strand=True,
                        disallow_mismatches=False, check_wt_primer=False):

    strand_sign = 1 if fwd_strand else -1

    tm_range = params.get('tm_range', (60, 65))
    spurious_tm_clip = params.get('spurious_tm_clip', 40)
    size_range = params.get('size_range', (18, 30))
    thermo_params = params.get('thermo_params')
    mv_conc = thermo_params.get('mv_conc', 50)
    dv_conc = thermo_params.get('dv_conc', 0)
    dntp_conc = thermo_params.get('dntp_conc', 0.8)
    dna_conc = thermo_params.get('dna_conc', 50)
    temp_c = thermo_params.get('temp_c', 37)
    max_loop = thermo_params.get('max_loop', 30)

    # Mismatch score weights from the 3' end of the primer
    weight_mismatch_by_idx = params.get('mismatch_weights',
                                        [5, 4, 4, 3, 3, 2, 1])

    # Score calculation functions for thermodynamic interactions
    f_weight_homo_dim_tm =  lambda x: -(1.0 / spurious_tm_clip)*x
    f_weight_hairpin_dim_tm = lambda x: -(1.0 / spurious_tm_clip)*x
    f_weight_hetero_dim_tm = lambda x: -2.0 * abs(x - 0.5 * (tm_range[0] +
                                                  tm_range[1]))/(tm_range[1] -
                                                  tm_range[0])

    if index - size_range[1] < 0 or index + size_range[1] > (len(genome)-1):
        return None

    # Initialize
    delta = size_range[0] - 1
    delta_lim = size_range[1] + 1
    primer = None 
    prev_score = -1000  

    if fwd_strand:
        root_idx  = index + 1  # + 1 for exclusive upper-bound indexing
        candidate_seq_area = genome[root_idx - delta_lim:root_idx]
        # If mismatches are not allowed, check the initial footprint
        if disallow_mismatches:
            if np.sum(mismatches_lut[root_idx - delta:root_idx]) > 0:
                # printPrimer(index, delta, index_lut, genome, ref_genome)
                return None
    else:
        root_idx = index
        candidate_seq_area = seqstr.reverseComplement(
                                genome[root_idx:root_idx + delta_lim])
        # If mismatches are not allowed, check the initial footprint
        if disallow_mismatches:
            if np.sum(mismatches_lut[root_idx:root_idx + delta]) > 0:
                # printPrimer(index, -delta, index_lut, genome, ref_genome)
                return None

    # Check the 3' end for GC content out 5 bases
    end3p_check = candidate_seq_area[-5:]
    if end3p_check.count(b"G") + end3p_check.count(b"C") > 3:
        return None

    # Candidate seq area is the total potential primer sequence from 5' to 3',
    # which we examine from the 3' end to the 5' end
    wt_primer = None
    while delta < delta_lim:
        delta += 1

        if disallow_mismatches:
            if mismatches_lut[root_idx - strand_sign * delta] == 1:
                return primer

        candidate = candidate_seq_area[-delta:]

        primer_tm = p3b.calcTm(
                seq=candidate,
                mv_conc=mv_conc,
                dv_conc=dv_conc,
                dntp_conc=dntp_conc,
                dna_conc=dna_conc,
                )

        if primer_tm < tm_range[0]:
            continue
        elif primer_tm > tm_range[1]:
            return primer

        hrp_res = p3b.calcHairpin(candidate,
                mv_conc=mv_conc,
                dv_conc=dv_conc,
                dntp_conc=dntp_conc,
                dna_conc=dna_conc,
                temp_c=temp_c,
                max_loop=max_loop)
        hrp_tm = hrp_res.tm if hrp_res is not None else 0.0
        if hrp_tm > spurious_tm_clip:
            return primer

        homo_res = p3b.calcHomodimer(candidate,
                mv_conc=mv_conc,
                dv_conc=dv_conc,
                dntp_conc=dntp_conc,
                dna_conc=dna_conc,
                temp_c=temp_c,
                max_loop=max_loop)
        homo_tm = homo_res.tm if homo_res is not None else 0.0
        if homo_tm > spurious_tm_clip:
            return primer

        wt_primer = None
        if check_wt_primer:
            g_idx = index - delta + 1 if fwd_strand else index
            wt_primer = generateWtPrimer(g_idx, delta, strand_sign, index_lut, ref_genome, params)
            if wt_primer.tm > tm_range[1]:
                return None
            elif wt_primer.tm < tm_range[0]:
                continue
            if max((wt_primer.tm_homo, wt_primer.tm_hairpin)) > spurious_tm_clip:
                return None

        local_mismatch_idxs = [0] * delta
        local_mismatch_score = 0

        if not disallow_mismatches:
            weight_mismatch_by_idx_len = len(weight_mismatch_by_idx)
            num_mismatches = 0
            for i in range(0, delta):
                i = i * strand_sign
                if edges_lut[index - i * strand_sign] == 1:
                    return primer
                if mismatches_lut[index - i * strand_sign] == 1:
                    local_mismatch_idxs[i] = 1
                    num_mismatches += 1
                    if abs(i) < weight_mismatch_by_idx_len:
                        local_mismatch_score += weight_mismatch_by_idx[i]
                    else:
                        local_mismatch_score += weight_mismatch_by_idx[-1]

            if num_mismatches < params.get('min_num_mismatches', 1):
                continue

        score = (f_weight_hetero_dim_tm(primer_tm) +
                 f_weight_hairpin_dim_tm(hrp_tm) +
                 f_weight_homo_dim_tm(homo_tm) +
                 local_mismatch_score)

        if score > prev_score:
            primer = CandidatePrimer(
                idx=index - delta + 1 if fwd_strand else index,
                seq=candidate,
                strand=strand_sign,
                length=delta,
                mismatch_idxs=local_mismatch_idxs,
                tm=primer_tm,
                tm_homo=homo_tm,
                tm_hairpin=hrp_tm,
                score=score
            )
            prev_score = score
    # end while
    return primer
