'''
primercandidate | primercandidate.py
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Methods for finding discriminatory and common primer candidates given
a target index and the necessary genome information and LUTs.

'''

import os

from collections import namedtuple

import primer3

LOCAL_DIR = os.path.dirname(os.path.realpath(__file__))

from libnano import seqstr

# ~~~~~~~~~~~~~~~~~~~~~ Function aliases for convenience ~~~~~~~~~~~~~~~~~~~~ #

rc = seqstr.reverseComplement
caclHetero = primer3.calcHeterodimer
calcHomo = primer3.calcHomodimer
calcHrp = primer3.calcHairpin
calcTm = primer3.calcTm


CandidatePrimer = namedtuple('CandidatePrimer',
        ['idx',                 # Index of the 5' most nt on the fwd strand
                                #   of the genome
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


def findDiscriminatoryPrimer(idx, strand, idx_lut, genome, ref_genome,
                             edges_lut, mismatches_lut, params):

    # Initial param reconciliation
    tm_range = params.get('tm_range', (60, 65))
    spurious_tm_clip = params.get('spurious_tm_clip', 40)
    size_range = params.get('size_range', (18, 30))
    thermo_params = params.get('thermo_params')
    min_num_mismatches = params.get('min_num_mismatches', 1)

    # Mismatch score weights from the 3' end of the primer
    weight_mismatch_by_idx = params.get('mismatch_weights',
                                        [5, 4, 4, 3, 3, 2, 1])

    # Score calculation functions for thermodynamic interactions
    f_weight_homo_dim_tm = lambda x: -(1.0 / spurious_tm_clip) * x
    f_weight_hairpin_dim_tm = lambda x: -(1.0 / spurious_tm_clip) * x
    f_weight_hetero_dim_tm = lambda x: -2.0 * abs(x - 0.5 * (tm_range[0] +
                                                  tm_range[1]))/(tm_range[1] -
                                                  tm_range[0])

    if idx - size_range[1] < 0 or idx + size_range[1] > (len(genome)-1):
        return None, None

    delta = size_range[0] - 1
    delta_lim = size_range[1] + 1
    prev_score = -1000

    # Candidate seq area is the total potential primer sequence from 5' to 3',
    # which we examine from the 3' end to the 5' end
    if strand == 1:
        root_idx  = idx + 1  # + 1 for exclusive upper-bound idxing
        candidate_seq_area = genome[root_idx - delta_lim:root_idx]
        wt_candidate_seq_area = ref_genome[idx_lut[root_idx] -
                                       delta_lim:idx_lut[root_idx]]
    else:
        root_idx = idx
        candidate_seq_area = rc(genome[root_idx:root_idx + delta_lim])
        wt_candidate_seq_area = rc(ref_genome[idx_lut[root_idx]:
                                              idx_lut[root_idx + delta_lim]])

    # Check the 3' end for high end stability
    end3p_check = candidate_seq_area[-5:]
    if end3p_check.count(b'G') + end3p_check.count(b'C') > 3:
        return None, None

    best_primer_pair = (None, None)
    while delta < delta_lim:
        delta += 1
        # This is the putative primer sequence given the current delta
        primer_candidate = candidate_seq_area[-delta:]
        wt_primer_candidate = wt_candidate_seq_area[-delta:]

        # ~~~~~~~~~~~ Check for compliance with hard thermo limits ~~~~~~~~~~ #
        primer_tm = calcTm(seq=primer_candidate, **thermo_params)
        wt_primer_tm = calcTm(seq=wt_primer_candidate, **thermo_params)

        if primer_tm < tm_range[0] or wt_primer_tm < tm_range[0]:
            continue
        elif primer_tm > tm_range[1] or wt_primer_tm > tm_range[1]:
            break

        hrp_tm = calcHrp(primer_candidate, **thermo_params).tm
        homo_tm = calcHomo(primer_candidate, **thermo_params).tm
        wt_hrp_tm = calcHrp(wt_primer_candidate, **thermo_params).tm
        wt_homo_tm = calcHomo(wt_primer_candidate, **thermo_params).tm

        if hrp_tm > spurious_tm_clip or homo_tm > spurious_tm_clip:
            break

        if wt_hrp_tm > spurious_tm_clip or wt_homo_tm > spurious_tm_clip:
            break

        # ~~~~~~~~~~~~~~~ Score putative discriminatory primer ~~~~~~~~~~~~~~ #
        local_mismatch_idxs = [0] * delta
        local_mismatch_score = 0

        weight_mismatch_by_idx_len = len(weight_mismatch_by_idx)
        num_mismatches = 0
        for i in range(0, delta):
            i = i * strand
            # If primer now contains an edge, scrap it and stop looking
            if edges_lut[idx - i * strand] == 1:
                break
            if mismatches_lut[idx - i * strand] == 1:
                local_mismatch_idxs[i] = 1
                num_mismatches += 1
                if abs(i) < weight_mismatch_by_idx_len:
                    local_mismatch_score += weight_mismatch_by_idx[i]
                else:
                    local_mismatch_score += weight_mismatch_by_idx[-1]

        # Insure that primer contains the minimum number of mismatches
        if num_mismatches < min_num_mismatches:
            continue

        score = (f_weight_hetero_dim_tm(primer_tm) +
                 f_weight_hairpin_dim_tm(hrp_tm) +
                 f_weight_homo_dim_tm(homo_tm) +
                 local_mismatch_score)

        if score > prev_score:
            primer = CandidatePrimer(
                idx=idx - delta + 1 if strand == 1 else idx,
                seq=primer_candidate,
                strand=strand,
                length=delta,
                mismatch_idxs=local_mismatch_idxs,
                tm=primer_tm,
                tm_homo=homo_tm,
                tm_hairpin=hrp_tm,
                score=score
            )
            wt_primer = CandidatePrimer(
                idx=idx_lut[idx - delta + 1 if strand == 1 else idx],
                seq=wt_primer_candidate,
                strand=strand,
                length=delta,
                mismatch_idxs=[0] * delta,
                tm=wt_primer_tm,
                tm_homo=wt_homo_tm,
                tm_hairpin=wt_hrp_tm,
                score=0
            )
            best_primer_pair = (primer, wt_primer)
        prev_score = score
    return best_primer_pair


def findCommonPrimer(idx, strand, idx_lut, genome, ref_genome,
                     edges_lut, mismatches_lut, params):

    # Initial param reconciliation
    tm_range = params.get('tm_range', (60, 65))
    spurious_tm_clip = params.get('spurious_tm_clip', 40)
    size_range = params.get('size_range', (18, 30))
    thermo_params = params.get('thermo_params')

    # Score calculation functions for thermodynamic interactions
    f_weight_homo_dim_tm = lambda x: -(1.0 / spurious_tm_clip) * x
    f_weight_hairpin_dim_tm = lambda x: -(1.0 / spurious_tm_clip) * x
    f_weight_hetero_dim_tm = lambda x: -2.0 * abs(x - 0.5 * (tm_range[0] +
                                                  tm_range[1]))/(tm_range[1] -
                                                  tm_range[0])

    if idx - size_range[1] < 0 or idx + size_range[1] > (len(genome)-1):
        return None

    delta = size_range[0] - 1
    delta_lim = size_range[1] + 1
    prev_score = -1000

    # Candidate seq area is the total potential primer sequence from 5' to 3',
    # which we examine from the 3' end to the 5' end
    if strand == 1:
        root_idx  = idx + 1  # + 1 for exclusive upper-bound idxing
        candidate_seq_area = genome[root_idx - delta_lim:root_idx]
    else:
        root_idx = idx
        candidate_seq_area = rc(genome[root_idx:root_idx + delta_lim])

    # Check the 3' end for high end stability
    end3p_check = candidate_seq_area[-5:]
    if end3p_check.count(b'G') + end3p_check.count(b'C') > 3:
        return None

    best_common_primer = None
    while delta < delta_lim:
        delta += 1

        if mismatches_lut[root_idx - strand * delta] == 1:
            return best_common_primer

        primer_candidate = candidate_seq_area[-delta:]

        # ~~~~~~~~~~~ Check for compliance with hard thermo limits ~~~~~~~~~~ #
        primer_tm = calcTm(seq=primer_candidate, **thermo_params)

        if primer_tm < tm_range[0]:
            continue
        elif primer_tm > tm_range[1]:
            break

        hrp_tm = calcHrp(primer_candidate, **thermo_params).tm
        homo_tm = calcHomo(primer_candidate, **thermo_params).tm

        if hrp_tm > spurious_tm_clip or homo_tm > spurious_tm_clip:
            break

        score = (f_weight_hetero_dim_tm(primer_tm) +
                 f_weight_hairpin_dim_tm(hrp_tm) +
                 f_weight_homo_dim_tm(homo_tm))

        if score > prev_score:
            best_common_primer = CandidatePrimer(
                idx=idx - delta + 1 if strand == 1 else idx,
                seq=primer_candidate,
                strand=strand,
                length=delta,
                mismatch_idxs=[0] * delta,
                tm=primer_tm,
                tm_homo=homo_tm,
                tm_hairpin=hrp_tm,
                score=score
            )
            prev_score = score
    return best_common_primer
