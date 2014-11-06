'''
primerpartition | primerpartition.py
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Contains methods for partitioning PrimerCandidates into the MASC PCR bins.

The methods are strand specific (partitionCandidatesFWD partitions forward 
strand candidates and partitionCandidatesREV partitions reverse strand
candidates)

'''

import primer3

import primercandidate


def partitionCandidatesFWD(
            primer_candidates,      # List of PCs ordered by idx
            start_idx,              # Start idx for MASC PCR design
            end_idx,                # End idx for MASC PCR design
            idx_lut,                # LUT mapping genome idx to ref_genome idx
            genome_str,             # Genome sequence string
            ref_genome_str,         # Reference genome sequence string
            edge_lut,               # Pre-populated edge LUT
            mismatch_lut,           # Pre-populated mismatch LUT
            params):

    spurious_tm_clip =      params['spurious_tm_clip']
    product_sizes =         params['product_sizes']
    size_tol =              params['product_size_tolerance']

    bin_size = (end_idx-start_idx)/len(product_sizes)
    size_range_iter = iter(product_sizes)
    masc_pc_bins = [[] for x in product_sizes]
    current_idx = start_idx
    current_size = next(size_range_iter)

    i = 0
    # primer_candidates already sorted by index from 5' - 3' on fwd strand
    for discriminatory_primer, wt_primer in primer_candidates:
        discriminatory_primer_end5p = discriminatory_primer.idx
        common_primer_end3p = (discriminatory_primer_end5p + current_size -
                               size_tol)
        if common_primer_end3p > current_idx + bin_size:
            if common_primer_end3p > end_idx:
                break
            if (current_idx < discriminatory_primer_end5p < current_idx +
                bin_size):
                # only one end out of bounds
                continue

            current_idx += bin_size
            try:
                current_size = next(size_range_iter)
            except StopIteration:
                continue
            common_primer_end3p = (discriminatory_primer_end5p + current_size -
                                   size_tol)
            i += 1
            if current_idx >= end_idx:
                break

        # check if edge in the middle of PCR product
        if (edge_lut[discriminatory_primer_end5p:
                     common_primer_end3p].count(1) == 0):
            # find the limit with no edges in it
            lim3p = min(common_primer_end3p + 2*size_tol + 1, end_idx + 1)
            max3p_no_edge = common_primer_end3p
            for max3p_no_edge in range(common_primer_end3p + 1, lim3p):
                if edge_lut[max3p_no_edge] != 0:
                    max3p_no_edge -= 1
                    break

            rev_candidate_best = None
            for j in range(common_primer_end3p, max3p_no_edge + 1):
                rev_candidate = primercandidate.findCommonPrimer(
                    j, -1, idx_lut, genome_str, ref_genome_str, edge_lut, mismatch_lut,
                    params)
                if rev_candidate is not None:
                    if primer3.calcHeterodimer(
                            discriminatory_primer.seq,
                            rev_candidate.seq,
                            **params['thermo_params']).tm > spurious_tm_clip:
                        break
                    if rev_candidate_best is None:
                        rev_candidate_best = rev_candidate
                    else:
                        if rev_candidate_best.score < rev_candidate.score:
                            rev_candidate_best = rev_candidate

            if rev_candidate_best is not None:
                masc_pc_bins[i].append((discriminatory_primer,
                                    wt_primer,
                                    rev_candidate_best))
    return masc_pc_bins


def partitionCandidatesREV(primer_candidates,
                           start_idx, end_idx,
                           idx_lut,
                           genome_str, ref_genome_str,
                           edge_lut, mismatch_lut,
                           params=None):
    '''
    primer_candidates: ordered by lowest to highest index of primer_candidates
    priming points
    '''
    spurious_tm_clip =      params['spurious_tm_clip']
    product_sizes =         params['product_sizes']
    size_tol =              params['product_size_tolerance']

    bin_size = (end_idx-start_idx)/len(product_sizes)
    size_range_iter = iter(product_sizes)
    masc_pc_bins = [[] for x in product_sizes]
    current_idx = start_idx
    current_size = next(size_range_iter)

    i = 0
    # primer_candidates already sorted by index from 5' - 3' on fwd strand
    for discriminatory_primer, wt_primer in primer_candidates:
        discriminatory_primer_end5p = (discriminatory_primer.idx +
                                       discriminatory_primer.length)
        common_primer_end3p = (discriminatory_primer_end5p - current_size +
                               size_tol)

        if discriminatory_primer_end5p > current_idx + bin_size:
            if discriminatory_primer_end5p > end_idx:
                break
            if  current_idx < common_primer_end3p < current_idx + bin_size:
                # one end out of bounds
                continue

            current_idx += bin_size
            try:
                current_size = next(size_range_iter)
            except StopIteration:
                continue
            common_primer_end3p = (discriminatory_primer_end5p - current_size +
                                   size_tol)
            i += 1
            if current_idx >= end_idx:
                break

        # check if edge in the middle of PCR product
        if common_primer_end3p > current_idx:
            # check if edge in the middle
            if (edge_lut[common_primer_end3p:
                         discriminatory_primer_end5p+1].count(1) == 0):
                # find the limit with no edges in it
                lim3p = max(common_primer_end3p-2*size_tol - 1, start_idx - 1)
                min3p_no_edge = common_primer_end3p
                for min3p_no_edge in range(common_primer_end3p-1, lim3p, -1):
                    if edge_lut[min3p_no_edge] != 0:
                        min3p_no_edge += 1
                        break

                fwd_candidate_best = None
                primercandidate.mismatch_count = 0
                for j in range(common_primer_end3p, min3p_no_edge - 1, -1):
                    fwd_candidate = primercandidate.findCommonPrimer(
                        j, 1, idx_lut, genome_str, ref_genome_str, edge_lut,
                        mismatch_lut, params)
                    if fwd_candidate is not None:
                        if (primer3.calcHeterodimer(
                                discriminatory_primer.seq,
                                fwd_candidate.seq,
                                **params['thermo_params']).tm > 
                                spurious_tm_clip):
                            break
                        if fwd_candidate_best is None:
                            fwd_candidate_best = fwd_candidate
                        else:
                            if fwd_candidate_best.score < fwd_candidate.score:
                                fwd_candidate_best = fwd_candidate
                if fwd_candidate_best is not None:
                    masc_pc_bins[i].append((discriminatory_primer,
                                        wt_primer,
                                        fwd_candidate_best))
    return masc_pc_bins
