'''
mascpcr | mascpcr.py

This is the main entry point for masc pcr design.

'''
from __future__ import print_function

import copy
import itertools
import multiprocessing as mp
import os
import sys

from collections import namedtuple

import primer3
from Bio import SeqIO

LOCAL_DIR = os.path.dirname(os.path.realpath(__file__))
CACHE_DIR = os.path.join(LOCAL_DIR, 'cache')

import indexing
import primercandidate
import genbankfeatures
import offtarget

from libnano import seqstr

try:
    CPU_COUNT = mp.cpu_count()
except NotImplementedError:
    CPU_COUNT = 2


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Default parameters for various parts of the design process.
#   ** Ranges are inclusive
#
# This dictionary is updated with user provided parameters, so the user does
# not need to provide an exhaustive / equivalent dictionary
# ~~

DEFAULT_PARAMS = {
    # Acceptable primer melting temperature range. The ideal melting
    # temperature will be considered to be at the center of this range
    'tm_range': (60, 65),
    # Acceptable primer size range in bp.
    'size_range': (18, 30),
    # Product sizes from 5' to 3' across the segment in bp
    'product_sizes': (700, 600, 500, 400, 300, 250, 200, 150),
    # 'product_sizes': (850, 700, 600, 500, 400, 300, 250, 200, 150, 100),
    # PCR product size tolerance (+/- this value in bp)
    'product_size_tolerance': 10,
    # minimum base distance of a primer to the boundaries of a segment
    'segment_edge_offset': 67,
    # The following parameters are related to thermodynamic calculations
    'thermo_params': {
        # [thermo] Monovalent cation concentration in mM
        'mv_conc': 50,
        # [thermo] Divalent cation concentration in mM
        'dv_conc': 1.5,
        # [thermo] dNTP concentration in mM
        'dntp_conc': 0.2,
        # [thermo] DNA concentration in nM
        'dna_conc': 200,
    },
    # Default output base fn / fp
    'output_fp_base': 'mascpcr_out',
    # Minimum number of mismatches for discriminatory primer
    'min_num_mismatches': 2,
}

MascPrimerSet = namedtuple('MascPrimerSet',
        ['score',           # Sum of the d_primer and c_primer scores
         'd_primer',        # Discriminatory primer
         'w_primer',        # Wild-type primer
         'c_primer'         # Common primer
         ])


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

# Command line / high level entry point for MASC PCR primer design
def doDesignPrimers(
                genome_gb,          # Genome sequence string
                ref_genome_gb,      # Reference genome sequence string
                start_idx,          # Start index for MASC PCR design
                end_idx,            # End index (inclusive) for MASC PCR design
                border_feature_types=None,
                border_feature_regexs=None,
                cache_luts=False,
                params=None):
    print('Reading in genome and reference genome genbank files...')
    genome_rec = SeqIO.read(genome_gb, 'gb')
    ref_genome_rec = SeqIO.read(ref_genome_gb, 'gb')
    genome_str = str(genome_rec.seq).encode('utf-8')
    ref_genome_str = str(ref_genome_rec.seq).encode('utf-8')

    # ~~~~~~~~~~~~~~~~~~~~~~~~~ Build LUTs as needed ~~~~~~~~~~~~~~~~~~~~~~~~ #

    cache_dir = None
    if cache_luts:
        cache_dir = CACHE_DIR
        try:
            os.makedirs(cache_dir)
        except OSError as e:
            if e.errno != 17:
                raise

    print('Building index lookup table with Mauve Aligner...')
    cache_dir = None
    if cache_luts:
        cache_dir = CACHE_DIR
        try:
            os.makedirs(cache_dir)
        except OSError as e:
            if e.errno != 17:
                raise
    idx_lut = indexing.buildIdxLUT(genome_gb, ref_genome_gb, genome_str,
                                   ref_genome_str, cache_dir=cache_dir)

    print('Building edge lookup table...')
    edge_lut = indexing.buildEdgeLUT(idx_lut, cache_dir=cache_dir)

    print('Building mismatch lookup table...')
    mismatch_lut = indexing.buildMismatchLUT(idx_lut, genome_str,
                                             ref_genome_str,
                                             cache_dir=cache_dir)

    border_arr = None
    if border_feature_types is not None:
        print('Finding border indices of feature type(s): {}...'.format(
              ','.join(border_feature_types)))
        border_arr = genbankfeatures.genBorderArray(genome_rec,
                                border_feature_types, border_feature_regexs)
        print('Found {} border indices'.format(border_arr.count(1)))

    return designPrimers(idx_lut, genome_str, ref_genome_str, start_idx,
                         end_idx, edge_lut, mismatch_lut, border_arr=border_arr,
                         cache_luts=cache_luts, params=params)


# Main entry point for MASC PCR primer design
def designPrimers(
            idx_lut,            # LUT mapping genome index to ref_genome index
            genome_str,         # Genome sequence string
            ref_genome_str,     # Reference genome sequence string
            start_idx,          # Start index for MASC PCR design
            end_idx,            # End index (inclusive) for MASC PCR design
            edge_lut,           # Pre-populated edge LUT
            mismatch_lut,       # Pre-populated mismatch LUT
            border_arr=None,    # Array of feature border indices to encompass
            cache_luts=False,   # Cache the LUTs w/ an md5 hash of the genome
            params=None):

    # ~~~~~~~~~~ Reconcile default params with user-provided params ~~~~~~~~~ #
    _params = copy.deepcopy(DEFAULT_PARAMS)
    _params.update(params or {})
    params = _params
    thermo_params = params.get('thermo_params')

    # ~~~~~~~ Generate potential primers for both fwd and rev strands ~~~~~~~ #
    fwd_strand_primer_candidates = []
    rev_strand_primer_candidates = []
    fwd_candidates_found = 0
    rev_candidates_found = 0
    candidates_checked = 0

    mismatch_count = mismatch_lut.count(1)
    print("Total Mismatches: {}, Percent Mismatches: {:.2f}".format(
          mismatch_count, mismatch_count/float(len(mismatch_lut))*100))
    print('Finding candidate primers...')

    for idx, is_mismatch in enumerate(mismatch_lut[start_idx:end_idx+1],
                                      start=start_idx):
        if is_mismatch:
            fwd_c, fwd_wt = primercandidate.findDiscriminatoryPrimer(
                                idx, 1, idx_lut, genome_str,
                                ref_genome_str, edge_lut,
                                mismatch_lut, params)
            rev_c, rev_wt = primercandidate.findDiscriminatoryPrimer(
                                idx, -1, idx_lut, genome_str,
                                ref_genome_str, edge_lut,
                                mismatch_lut, params)
            if fwd_c is not None:
                fwd_strand_primer_candidates.append((fwd_c, fwd_wt))
                fwd_candidates_found += 1
            if rev_c is not None:
                rev_strand_primer_candidates.append((rev_c, rev_wt))
                rev_candidates_found += 1
            candidates_checked += 1
            sys.stdout.write('Mismatches checked: {:5d} | Primer candidates'
                             ' found - fwd: {:3d}, rev: {:3d}\r'.format(
                             candidates_checked, fwd_candidates_found,
                             rev_candidates_found))
            sys.stdout.flush()
    print('Mismatches checked: {:5d} | Primer candidates found -> fwd: '
          '{:3d}, rev: {:3d}'.format(candidates_checked,
          fwd_candidates_found, rev_candidates_found))

    # ~~~~~~~~~~~~~~ Partition primer candidates into MASC bins ~~~~~~~~~~~~~ #
    print('Partitioning primers into fwd and rev bins...')
    fwd_bins_l2s = partitionCandidatesFWD(fwd_strand_primer_candidates,
                                          start_idx, end_idx,
                                          idx_lut,
                                          genome_str, ref_genome_str,
                                          edge_lut, mismatch_lut,
                                          params=DEFAULT_PARAMS,
                                          is_long_to_short=True)
    print("Fwd bin counts (l2s):")
    for i in range(len(fwd_bins_l2s)):
        print("Bin %d: %d" % (i, len(fwd_bins_l2s[i])))


    rev_bins_l2s = partitionCandidatesREV(rev_strand_primer_candidates,
                            start_idx, end_idx,
                            idx_lut,
                            genome_str, ref_genome_str,
                            edge_lut, mismatch_lut,
                            params=DEFAULT_PARAMS,
                            is_long_to_short=True)
    print("\nRev bin counts (l2s)")
    for i in range(len(rev_bins_l2s)):
        print("Bin %d: %d" % (i, len(rev_bins_l2s[i])))

    # ~~~~~~~~~~~~~~~~~~~~ Combine and score primer pairs ~~~~~~~~~~~~~~~~~~~ #
    print('Combining bins and scoring primer pairs...')
    combined_bins = [bf+br for bf, br in zip(fwd_bins_l2s, rev_bins_l2s)]
    for bin_idx, bin in enumerate(combined_bins):
        for pair_idx, pair in enumerate(bin):
            f_idx = pair[0].idx
            r_idx = pair[2].idx
            idx1, idx2 = (f_idx, r_idx) if f_idx < r_idx else (r_idx, f_idx)
            score = pair[0].score + pair[2].score
            if border_arr:
                # Add an arbitrarily high value to the score for each border
                # in the primer pair footprint
                score += (border_arr[idx1:idx2].count(1) * 50)
            combined_bins[bin_idx][pair_idx] = MascPrimerSet(score, pair[0],
                                                             pair[1], pair[2])
        # Sort each bin after scoring (highest score first)
        combined_bins[bin_idx].sort(key=lambda pset: -pset.score)

    print('Finding optimal set of MASC primers...')

    def checkSetForHeterodimers(set_of_primer_sets, tm_max=40):
        all_primers = []
        for bin_idx, primer_pair_idx in enumerate(set_of_primer_sets):
            primer_set_obj = combined_bins[bin_idx][primer_pair_idx]
            all_primers.append((bin_idx, primer_set_obj.d_primer))
            all_primers.append((bin_idx, primer_set_obj.w_primer))
            all_primers.append((bin_idx, primer_set_obj.c_primer))
        for p1, p2 in itertools.combinations(all_primers, 2):
            if (primer3.calcHeterodimer(p1[1].seq, p2[1].seq,
                                        **thermo_params).tm > tm_max):
                    return (False, p1[0], p2[0])
        return (True,)


    def checkSetForOffTarget(set_of_primer_sets, tm_max=40):
        for bin_idx, primer_pair_idx in enumerate(set_of_primer_sets):
            primer_set_obj = combined_bins[bin_idx][primer_pair_idx]
            d_primer = primer_set_obj.d_primer
            w_primer = primer_set_obj.w_primer
            c_primer = primer_set_obj.c_primer
            cached_tms = offtarget_tm_cache[bin_idx].get(primer_pair_idx)
            if cached_tms is not None:
                d_tm, w_tm, c_tm = cached_tms
            else:
                d_tm = offtarget.checkOffTarget(d_primer.seq, genome_str,
                                                 d_primer.idx)
                w_tm = offtarget.checkOffTarget(w_primer.seq, ref_genome_str,
                                                 w_primer.idx)
                c_tm = offtarget.checkOffTarget(c_primer.seq, genome_str,
                                                 c_primer.idx)
                offtarget_tm_cache[bin_idx][primer_pair_idx] = (d_tm, w_tm,
                                                                c_tm)
            if max(d_tm, w_tm, c_tm) > tm_max:
                print(d_tm, w_tm, c_tm)
                unacceptable_primers[bin_idx].append(primer_pair_idx)
                return False
        return True

    num_bins = len(combined_bins)
    bin_lengths = [len(bin) for bin in combined_bins]
    unacceptable_primers = [[] for _ in range(num_bins)]
    heterodimer_conflicts = [[[] for _ in range(bin_lengths[bin_idx])]
                             for bin_idx in range(num_bins)]
    offtarget_tm_cache = [{} for _ in range(num_bins)]

    set_of_primer_sets = [0] * num_bins


    def checkForHeterdimerConflicts(bin1_idx, set_of_primer_sets):
        ps1_idx = set_of_primer_sets[bin1_idx]
        for bin2_idx in range(num_bins):
            if bin1_idx == bin2_idx:
                continue
            ps2_idx = set_of_primer_sets[bin2_idx]
            if (bin2_idx, ps2_idx) in heterodimer_conflicts[bin1_idx][ps1_idx]:
                return False
        return True

    while True:
        for bin_idx in range(num_bins):
            idx_clean = False
            while not idx_clean:
                # Make sure that we aren't checking a primer set known to have
                # mispriming problems
                if len(unacceptable_primers[bin_idx]) == bin_lengths[bin_idx]:
                    print('Failed to find MASC PCR primer set due to total '
                          'mispriming of primer sets in bin {}'.format(bin_idx))
                while (set_of_primer_sets[bin_idx] in
                       unacceptable_primers[bin_idx]):
                    if set_of_primer_sets[bin_idx] < bin_lengths[bin_idx] - 1:
                        set_of_primer_sets[bin_idx] += 1
                    else:
                        set_of_primer_sets[bin_idx] = 0
                # Make sure that we don't have known primer compatibility
                # problems within the set
                while not checkForHeterdimerConflicts(bin_idx,
                                                      set_of_primer_sets):
                    if set_of_primer_sets[bin_idx] < bin_lengths[bin_idx] - 1:
                        set_of_primer_sets[bin_idx] += 1
                    else:
                        set_of_primer_sets[bin_idx] = 0
                else:
                    idx_clean = True

        print('Checking binned primer sets: {}'.format(set_of_primer_sets))
        hetero_res = checkSetForHeterodimers(set_of_primer_sets, 40)
        if not hetero_res[0]:
            bin1_idx, bin2_idx = hetero_res[1], hetero_res[2]
            ps1_idx = set_of_primer_sets[bin1_idx]
            ps2_idx = set_of_primer_sets[bin2_idx]
            print('Heterodimer conflict between {}.{}, {}.{}'.format(bin1_idx,
                    ps1_idx, bin2_idx, ps2_idx))
            try:
                bin1_delta = (combined_bins[bin1_idx][ps1_idx + 1][0] -
                              combined_bins[bin1_idx][ps1_idx][0])
            except IndexError:
                bin1_delta = -1
            try:
                bin2_delta = (combined_bins[bin2_idx][ps2_idx + 1][0] -
                              combined_bins[bin2_idx][ps2_idx][0])
            except IndexError:
                bin2_delta = -1
            if bin1_delta == -1:
                if bin2_delta == -1:
                    print('Failed to find MASC PCR primer set due to total '
                          'heterodimer incompatibility between bins {} and '
                          '{}'.format(bin1_idx, bin2_idx))
                    sys.exit(1)
                else:
                    set_of_primer_sets[bin2_idx] += 1
                    set_of_primer_sets[bin1_idx] = 0
            elif bin2_delta == -1:
                set_of_primer_sets[bin1_idx] += 1
                set_of_primer_sets[bin2_idx] = 0
            elif bin1_delta < bin2_delta:
                set_of_primer_sets[bin1_idx] += 1
                set_of_primer_sets[bin2_idx] = 0
            else:
                set_of_primer_sets[bin2_idx] += 1
                set_of_primer_sets[bin1_idx] = 0
            heterodimer_conflicts[bin1_idx][ps1_idx].append((bin2_idx, ps2_idx))
            heterodimer_conflicts[bin2_idx][ps1_idx].append((bin1_idx, ps1_idx))
            continue
        if not checkSetForOffTarget(set_of_primer_sets, 40):
            print('-> Failed mispriming check')
            continue
        else:
            break

    print('Found set of primer sets: ', set_of_primer_sets)

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~ Write CSV output ~~~~~~~~~~~~~~~~~~~~~~~~~~ #

    def getPrimerPairDelta(pc1, pc2):
        pc1_5p_idx = pc1.idx
        pc2_5p_idx = pc2.idx
        if pc1.strand == -1:
            pc1_5p_idx += pc1.length
        if pc2.strand == -1:
            pc2_5p_idx += pc2.length
        return abs(pc1.idx - pc2.idx)


    def writePrimerRec(fd, primer_obj, recoded_idx=True, gen_wt_idx=False):
        po = primer_obj
        if recoded_idx:
            if gen_wt_idx:
                po_idx = '{},{}'.format(po.idx, idx_lut[po.idx])
            else:
                po_idx = '{},N/A'.format(po.idx, idx_lut[po.idx])
        else:
            po_idx = 'N/A,{}'.format(po.idx)
        fd.write(','.join([str(x) for x in [
            po.strand, po_idx, po.length, po.seq, sum(po.mismatch_idxs),
            max(0, po.tm), max(0, po.tm_hairpin), max(0, po.tm_homo)]]) + '\n')

    # Get the output filename from the parameters dictionary. Default to
    # "mascpcr_out" in the local directory
    output_fp = params.get('output_fp_base', 'mascpcr_out') + '.csv'

    with open(output_fp, 'w') as fd_out:
        fd_out.write(','.join([
            'Bin number',
            'Potential pairs found',
            'Pair chosen',
            'Ideal PCR product size',
            'Actual PCR product size',
            'Pair score',
            'Included # of gen9 fragment junctions',
            'Discriminatory power (# mismatches)'
            ]) + '\n')
        for bin_idx, ps_idx in enumerate(set_of_primer_sets):
            primer_set = combined_bins[bin_idx][ps_idx]
            fp = primer_set.d_primer
            rp = primer_set.c_primer
            fd_out.write(','.join([str(x) for x in [
                bin_idx + 1,
                len(combined_bins[bin_idx]),
                ps_idx,
                params['product_sizes'][bin_idx],
                getPrimerPairDelta(fp, rp),
                primer_set[0],
                int(primer_set[0])/50,
                sum(fp.mismatch_idxs)
                ]]) + '\n')
        fd_out.write('\n' + ','.join([
            'Bin number',
            'Primer type',
            'Strand',
            'Recoded Index',
            'Wildtype Index',
            'Length',
            'Sequence',
            'Number of mismatches',
            'Tm (C)',
            'Hairpin Tm (C)',
            'Homodimer Tm (C)'
            ]) + '\n')
        for bin_idx, ps_idx in enumerate(set_of_primer_sets):
            primer_set = combined_bins[bin_idx][ps_idx]
            d_primer = primer_set.d_primer
            w_primer = primer_set.w_primer
            c_primer = primer_set.c_primer
            fd_out.write('{},Discriminatory primer,'.format(bin_idx+1))
            writePrimerRec(fd_out, d_primer)
            fd_out.write(',Wildtype primer,')
            writePrimerRec(fd_out, w_primer, recoded_idx=False)
            fd_out.write(',Common primer,')
            writePrimerRec(fd_out, c_primer, gen_wt_idx=True)
        fd_out.write('\nDESIGN PARAMETERS' + '\n' + '\n'.join(
            ['{},{}'.format(k, repr(v).replace(',', ';')) for k, v in
             params.items()]))


def partitionCandidatesFWD(discriminatory_primers,
                        start_idx, end_idx,
                        idx_lut,
                        genome_str, ref_genome_str,
                        edge_lut, mismatch_lut,
                        params=None,
                        is_long_to_short=True):
    """
    discriminatory_primers: ordered by lowest to highest index of discriminatory_primers
    priming points
    is_long_to_short: are the primers sets ordered longest to shortest
    across the fragment
    """
    tm_clip = 40
    product_sizes = params['product_sizes']
    bin_size = (end_idx-start_idx)/len(product_sizes)
    size_tol = params['product_size_tolerance']
    if is_long_to_short:
        size_range_iter = iter(product_sizes)

    out_bins = [[] for x in product_sizes]

    current_idx = start_idx
    current_size = next(size_range_iter)
    i = 0
    # discriminatory_primers already sorted by index from 5' - 3' on fwd strand
    for discriminatory_primer, wt_primer in discriminatory_primers:
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

            # increment bin and size
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
                            **params['thermo_params']).tm > tm_clip:
                        break
                    if rev_candidate_best is None:
                        rev_candidate_best = rev_candidate
                    else:
                        if rev_candidate_best.score < rev_candidate.score:
                            rev_candidate_best = rev_candidate

            if rev_candidate_best is not None:
                out_bins[i].append((discriminatory_primer,
                                    wt_primer,
                                    rev_candidate_best))
    return out_bins
# end def

def partitionCandidatesREV(discriminatory_primers,
                           start_idx, end_idx,
                           idx_lut,
                           genome_str, ref_genome_str,
                           edge_lut, mismatch_lut,
                           params=None,
                           is_long_to_short=True):
    """
    discriminatory_primers: ordered by lowest to highest index of discriminatory_primers
    priming points
    is_long_to_short: are the primers sets ordered longest to shortest
    across the fragment
    """
    tm_clip = 40
    product_sizes = params['product_sizes']
    bin_size = (end_idx-start_idx)/len(product_sizes)
    size_tol = params['product_size_tolerance']
    if is_long_to_short:
        size_range_iter = iter(product_sizes)
    else:
        size_range_iter = reversed(product_sizes)

    out_bins = [[] for x in product_sizes]

    current_idx = start_idx
    current_size = next(size_range_iter)
    i = 0
    for discriminatory_primer, wt_primer in discriminatory_primers:
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

            # increment bin and size
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
                        if primer3.calcHeterodimer(
                                discriminatory_primer.seq,
                                fwd_candidate.seq,
                                **params['thermo_params']).tm > tm_clip:
                            break
                        if fwd_candidate_best is None:
                            fwd_candidate_best = fwd_candidate
                        else:
                            if fwd_candidate_best.score < fwd_candidate.score:
                                fwd_candidate_best = fwd_candidate
                if fwd_candidate_best is not None:
                    out_bins[i].append((discriminatory_primer,
                                        wt_primer,
                                        fwd_candidate_best))
    return out_bins
# end def


if __name__ == '__main__':
    GENOME = os.path.join(LOCAL_DIR, 'test_data', '2014_02_18_gen9_54_failed_seg_fixed_with_gc_fixes_and_MLJ_fasta_edits.gb')
    REF_GENOME = os.path.join(LOCAL_DIR, 'test_data', 'mds42_full.gb')

    params = {
        'output_fp_base': 'seg23_masc_pcr_design'
    }
    doDesignPrimers(
                genome_gb=GENOME,
                ref_genome_gb=REF_GENOME,
                start_idx=1039358,        # Start index for seg 23
                end_idx=1084631,          # End index for seg 23
                border_feature_types=['gen9_fragment'],
                border_feature_regexs=None,
                cache_luts=True,
                params=params)
