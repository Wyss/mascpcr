'''
mascpcr | mascpcr.py

This is the main entry point for masc pcr design.

'''
from __future__ import print_function

import copy
import functools
import hashlib
import itertools
import multiprocessing as mp
import os
import operator
import sys
import six

import primer3
import bitarray
import numpy as np
from Bio import SeqIO

LOCAL_DIR = os.path.dirname(os.path.realpath(__file__))
CACHE_DIR = os.path.join(LOCAL_DIR, 'cache')

import buildindex
import primercandidate
import genbankfeatures
import offtarget
import indexutils

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
    # # [thermo] Max run of A bases allowed in a primer
    # 'max_a': 4,
    # # [thermo] Max run of T bases allowed in a primer
    # 'max_t': 4,
    # # [thermo] Max run of G bases allowed in a primer
    # 'max_g': 3,
    # # [thermo] Max run of C bases allowed in a primer
    # 'max_c': 3,
    # # [thermo] Max run of mixed A + T bases allowed in a primer
    # 'max_at': 4,
    # # [thermo] Max run of mixed G + C bases allowed in a primer
    # 'max_gc': 3,
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

def md5fn(string):
    raw_hash = hashlib.md5()
    if isinstance(string, six.string_types):
        raw_hash.update(string.encode('utf-8'))
    else:   # assume bytes (PY3) or str (PY2)
        raw_hash.update(string)
    raw_fn = raw_hash.hexdigest()
    # Remove forbidden / problematic characters
    raw_fn.replace('/', '_').replace('\\', '_').replace(':', '-')
    fn = 'cache' + raw_fn
    return fn


def genCacheFp(string, ext=''):
    return os.path.join(CACHE_DIR, md5fn(string)) + ext


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
    genome_seq = str(genome_rec.seq).encode('utf-8')
    ref_genome_seq = str(ref_genome_rec.seq).encode('utf-8')
    print('Building index with Mauve...')
    idx_lut = None
    if cache_luts:
        idx_lut_fp = genCacheFp(repr(locals) + genome_gb + ref_genome_gb,
                                '.npy')
        if os.path.isfile(idx_lut_fp):
            idx_lut = np.load(idx_lut_fp)
    if idx_lut is None:
        idx_lut = buildindex.buildIndex(genome_gb, ref_genome_gb,
                                        genome_seq, ref_genome_seq)
    if cache_luts:
        np.save(idx_lut_fp, idx_lut)
    border_arr = None
    if border_feature_types is not None:
        print('Finding border indices of feature type(s) {}...'.format(
              border_feature_types))
        border_arr = genbankfeatures.genBorderArray(genome_rec,
                                border_feature_types, border_feature_regexs)
        print('Found {} border indices'.format(border_arr.count(1)))
    return designPrimers(idx_lut, genome_seq, ref_genome_seq, start_idx,
                         end_idx, border_arr=border_arr,
                         cache_luts=cache_luts, params=params)



# Main entry point for MASC PCR primer design
def designPrimers(
            index_lut,          # LUT mapping genome index to ref_genome index
            genome,             # Genome sequence string
            ref_genome,         # Reference genome sequence string
            start_idx,          # Start index for MASC PCR design
            end_idx,            # End index (inclusive) for MASC PCR design
            edge_lut=None,      # Pre-populated edge LUT
            mismatch_lut=None,  # Pre-populated mismatch LUT
            border_arr=None,    # Array of feature border indices to encompass
            cache_luts=False,   # Cache the LUTs w/ an md5 hash of the genome
            params=None):

    # ~~~~~~~~~~ Reconcile default params with user-provided params ~~~~~~~~~ #
    _params = copy.deepcopy(DEFAULT_PARAMS)
    _params.update(params or {})
    params = _params

    # ~~~~~~~~~~~~~~~~~~~~~~~~~ Build LUTs if needed ~~~~~~~~~~~~~~~~~~~~~~~~ #
    if cache_luts:
        cache_dir = os.path.join(LOCAL_DIR, 'cache')
        try:
            os.makedirs(cache_dir)
        except OSError as e:
            if e.errno != 17:
                raise
        genome_hash = md5fn(genome)

        print('Building edge lookup table...')
        cache_fn = genome_hash + '_edgelut'
        full_cache_fp = os.path.join(cache_dir, cache_fn + '.bitarray')
        if os.path.isfile(full_cache_fp) and edge_lut is None:
            edge_lut = bitarray.bitarray()
            with open(full_cache_fp, 'rb') as fd:
                edge_lut.fromfile(fd)
        else:
            edge_lut = edge_lut or indexutils.findEdges(index_lut)
            with open(full_cache_fp, 'wb') as fd:
                edge_lut.tofile(fd)

        print('Building mismatch lookup table...')
        cache_fn = genome_hash + '_mismatchlut'
        full_cache_fp = os.path.join(cache_dir, cache_fn + '.bitarray')
        if os.path.isfile(full_cache_fp) and mismatch_lut is None:
            mismatch_lut = bitarray.bitarray()
            with open(full_cache_fp, 'rb') as fd:
                mismatch_lut.fromfile(fd)
        else:
            mismatch_lut = mismatch_lut or (mismatch_lut or
                    indexutils.findMismatches(index_lut, genome, ref_genome))
            with open(full_cache_fp, 'wb') as fd:
                mismatch_lut.tofile(fd)

    else:
        print('Building edge lookup table...')
        edge_lut = edge_lut or indexutils.findEdges(index_lut)
        print('Building mismatch lookup table...')
        mismatch_lut = (mismatch_lut or
                        indexutils.findMismatches(index_lut, genome, ref_genome))

    # ~~~~~~~ Generate potential primers for both fwd and rev strands ~~~~~~~ #
    fwd_strand_primer_candidates = []
    rev_strand_primer_candidates = []
    fwd_candidates_found = 0
    rev_candidates_found = 0
    candidates_checked = 0

    mismatch_count = mismatch_lut.count(1)
    print("Total Mismatches: {}, Percent Mismatches: {:.2f}".format(mismatch_count,
          mismatch_count/float(len(mismatch_lut))*100))
    print('Finding candidate primers...')
    for idx, is_mismatch in enumerate(mismatch_lut[start_idx:end_idx+1],
                                      start=start_idx):
        if is_mismatch:
            fc = primercandidate.findCandidateSequence(idx, index_lut, genome,
                                                       ref_genome, edge_lut,
                                                       mismatch_lut, params,
                                                       check_wt_primer=True)
            rc = primercandidate.findCandidateSequence(idx, index_lut, genome,
                                                       ref_genome, edge_lut,
                                                       mismatch_lut, params,
                                                       fwd_strand=False,
                                                       check_wt_primer=True)
            if fc is not None:
                fwd_strand_primer_candidates.append(fc)
                fwd_candidates_found += 1
            if rc is not None:
                rev_strand_primer_candidates.append(rc)
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
                            index_lut,
                            genome, ref_genome,
                            edge_lut, mismatch_lut,
                            params=DEFAULT_PARAMS,
                            is_long_to_short=True)
    print("fwd bin counts (l2s)")
    for i in range(len(fwd_bins_l2s)):
        print("bin %d: %d" % (i, len(fwd_bins_l2s[i])))


    rev_bins_l2s = partitionCandidatesREV(rev_strand_primer_candidates,
                            start_idx, end_idx,
                            index_lut,
                            genome, ref_genome,
                            edge_lut, mismatch_lut,
                            params=DEFAULT_PARAMS,
                            is_long_to_short=True)
    print("\nrev bin counts (l2s)")
    for i in range(len(rev_bins_l2s)):
        print("bin %d: %d" % (i, len(rev_bins_l2s[i])))

    # ~~~~~~~~~~~~~~~~~~~~ Combine and score primer pairs ~~~~~~~~~~~~~~~~~~~ #
    print('Combining bins and scoring primer pairs...')
    combined_bins = [bf+br for bf, br in zip(fwd_bins_l2s, rev_bins_l2s)]
    for bin_idx, bin in enumerate(combined_bins):
        for pair_idx, pair in enumerate(bin):
            f_idx = pair[0].idx
            r_idx = pair[1].idx
            score = pair[0].score + pair[1].score
            if border_arr:
                # Add an arbitrarily high value to the score for each border
                # in the primer pair footprint
                score += (border_arr[f_idx:r_idx].count(1) * 50)
            combined_bins[bin_idx][pair_idx] = (score, pair[0], pair[1])
        # Sort each bin after scoring (highest score first)
        combined_bins[bin_idx].sort(key=lambda rec: -rec[0])
        best_pair = combined_bins[bin_idx][0]
        fwd_primer_seq = genome[best_pair[1].idx - best_pair[1].length:
                                best_pair[1].idx]
        fwd_primer_start_idx = best_pair[1].idx - best_pair[1].length
        # print('Bin {} | Best primer pair score: {}'.format(bin_idx,
        #                                                    best_pair[0]))
        # print('      | Fwd off-target tm: {}'.format(
        #         offtarget.checkOffTarget(fwd_primer_seq, genome,
        #                                  fwd_primer_start_idx)[0]))

    # ~~~~~~~~~~~~~~~ Find the best set of compatible primers ~~~~~~~~~~~~~~~ #
    # A simple BFS should find the best possible set of primer sets given our
    # constraints of:
    #   1) No heterodimerization within a set
    #   2) No off-target binding within the genome
    #
    print('Finding optimal set of MASC primers...')


    def getPrimerSeq(primer_candidate):
        pc_seq = genome[primer_candidate.idx:primer_candidate.idx +
                          primer_candidate.length]
        if primer_candidate.strand == 1:
            return pc_seq
        else:
            return seqstr.reverseComplement(pc_seq)

    #             if sum([bin_depth in up for bin_depth, up in zip(set_of_sets,
    #                     unacceptable_primers)]) > 0:

    def checkSetForHeterodimers(set_of_primer_sets, tm_max=40):
        for bin_idx, primer_pair_idx in enumerate(set_of_primer_sets):
            set1_primer_tuple = combined_bins[bin_idx][primer_pair_idx]
            set1_f_seq = getPrimerSeq(set1_primer_tuple[1])
            set1_r_seq = getPrimerSeq(set1_primer_tuple[2])
            for bin2_idx, primer_pair_idx2 in enumerate(set_of_primer_sets):
                set2_primer_tuple = combined_bins[bin2_idx][primer_pair_idx2]
                set2_f_seq = getPrimerSeq(set2_primer_tuple[1])
                set2_r_seq = getPrimerSeq(set2_primer_tuple[2])
                max_tm = max([
                    primer3.calcHeterodimer(p1, p2,
                        **params['thermo_params']).tm
                    for p1, p2 in itertools.product([set1_f_seq, set1_r_seq],
                                                    [set2_f_seq, set2_r_seq])
                ])
                if max_tm > tm_max:
                    print(set1_f_seq, set1_r_seq, set2_f_seq, set2_r_seq)
                    return (False, bin_idx, bin2_idx)
        return (True,)


    def checkSetForOffTarget(set_of_primer_sets, tm_max=40):
        for bin_idx, primer_pair_idx in enumerate(set_of_primer_sets):
            primer_tuple = combined_bins[bin_idx][primer_pair_idx]
            f_seq = getPrimerSeq(primer_tuple[1])
            f_idx = primer_tuple[1].idx
            r_seq = getPrimerSeq(primer_tuple[2])
            r_idx = primer_tuple[2].idx
            cached_tms = offtarget_tm_cache[bin_idx].get(primer_pair_idx)
            if cached_tms is not None:
                f_tm, r_tm = cached_tms
            else:
                f_res = offtarget.checkOffTarget(f_seq, genome, f_idx)
                r_res = offtarget.checkOffTarget(r_seq, genome, r_idx)
                f_tm, r_tm = f_res[0], r_res[0]
                offtarget_tm_cache[bin_idx][primer_pair_idx] = (f_tm, r_tm)
            if max(f_tm, r_tm) > tm_max:
                unacceptable_primers[bin_idx].append(primer_pair_idx)
                return False
        return True


    # def calculateSetScore(set_of_primer_sets):
    #     set_of_sets_score = 0
    #     for bin_idx, primer_pair_idx in enumerate(set_of_primer_sets):
    #         set_of_sets_score += combined_bins[bin_idx][primer_pair_idx][0]
    #     return set_of_sets_score

    num_bins = len(combined_bins)
    bin_lengths = [len(bin) for bin in combined_bins]
    unacceptable_primers = [[] for _ in range(num_bins)]
    heterodimer_conflicts = [[[] for _ in range(bin_lengths[bin_idx])]
                             for bin_idx in range(num_bins)]
    offtarget_tm_cache = [{} for _ in range(num_bins)]

    set_of_primer_sets = [0] * num_bins
    combinations_tried = 0
    num_possible_combinations = functools.reduce(operator.mul,
                                                 bin_lengths, 1)

    # for ps in combined_bins[0]:
    #     print(ps[1].strand)
    #     print(ps[1].seq)
    #     for midx in index_lut[ps[1].idx:ps[1].idx+ps[1].length]:
    #         print('*' if midx == -1 else ' ', end='')
    #     print('')
    #     print(index_lut[ps[1].idx:ps[1].idx+ps[1].length])
    #     print(ps[1].mismatch_idxs)
    #     print(getPrimerSeq(ps[1]))
    #     print(ps[2].strand)
    #     print(ps[2].seq)
    #     print(getPrimerSeq(ps[2]))
    #     print(genome[ps[2].idx-29:ps[2].idx+30])
    #     print(ps[1].idx-ps[2].idx)
    #     print()

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
            print('Heterodimer conflict between {}.{}, {}.{}'.format(bin1_idx, ps1_idx, bin2_idx, ps2_idx))
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

    def generateWtPrimer(discriminatory_primer):
        dp = discriminatory_primer
        wt_seq = ref_genome[index_lut[dp.idx]:index_lut[dp.idx + dp.length]]
        if dp.strand == -1:
            wt_seq = seqstr.reverseComplement(wt_seq)
        return primercandidate.CandidatePrimer(
            idx=index_lut[dp.idx],
            seq=wt_seq,
            strand=dp.strand,
            length=dp.length,
            mismatch_idxs=[0] * dp.length,
            tm=primer3.calcTm(wt_seq, **params['thermo_params']),
            tm_homo=primer3.calcHomodimer(wt_seq, **params['thermo_params']).tm,
            tm_hairpin=primer3.calcHairpin(wt_seq, **params['thermo_params']).tm,
            score=0
        )

    def writePrimerRec(fd, primer_obj, recoded_idx=True, gen_wt_idx=False):
        po = primer_obj
        if recoded_idx:
            if gen_wt_idx:
                po_idx = '{},{}'.format(po.idx, index_lut[po.idx])
            else:
                po_idx = '{},N/A'.format(po.idx, index_lut[po.idx])
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
            fp = primer_set[1]
            rp = primer_set[2]
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
            fp = primer_set[1]
            rp = primer_set[2]
            wt = generateWtPrimer(fp)
            fd_out.write('{},Discriminatory primer,'.format(bin_idx+1))
            writePrimerRec(fd_out, fp)
            fd_out.write(',Wildtype primer,')
            writePrimerRec(fd_out, wt, recoded_idx=False)
            fd_out.write(',Common primer,')
            writePrimerRec(fd_out, rp, gen_wt_idx=True)
        fd_out.write('\nDESIGN PARAMETERS' + '\n' + '\n'.join(
            ['{},{}'.format(k, repr(v).replace(',', ';')) for k, v in
             params.items()]))


def printPrimerPair(primerpair, idx_lut, genome, ref_genome):
    idx0 = primerpair[0].idx
    l0 = primerpair[0].length
    idx1 = primerpair[1].idx
    l1 = primerpair[1].length
    print(genome[idx0-l0:idx0+1], genome[idx1-l1:idx1+1])
    idxmap0 = idx_lut[idx0]
    idxmap1 = idx_lut[idx1]
    print(ref_genome[idxmap0-l0:idxmap0+1], ref_genome[idxmap1-l1:idxmap1+1])


def partitionCandidatesFWD(discriminatory_primers,
                        start_idx, end_idx,
                        index_lut,
                        genome, ref_genome,
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
    for discriminatory_primer in discriminatory_primers:
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
                rev_candidate = primercandidate.findCandidateSequence(
                    j, index_lut, genome, ref_genome, edge_lut, mismatch_lut,
                    params, fwd_strand=False,
                    disallow_mismatches=True)
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
                out_bins[i].append((discriminatory_primer, rev_candidate_best))
    return out_bins
# end def

def partitionCandidatesREV(discriminatory_primers,
                           start_idx, end_idx,
                           index_lut,
                           genome, ref_genome,
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
    for discriminatory_primer in discriminatory_primers:
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
                    fwd_candidate = primercandidate.findCandidateSequence(
                        j, index_lut, genome, ref_genome, edge_lut,
                        mismatch_lut, params, fwd_strand=True,
                        disallow_mismatches=True)
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
