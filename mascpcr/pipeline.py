'''
pipeline | pipeline.py
~~~~~~~~~~~~~~~~~~~~~~

This is the main entry point for MASC PCR primer design. It contains the core
`findMascPrimers` method as well as some helper methods related to LUT 
generation (if you're generating multiple MASC PCR primer sets from the same
genome/reference genome it's best to just make these calls once).

'''
from __future__ import print_function

import copy
import itertools
import os
import sys

from collections import namedtuple

import primer3
from Bio import SeqIO

LOCAL_DIR = os.path.dirname(os.path.realpath(__file__))
CACHE_DIR = os.path.join(LOCAL_DIR, 'cache')
CWD = os.getcwd()

import indexing
import ioutil
import primercandidate
import primerpartition
import genbankfeatures
import offtarget


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
    # Tm above which primer candidates will be thrown out for 
    # homo/heterodimerization, hairpin formation, mispriming, etc.
    'spurious_tm_clip': 40,
    # Filepath at which to write output files, defaults to current working dir
    'output_fp': CWD,
    # Default output basename for files and sequences
    'output_basename': 'mascpcr_out',
    # Minimum number of mismatches for discriminatory primer
    'min_num_mismatches': 2,
    # Whether or not to output a pipeline report
    'output_report': True,
    # Whether or not to dump these parameters into the output report
    'dump_params': True,
    # Whether or not to output an .xlsx file for IDT orders
    'output_xlsx': True
}

MascPrimerSet = namedtuple('MascPrimerSet',
        ['score',           # Sum of the d_primer and c_primer scores
         'd_primer',        # Discriminatory primer
         'w_primer',        # Wild-type primer
         'c_primer'         # Common primer
         ])


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Helper functions ~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
def generateLUTs(
            genome_gb,          # Genome sequence genbank file
            ref_genome_gb,      # Reference genome genbank file
            start_idx,          # Start index for MASC PCR design
            end_idx,            # End index (inclusive) for MASC PCR design
            border_feature_types=None, # List of genbank feature types
            border_feature_regexs=None, # List of qualifier: regex keys
            cache_luts=True,    # Whether or not to cache the LUTs
            cache_dir=CACHE_DIR):   # Directory in which to cache LUTs

    print('Reading in genome and reference genome genbank files...')
    genome_rec = SeqIO.read(genome_gb, 'gb')
    ref_genome_rec = SeqIO.read(ref_genome_gb, 'gb')
    genome_str = str(genome_rec.seq).encode('utf-8')
    ref_genome_str = str(ref_genome_rec.seq).encode('utf-8')

    # ~~~~~~~~~~~~~~~~~~~~~~~~~ Build LUTs as needed ~~~~~~~~~~~~~~~~~~~~~~~~ #

    if cache_luts:
        try:
            os.makedirs(cache_dir)
        except OSError as e:
            if e.errno != 17:
                raise
    else:
        cache_dir = None

    print('Building index lookup table with Mauve Aligner...')
    idx_lut = indexing.buildIdxLUT(genome_gb, ref_genome_gb, genome_str,
                                   ref_genome_str, cache_dir=cache_dir)

    print('Building edge lookup table...')
    edge_lut = indexing.buildEdgeLUT(idx_lut, cache_dir=cache_dir)

    print('Building mismatch lookup table...')
    mismatch_lut = indexing.buildMismatchLUT(idx_lut, genome_str,
                                             ref_genome_str,
                                             cache_dir=cache_dir)

    border_lut = None
    if border_feature_types is not None:
        print('Finding border indices of feature type(s): {}...'.format(
              ','.join(border_feature_types)))
        border_lut = genbankfeatures.buildBorderLUT(genome_rec,
                                                    border_feature_types,
                                                    border_feature_regexs)
        print('Found {} border indices'.format(border_lut.count(1)))

    return idx_lut, edge_lut, mismatch_lut, border_lut


# ~~~~~~~~~~~~~~~~~~~~~~~~ Main pipeline entry point ~~~~~~~~~~~~~~~~~~~~~~~~ #

# Main entry point for MASC PCR primer design
def findMascPrimers(
            idx_lut,            # LUT mapping genome index to ref_genome index
            genome_str,         # Genome sequence string
            ref_genome_str,     # Reference genome sequence string
            start_idx,          # Start index for MASC PCR design
            end_idx,            # End index (inclusive) for MASC PCR design
            edge_lut,           # Pre-populated edge LUT
            mismatch_lut,       # Pre-populated mismatch LUT
            border_lut=None,    # Array of feature border indices to overlap
            params=None):       # General design params (see DEFAULT_PARAMS)

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
    fwd_bins_l2s = primerpartition.partitionCandidatesFWD(
                            fwd_strand_primer_candidates, start_idx, end_idx,
                            idx_lut, genome_str, ref_genome_str, edge_lut,
                            mismatch_lut, params=DEFAULT_PARAMS)
    print("Fwd bin counts (l2s):")
    for i in range(len(fwd_bins_l2s)):
        print("\tBin %d: %d" % (i, len(fwd_bins_l2s[i])))


    rev_bins_l2s = primerpartition.partitionCandidatesREV(
                            rev_strand_primer_candidates, start_idx, end_idx,
                            idx_lut, genome_str, ref_genome_str, edge_lut,
                            mismatch_lut, params=DEFAULT_PARAMS)
    print("\nRev bin counts (l2s)")
    for i in range(len(rev_bins_l2s)):
        print("\tBin %d: %d" % (i, len(rev_bins_l2s[i])))

    # ~~~~~~~~~~~~~~~~~~~~ Combine and score primer pairs ~~~~~~~~~~~~~~~~~~~ #
    print('Combining bins and scoring primer pairs...')
    combined_bins = [bf + br for bf, br in zip(fwd_bins_l2s, rev_bins_l2s)]
    for bin_idx, bin in enumerate(combined_bins):
        for pair_idx, pair in enumerate(bin):
            f_idx = pair[0].idx
            r_idx = pair[2].idx
            idx1, idx2 = (f_idx, r_idx) if f_idx < r_idx else (r_idx, f_idx)
            score = pair[0].score + pair[2].score
            if border_lut:
                # Add an arbitrarily high value to the score for each border
                # in the primer pair footprint
                score += (border_lut[idx1:idx2].count(1) * 50)
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

        print('Checking binned primer sets: {}'.format(set_of_primer_sets),
              end='')
        hetero_res = checkSetForHeterodimers(set_of_primer_sets, 40)
        if not hetero_res[0]:
            bin1_idx, bin2_idx = hetero_res[1], hetero_res[2]
            ps1_idx = set_of_primer_sets[bin1_idx]
            ps2_idx = set_of_primer_sets[bin2_idx]
            print('->Heterodimer conflict between {}.{}, {}.{}'.format(bin1_idx,
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
                    print('\nFailed to find MASC PCR primer set due to total '
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

    print('\nFound set of primer sets: ', set_of_primer_sets)

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~ Write output files ~~~~~~~~~~~~~~~~~~~~~~~~~ #
    if params['output_report']:
        ioutil.writeCsvReport(idx_lut, set_of_primer_sets, combined_bins, 
                              params)
    if params['output_xlsx']:
        ioutil.writeIdtXslxFile(idx_lut, set_of_primer_sets, combined_bins, 
                                params)

