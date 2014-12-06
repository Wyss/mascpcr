# Copyright (C) 2014. Ben Pruitt & Nick Conway
# See LICENSE for full GPLv2 license.
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
"""
mascpcr.pipeline
~~~~~~~~~~~~~~~~

This is the main entry point for MASC PCR primer design. It contains the core
`findMascPrimers` method as well as some helper methods related to LUT
generation (if you're generating multiple MASC PCR primer sets from the same
genome/reference genome it's best to just make these calls once).

"""
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
# ~~~ #

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
    # minimum base distance of a primer to the boundaries of a bin
    # (not implemented)
    'bin_edge_offset': 67,
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
    'output_xlsx': True,
    # Format of the .xlsx file (``continuous`` or ``by_row``)
    # If the format is "continuous", the primers will be loaded in row order
    # with no gaps (A1: bin 1 d_primer, A2: bin 2 d_primer, etc)
    # If the format is "by_row", primers will be separated in rows by type.
    # If possible, the rows of primers will be separated by empty rows
    'xlsx_format': 'by_row'
}

MascPrimerSet = namedtuple('MascPrimerSet',
        ['score',           # Sum of the d_primer and c_primer scores
         'd_primer',        # Discriminatory primer
         'w_primer',        # Wild-type primer
         'c_primer'         # Common primer
         ])


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Helper functions ~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
def generateLUTs(genome_fp, ref_genome_fp, start_idx, end_idx,
                 border_feature_types=None, border_qualifier_regexs=None,
                 cache_luts=True, cache_dir=CACHE_DIR):
    """Generate the LUT datastructures necessary for the mascpcr pipeline

    Builds the index, edge, mismatch, and (optionally) feature border lookup
    tables. All LUTs are of equal length to the recoded / modified genome.

        *index (idx_lut)*: LUT mapping recoded / modified genome indices to
                           wildtype genome indices. The numpy array is indexed
                           by recoded genome index and the value at each
                           position is the corresponding wildtype genome index.
                           If the value at a given position is -1, then that
                           base in the recoded genome did not map to the
                           wildtype genome.

        *edge (edge_lut)*: LUT of 'edges' in the index mapping. Edges are
                           discontinuities in the mapping where a continuous
                           region of recoded -> wildtype mapping is
                           interrupted.

        *mismatch (mismatch_lut)*: LUT of localized base mismatches used to
                                   identify good points of discrimination
                                   between the recoded and wildtype genomes.

        *border (border_lut)*: LUT of the border indices of features matching
                               user-profided type and qualifier specifications.
                               Used to prioritize primer sets that span such
                               borders as a means of validating the success
                               of an assembly.

    Args:
        genome_fp (str)                         : filepath to recoded /
                                                  modified genome genbank file
        ref_genome_fp (str)                     : filepath to wt / reference
                                                  genome genbank file
        start_idx (int)                         : start index of recoded genome
                                                  for primer design
        end_idx (int)                           : end index of recoded genome
                                                  for primer design

        border_feature_types (list, optional)   : list of feature types
                                                  (e.g., ['gene', 'CDS'])
        border_qualifier_regexs (dict, optional): dict of <field name>:
                                                  <value regex> entries
        cache_luts (bool, optional)             : whether or not to cache the
                                                  LUTs on disk for future runs
        cache_dir (str, optional)               : the directory in which to
                                                  cache the LUTs (defaults to
                                                  the ``cache/`` directory
                                                  inside of the module
                                                  directory)

    Returns:
        Genome sequence (str), reference genome sequence (str), index_lut
        (numpy array), edge_lut (binary bitarray), mismatch_lut
        (binary bitarray), (optionally) border_lut (binary bitarray) -- None if
        not generated.

    Raises:
        ``OSError``

    """

    print('Reading in genome and reference genome genbank files...')
    genome_rec = SeqIO.read(genome_fp, 'gb')
    ref_genome_rec = SeqIO.read(ref_genome_fp, 'gb')
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
    idx_lut = indexing.buildIdxLUT(genome_fp, ref_genome_fp, genome_str,
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
                                                    border_qualifier_regexs)
        print('Found {} border indices'.format(border_lut.count(1)))

    return (genome_str, ref_genome_str, idx_lut, edge_lut, mismatch_lut,
            border_lut)


# ~~~~~~~~~~~~~~~~~~~~~~~~ Main pipeline entry point ~~~~~~~~~~~~~~~~~~~~~~~~ #

def findMascPrimers(idx_lut, genome_str, ref_genome_str, start_idx, end_idx,
                    edge_lut, mismatch_lut, border_lut=None, params=None):
    """Find a set of MASC PCR primer sets for the given genomic region

    This is the main entry point for the ``mascpcr`` design pipeline. It
    requires the the underlying LUTs have been prebuilt, and performs all
    primer finding and output generating tasks as specified in the user-
    provided parameters (``params``).

    Args:
        idx_lut (``numpy.ndarray``)         : numpy array of index mappings
                                              between two genomes (see
                                              :func:``buildIdxLUT``)
        genome_str (str)                    : recoded/modified genome sequence
        ref_genome_str (str)                : reference genome sequence
        start_idx (int)                     : start index of recoded genome for
                                              primer design
        end_idx (int)                       : end index of recoded genome for
                                              primer design
        edge_lut (``bitarray.bitarray``)    : bitarray LUT of discontinuities
                                              in the index mapping
        mismatch_lut (``bitarray.bitarray``): bitarray LUT of localized
                                              mismatches in the index mapping

        border_lut (``bitarray.bitarray``, optional): bitarray LUT of
                                              user-specified feature border
                                              indices
        params (dict, optional)             : user-specified pipeline
                                              parameters to override the
                                              defaults (see ``DEFAULT_PARAMS``)

    Returns:
        ``None`` (prints status updates to stdout and writes output to files
        as specified in the :param:``params``)

    Raises:
        ``OSError``

    """

    # ~~~~~~~~~~ Reconcile default params with user-provided params ~~~~~~~~~ #
    _params = copy.deepcopy(DEFAULT_PARAMS)
    _params.update(params or {})
    params = _params
    thermo_params = params.get('thermo_params')

    # Coerce off-by-one index issues
    if end_idx == len(genome_str):
        end_idx = len(genome_str) - 1

    # ~~~~~~~ Generate potential primers for both fwd and rev strands ~~~~~~~ #
    fwd_strand_primer_candidates = []
    rev_strand_primer_candidates = []
    fwd_candidates_found = 0
    rev_candidates_found = 0
    candidates_checked = 0
    print('Finding MASC PCR primers for target with output basename: ',
          params['output_basename'])
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

    # Get target indices / seq area for local mispriming checks
    target_area_offset = min(1000, (len(genome_str) -
                            (end_idx - start_idx) / 2))
    target_start_idx = start_idx - target_area_offset
    target_end_idx = end_idx + target_area_offset

    target_area = genome_str[max(0, target_start_idx):
                             min(target_end_idx, len(genome_str))]

    if target_start_idx < 0:
        target_area = genome_str[target_start_idx:] + target_area
    if target_end_idx > len(genome_str):
        target_area += genome_str[:target_end_idx - len(genome_str)]

    ref_target_start_idx = idx_lut[start_idx]
    ref_target_end_idx = idx_lut[end_idx]
    ref_start_offset = 0
    ref_end_offset = 0

    if ref_target_start_idx == -1:
        ref_start_offset = 1
        while (ref_target_start_idx == -1 and start_idx -
                ref_start_offset > 0):
            ref_target_start_idx = idx_lut[start_idx - ref_start_offset]
            ref_start_offset += 1
    if ref_target_end_idx == -1:
        ref_end_offset = 1
        while (ref_target_end_idx == -1 and
               target_end_idx + ref_end_offset < len(genome_str)):
            ref_target_end_idx = idx_lut[end_idx + ref_end_offset]
            ref_end_offset += 1

    ref_target_start_idx -= target_area_offset
    ref_target_end_idx += target_area_offset

    ref_target_area = ref_genome_str[max(0, ref_target_start_idx):
                                     min(ref_target_end_idx,
                                         len(ref_genome_str))]

    if ref_target_start_idx < 0:
        ref_target_area = ref_genome_str[ref_target_start_idx:] + \
                          ref_target_area
    if ref_target_end_idx > len(ref_genome_str):
        ref_target_area += ref_genome_str[:ref_target_end_idx - \
                           len(ref_genome_str)]


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
                d_tm = offtarget.checkOffTarget(d_primer.seq, target_area,
                                                (d_primer.idx - start_idx +
                                                target_area_offset),
                                                params)
                w_tm = offtarget.checkOffTarget(w_primer.seq, ref_target_area,
                                                (w_primer.idx -
                                                idx_lut[start_idx] +
                                                target_area_offset -
                                                ref_start_offset), params)
                c_tm = offtarget.checkOffTarget(c_primer.seq, target_area,
                                                (c_primer.idx - start_idx +
                                                target_area_offset), params)
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
                    sys.exit(1)
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
                bin1_delta = (combined_bins[bin1_idx][ps1_idx][0] -
                              combined_bins[bin1_idx][ps1_idx + 1][0])
            except IndexError:
                bin1_delta = -1
            try:
                bin2_delta = (combined_bins[bin2_idx][ps2_idx][0] -
                              combined_bins[bin2_idx][ps2_idx + 1][0])
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
                    if bin1_idx != bin2_idx:
                        set_of_primer_sets[bin1_idx] = 0
            elif bin2_delta == -1:
                set_of_primer_sets[bin1_idx] += 1
                if bin1_idx != bin2_idx:
                    set_of_primer_sets[bin2_idx] = 0
            elif bin1_delta < bin2_delta:
                set_of_primer_sets[bin1_idx] += 1
                if bin1_idx != bin2_idx:
                    set_of_primer_sets[bin2_idx] = 0
            else:
                set_of_primer_sets[bin2_idx] += 1
                if bin1_idx != bin2_idx:
                    set_of_primer_sets[bin1_idx] = 0
            heterodimer_conflicts[bin1_idx][ps1_idx].append((bin2_idx, ps2_idx))
            heterodimer_conflicts[bin2_idx][ps1_idx].append((bin1_idx, ps1_idx))
            continue
        if not checkSetForOffTarget(set_of_primer_sets,
                                    params['spurious_tm_clip']):
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

