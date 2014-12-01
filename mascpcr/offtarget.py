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
mascpcr.offtarget
~~~~~~~~~~~~~~~~~

Check a primer sequence for potential mispriming events against a genome
sequence. Methods are provided for both full hybridizatino-based checks and 
3-prime end stability checks. As of now are are using the full hybridization-
based checks in the pipeline.

"""
from __future__ import print_function

import multiprocessing as mp

import numpy as np
import primer3

from libnano import seqstr

import sys


def checkOffTarget(primer_str, genome_str, primer_idx, params,
                   hamming_percentile=0.05, genome_rc_str=None):
    """Return the tm, idx, and strand of the strongest off-target hybridization

    The 5' index of the primer on the fwd strand must be provided for masking
    purposes. The reverse complement of the genome, ``genome_rc_str``, may be
    provided as a performance optimization.

    Under the hood, :func:``checkOffTarget`` calculates the hamming distance 
    between the primer and the respective underlying sequence at each index
    of the genome and its reverse complement. As a means of optimization, 
    only the bottom ``hamming_percentile`` of hamming distance indices will
    also be screened with a thermodynamic alignment. For example, a 
    ``hamming_percentile`` of 0.05 will result thermodynamic alignments at
    the indices of the genome with hamming distances from the ``primer_str``
    in the bottom 0.05 percent.

    Args:
        primer_str (str)            : primer sequence string
        genome_str (str)            : genome sequence string
        primer_idx (int)            : 5'-most index of a primer sequence on the 
                                      forward strand, used to mask the binding 
                                      region
        params (dict)               : parameters dictionary used throughout the 
                                      pipeline (see 
                                      :module:``mascpcr.pipeline``)

        hamming_percentile (float)  : Hamming distance percentile (0-100) 
                                      below which regions surrounding the 
                                      respective indices will be subject to 
                                      interrogation by thermodynamic alignment
        genome_rc_str (str)         : reverse complement of the genome 
                                      (optimization to minimize the number of 
                                      times this operation must be performed / 
                                      number of memory copies)

    Returns:
        The tm (deg. C), index, and strand of the strongest off-target 
        hybridization.

    Raises:
        None

    """
    genome_rc_str = genome_rc_str or seqstr.reverseComplement(genome_str)
    primer_length = len(primer_str)

    strand_results = mp.Queue()

    def _fwdStrand():
        fwd_hamming_distances = seqstr.rollingHammingDistance(primer_str, 
                                                              genome_rc_str)
        fwd_hd_thresh = np.percentile(fwd_hamming_distances, hamming_percentile)
        fwd_primer_footprint = (-(primer_idx+primer_length), (-primer_idx))
        fwd_hamming_distances[fwd_primer_footprint[0]: \
                              fwd_primer_footprint[1]] = primer_length
        fwd_hotspots, = np.where((fwd_hamming_distances < fwd_hd_thresh))
        highest_tm_idx = None
        highest_tm = -100
        for idx in fwd_hotspots:
            tm = primer3.calcHeterodimerTm(
                primer_str, genome_str[-(idx+primer_length):-idx],
                **params['thermo_params'])
            if tm > highest_tm:
                highest_tm_idx = idx
                highest_tm = tm
        strand_results.put((highest_tm, highest_tm_idx, 1))

    def _revStrand():
        rev_hamming_distances = seqstr.rollingHammingDistance(primer_str, 
                                                              genome_str)
        rev_hd_thresh = np.percentile(rev_hamming_distances, hamming_percentile)
        rev_primer_footprint = ((primer_idx), (primer_idx+primer_length))
        rev_hamming_distances[rev_primer_footprint[0]: \
                              rev_primer_footprint[1]] = primer_length
        rev_hotspots, = np.where((rev_hamming_distances < rev_hd_thresh))

        highest_tm_idx = None
        highest_tm = -100
        for idx in rev_hotspots:
            tm = primer3.calcHeterodimerTm(
                primer_str, genome_rc_str[idx:idx+primer_length],
                **params['thermo_params'])
            if tm > highest_tm:
                highest_tm_idx = idx
                highest_tm = tm
        strand_results.put((highest_tm, highest_tm_idx, 0))

    fwd_proc = mp.Process(target=_fwdStrand)
    rev_proc = mp.Process(target=_revStrand)
    fwd_proc.start()
    rev_proc.start()
    fwd_proc.join()
    rev_proc.join()

    res1 = strand_results.get()
    res2 = strand_results.get()

    return max(res1[0], res2[0])


def checkOffTarget3p(primer_str, genome_str, primer_idx, params,
                     hamming_percentile=0.05, genome_rc_str=None):
    """Return the dG, idx, and strand of the strongest off-target 3p 
    hybridization event. 

    The 5' index of the primer on the fwd strand must be provided for masking
    purposes. The reverse complement of the genome, ``genome_rc_str``, may be
    provided as a performance optimization.

    Under the hood, :func:``checkOffTarget`` calculates the hamming distance 
    between the primer and the respective underlying sequence at each index
    of the genome and its reverse complement. As a means of optimization, 
    only the bottom ``hamming_percentile`` of hamming distance indices will
    also be screened with a thermodynamic alignment. For example, a 
    ``hamming_percentile`` of 0.05 will result thermodynamic alignments at
    the indices of the genome with hamming distances from the ``primer_str``
    in the bottom 0.05 percent.

    Args:
        primer_str (str)                : primer sequence string
        genome_str (str)                : genome sequence string
        primer_idx (int)                : 5'-most index of a primer sequence on 
                                          the forward strand, used to mask the 
                                          binding region
        params (dict)                   : parameters dictionary used throughout 
                                          the pipeline (see 
                                          :module:``mascpcr.pipeline``)

        hamming_percentile (float, optional)    : Hamming distance percentile 
                                         (0-100) below which regions 
                                         surrounding the respective indices 
                                         will be subject to interrogation by 
                                         thermodynamic alignment
        genome_rc_str (str, optional)   : reverse complement of the genome 
                                          (optimization to minimize the number 
                                          of times this operation must be 
                                          performed / number of memory copies)

    Returns:
        The tm (deg. C), index, and strand of the strongest off-target 
        hybridization.

    Raises:
        None

    """
    genome_rc_str = genome_rc_str or seqstr.reverseComplement(genome_str)
    primer_length = len(primer_str)

    strand_results = mp.Queue()

    def _fwdStrand():
        fwd_hamming_distances = seqstr.rollingHammingDistance(primer_str[12:], 
                                                              genome_rc_str)
        fwd_hd_thresh = np.percentile(fwd_hamming_distances, hamming_percentile)
        fwd_primer_footprint = (-(primer_idx+primer_length), (-primer_idx))
        fwd_hamming_distances[fwd_primer_footprint[0]: \
                              fwd_primer_footprint[1]] = primer_length
        fwd_hotspots, = np.where((fwd_hamming_distances < fwd_hd_thresh))
        highest_dg_idx = None
        highest_dg = -100

        for idx in fwd_hotspots:
            dg = primer3.calcEndStability(
                primer_str, genome_str[-(idx+primer_length):-idx],
                **params['thermo_params']).dg
            if dg > highest_dg:
                highest_dg_idx = idx
                highest_dg = dg
        
        sys.stdout.flush()

        strand_results.put((highest_dg, highest_dg_idx, 1))

    def _revStrand():
        rev_hamming_distances = seqstr.rollingHammingDistance(primer_str[12:], 
                                                              genome_str)
        rev_hd_thresh = np.percentile(rev_hamming_distances, hamming_percentile)
        rev_primer_footprint = ((primer_idx), (primer_idx+primer_length))
        rev_hamming_distances[rev_primer_footprint[0]: \
                              rev_primer_footprint[1]] = primer_length
        rev_hotspots, = np.where((rev_hamming_distances < rev_hd_thresh))

        highest_dg_idx = None
        highest_dg = -100

        for idx in rev_hotspots:
            dg = primer3.calcEndStability(
                primer_str, genome_rc_str[idx:idx+primer_length],
                **params['thermo_params']).dg
            if dg > highest_dg:
                highest_dg_idx = idx
                highest_dg = dg

        strand_results.put((highest_dg, highest_dg_idx, 0))

    fwd_proc = mp.Process(target=_fwdStrand)
    rev_proc = mp.Process(target=_revStrand)
    fwd_proc.start()
    rev_proc.start()
    fwd_proc.join()
    rev_proc.join()

    res1 = strand_results.get()
    res2 = strand_results.get()

    return max(res1[0], res2[0])
