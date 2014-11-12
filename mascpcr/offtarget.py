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

import multiprocessing as mp

import numpy as np
import primer3

from libnano import seqstr


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
        ``primer_str``: primer sequence string
        ``genome_str``: genome sequence string
        ``primer_idx``: 5'-most index of a primer sequence on the forward 
                        strand, used to mask the binding region
        ``params``: parameters dictionary used throughout the pipeline (see
                    :module:``mascpcr.pipeline``)

    Kwargs:
        ``hamming_percentile``: Hamming distance percentile (0-100, float) 
                                below which regions surrounding the respective 
                                indices will be subject to interrogation by
                                thermodynamic alignment 
        ``genome_rc_str``: reverse complement of the genome (optimization to
                           minimize the number of times this operation must
                           be performed / number of memory copies)

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
        rev_primer_footprint = ((primer_idx)-primer_length,
                                (primer_idx+primer_length))
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


def checkOffTarget3p(primer_str, genome_str, primer_idx,
                     num_bases=10, thresholds=None, genome_rc_str=None):
    """Check for primer mispriming using 3 primer end homology

    This function is not production-ready, but represents an alternative
    to :func:``checkOffTarget`` that uses a lighter-weight heurstic that is
    solely based on 3-prime end homology between the primer and the genome. 
    See the ``libnano`` documentation for ``libnano.seqstr.can3pMisprime`` for
    more information regarding the 3 prime mispriming algorithm. 

    Args:
        ``primer_str``: primer sequence string
        ``genome_str``: genome sequence string
        ``primer_idx``: 5-prime most index of the primer binding site

    Kwargs:
        ``num_bases``: number of bases to consider at the 3 prime end of the 
                       primer sequence
        ``thresholds``: list of mismatch thresholds, keyed by hamming distance;
                        see ``libnano.seqstr.can3pMisprime`` documentation
        ``genome_rc_str``: reverse complement of the genome (optimization to
                           minimize the number of times this operation must
                           be performed / number of memory copies)

    Returns:
        ``True`` if the 3 prime end can misprime on either strand of the 
        genome, else ``False``

    Raises:
        ``ValueError``

    """
    genome_rc_str = genome_rc_str or seqstr.reverseComplement(genome_str)
    primer_length = len(primer_str)
    subprimer = primer_str[-num_bases:]

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
            raise ValueError("``num_bases`` too low for default thresholds")
    elif len(thresholds) != num_bases:
        raise ValueError('``thresholds`` must be the same length as '
                         '``num_bases``')

    strand_results = mp.Queue()

    def _fwdStrand():
        a = 1
        if num_bases < primer_idx:
            a = seqstr.can3pMisprime(subprimer,
                                        genome_rc_str[:primer_idx],
                                        thresholds,
                                        1)
        b = seqstr.can3pMisprime(subprimer,
                                    genome_rc_str[primer_idx+primer_length:],
                                    thresholds,
                                    1)
        strand_results.put(a and b)

    def _revStrand():
        a = 1
        if num_bases < primer_idx:
            a = seqstr.can3pMisprime(subprimer,
                                        genome_str[:primer_idx],
                                        thresholds,
                                        1)
        b = seqstr.can3pMisprime(subprimer,
                                    genome_str[primer_idx+primer_length:],
                                    thresholds,
                                    1)
        strand_results.put(a and b)

    fwd_proc = mp.Process(target=_fwdStrand)
    rev_proc = mp.Process(target=_revStrand)
    fwd_proc.start()
    rev_proc.start()
    fwd_proc.join()
    rev_proc.join()

    res1 = strand_results.get()
    res2 = strand_results.get()

    return res1 and res2
