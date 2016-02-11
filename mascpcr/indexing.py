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
mascpcr.indexing
~~~~~~~~~~~~~~~~

Methods related to generating index lookup tables (LUTs) between two related
sequences (e.g., a recoded and reference genome).

Also includes methods for generating arrays of indexing "edges"
(incongruencies in an index LUT, i.e., borders of insertions, deletions,
rearrangements, etc.) and single base (non-inseration/deletion/rearrangement)
discrepencies.

"""
import hashlib
import os
import six

import numpy as np
import mauve


def sha256fn(string):
    """Generate a filename-friendly sha256 hash of the provided string."""
    raw_hash = hashlib.sha256()
    if isinstance(string, six.string_types):
        raw_hash.update(string.encode('utf-8'))
    else:   # assume bytes (PY3) or str (PY2)
        raw_hash.update(string)
    fn = raw_hash.hexdigest()
    # Remove forbidden / problematic characters for filenames
    fn.replace('/', '_').replace('\\', '_').replace(':', '-')
    return fn


def getCachedNumpyArray(cache_dir, hash_str):
    """Load a cached numpy array using a sha256 hash of the provided string.

    Uses :func:``sha256fn`` to generate a filename-friendly sha256 hash string.

    Args:
        cache_dir (str)     : string filepath to the cache directory
        hash_str (str)      : string to run through :func:``sha256fn`` to
                              generate the filename hash

    Returns:
        Numpy array object if the cached array is successfully loaded,
        otherwise, ``None``.

    Raises:
        ``OsError``

    """
    cached_fn = sha256fn(hash_str) + '.npy'
    cached_fp = os.path.join(cache_dir, cached_fn)
    np_arr = None
    if os.path.isfile(cached_fp):
        try:
            np_arr = np.load(cached_fp)
        except OSError:
            pass
    return np_arr


def saveNumpyArrayToCache(np_arr, cache_dir, hash_str):
    """Save a numpy array to disk using a sha256 hash of the provided string.

    Uses :func:``sha256fn`` to generate a filename-friendly sha256 hash string.
    Expects that ``cache_dir`` already exists, raises ``OSError`` if it does
    not.

    Args:
        np_arr (``numpy.ndarray``)  : numpy array to save to disk
        cache_dir (str)             : string filepath to the cache directory
        hash_str (str)              : string to run through :func:``sha256fn``
                                      to generate the filename hash

    Returns:
        ``None``

    Raises:
        ``OSError``

    """
    if not os.path.isdir(cache_dir):
        raise OSError('%s does not exist or is not a valid path' % cache_dir)
    cached_fn = sha256fn(hash_str) + '.npy'
    cached_fp = os.path.join(cache_dir, cached_fn)
    np.save(cached_fp, np_arr)


def buildIdxLUT(genome_fp, ref_genome_fp, genome_str=None, ref_genome_str=None,
                cache_dir=None):
    """Build a LUT that maps the indices of `genome` to `ref_genome`.

    Uses the mauve-py package (Python bindings for the mauve aligner) to
    generate the LUT and optionally handles caching of the LUT.

    Args:
        genome_fp (str)                 : filepath to the recoded/modified
                                          genome
        ref_genome_fp (str)             : filepath to the reference genome

        genome_str (str, optional)      : genome sequence (if already
                                          available) -- optimization
        ref_genome_str (str, optional)  : reference genome sequence (if already
                                          available) -- optimization
        cache_dir (str, optional)       : if not ``None``, the directory where
                                          a cached version of the LUT will be
                                          saved
    Returns:
        Numpy array mapping indices in the recoded/modified genome sequence
        to indicies in the reference genome sequence. All indices greater than
        -1 are properly mapped and will be unique. Indices with a value of -1
        did not map (and likely represent an insertion or region of the
        genome sequence that is not shared with the reference genome)

    Raises:
        ``OSError``
    """
    idx_lut = None
    if cache_dir is not None:
        cache_str = (genome_fp + str(os.path.os.path.getmtime(genome_fp)) +
                    ref_genome_fp + str(os.path.os.path.getmtime(ref_genome_fp))
                    + 'idxLUT')
        idx_lut = getCachedNumpyArray(cache_dir, cache_str)
    if idx_lut is None:
        idx_lut = mauve.buildIndex(genome_fp, ref_genome_fp, genome_str,
                                   ref_genome_str)
        if cache_dir is not None:
            saveNumpyArrayToCache(idx_lut, cache_dir, cache_str)
    return idx_lut


def buildEdgeLUT(idx_lut, cache_dir=None):
    r"""Build a binary lookup table of 'edges' in the ``idx_lut``.

    "Edges" are points at which there is a discontinuity in the index
    mapping. E.g., there is an "edge" at the star below:

        Genome idx:         200, 201, 202, 203, 204, 205, 206
        Ref. idx:           143, 144, 145, 189, 190, 191, 192
                                            *

    Args:
        idx_lut (``numpy.ndarray``) : numpy array of index mappings between two
                                      genomes (see :func:``buildIdxLUT``)

        cache_dir (str, optional)   : if not ``None``, the directory where a
                                      cached version of the LUT will be saved

    Returns:
        Array of equal length to ``idx_lut`` in which values of 1 represent
        "edges".

    """
    edge_lut = None
    if cache_dir is not None:
        cache_str = str(idx_lut) + 'edgeLUT'
        edge_lut = getCachedNumpyArray(cache_dir, cache_str)
    if edge_lut is None:
        edge_lut = mauve.indexutils.findEdges(idx_lut)
        if cache_dir is not None:
            saveNumpyArrayToCache(edge_lut, cache_dir, cache_str)
    return edge_lut


def buildMismatchLUT(idx_lut, genome_str, ref_genome_str, cache_dir=None):
    """Build a LUT of localized base mismatches in an inter-genome mapping.

    These are single base-pair mismatches or highly-localized mismatches
    (i.e., a run of 5 mismatches) that are useful for the purposes of
    identifying good points of sequence discrimination.

    Args:
        idx_lut (``numpy.ndarray``) : numpy array of index mappings between
                                      two genomes (see :func:``buildIdxLUT``)
        genome_str (str)            : recoded/modified genome sequence
        ref_genome_str (str)        : reference genome sequence

        cache_dir (str, optional)   : if not ``None``, the directory where a
                                      cached version of the LUT will be saved

    Returns:
        Array of equal length to ``idx_lut`` in which values of 1 represent
        localized mismatches.

    """
    mismatch_lut = None
    if cache_dir is not None:
        cache_str = str(genome_str + ref_genome_str) + 'mismatchLUT'
        mismatch_lut = getCachedNumpyArray(cache_dir, cache_str)
    if mismatch_lut is None:
        mismatch_lut = mauve.indexutils.findMismatches(idx_lut, genome_str,
                                                       ref_genome_str)
        if cache_dir is not None:
            saveNumpyArrayToCache(mismatch_lut, cache_dir, cache_str)
    return mismatch_lut
