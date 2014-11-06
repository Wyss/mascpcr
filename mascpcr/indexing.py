'''
indexing | indexing.py
~~~~~~~~~~~~~~~~~~~~~~

Contains methods related to generating index LUTs between two sequences as well
as the respective bitarrays that represent "edges" and single-base mismatches
between the two sequences. 

'''
import hashlib
import os
import six

import bitarray
import numpy as np
import mauve


def sha256fn(string):
    ''' Generate a filename-friendly sha256 hash of the provided string '''
    raw_hash = hashlib.sha256()
    if isinstance(string, six.string_types):
        raw_hash.update(string.encode('utf-8'))
    else:   # assume bytes (PY3) or str (PY2)
        raw_hash.update(string)
    raw_fn = raw_hash.hexdigest()
    # Remove forbidden / problematic characters for filenames
    raw_fn.replace('/', '_').replace('\\', '_').replace(':', '-')
    fn = 'cache' + raw_fn
    return fn


def getCachedNumpyArray(cache_dir, hash_str):
    ''' Look to see if `cache_dir` contains a file with a filename that is
    the sha256 hash of `hash_str` with a ".npy" extension. If it exists,
    attempt to load the data into a numpy array. If the loading operation
    fails or the file does not exist, return `None`.
    '''
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
    ''' Save a numpy array to `cache_dir` by generating a filename that
    consists of the sha256 hash of `hash_str` with a ".npy" extension.
    '''
    cached_fn = sha256fn(hash_str) + '.npy'
    cached_fp = os.path.join(cache_dir, cached_fn)
    np.save(cached_fp, np_arr)  


def getCachedBitArray(cache_dir, hash_str):
    ''' Look to see if `cache_dir` contains a file with a filename that is
    the sha256 hash of `hash_str` with a ".bt" extension. If it exists,
    attempt to load the data into a bitarray. If the loading operation
    fails or the file does not exist, return `None`.
    '''
    cached_fn = sha256fn(hash_str) + '.bt'
    cached_fp = os.path.join(cache_dir, cached_fn)
    bt_arr = None
    if os.path.isfile(cached_fp):
        try:
            bt_arr = bitarray.bitarray()
            with open(cached_fp, 'rb+') as cached_fd:
                bt_arr.fromfile(cached_fd)
        except OSError:
            bt_arr = None
            pass
    return bt_arr


def saveBitArrayToCache(bt_arr, cache_dir, hash_str):
    ''' Save a numpy array to `cache_dir` by generating a filename that
    consists of the sha256 hash of `hash_str` with a ".bt" extension.
    '''
    cached_fn = sha256fn(hash_str) + '.bt'
    cached_fp = os.path.join(cache_dir, cached_fn)
    with open(cached_fp, 'wb+') as cached_fd:
        bt_arr.tofile(cached_fd)  


def buildIdxLUT(genome_gb, ref_genome_gb, genome_seq=None, ref_genome_seq=None,
                cache_dir=None):
    ''' Build a LUT that maps the indices of `genome` to `ref_genome` using
    the mauve-py package (Python bindings for the mauve aligner).
    '''
    idx_lut = None
    if cache_dir is not None:
        cache_str = (genome_gb + str(os.path.os.path.getmtime(genome_gb)) + 
                    ref_genome_gb + str(os.path.os.path.getmtime(ref_genome_gb))
                    + 'idxLUT')
        idx_lut = getCachedNumpyArray(cache_dir, cache_str)
    if idx_lut is None:
        idx_lut = mauve.buildIndex(genome_gb, ref_genome_gb, genome_seq, 
                                   ref_genome_seq)
        if cache_dir is not None:
            saveNumpyArrayToCache(idx_lut, cache_dir, cache_str)
    return idx_lut


def buildEdgeLUT(idx_lut, cache_dir=None):
    ''' Build a binary lookup table of "edges" in the `idx_lut`.

    "Edges" are points at which there is a discontinuity in the index 
    mapping. E.g., there is an "edge" at the star below:

        Genome idx:         200, 201, 202, 203, 204, 205, 206
        Ref. idx:           143, 144, 145, 189, 190, 191, 192
                                            *
    '''
    edge_lut = None
    if cache_dir is not None:
        cache_str = str(idx_lut) + 'edgeLUT'
        edge_lut = getCachedBitArray(cache_dir, cache_str)
    if edge_lut is None:
        edge_lut = mauve.indexutils.findEdges(idx_lut)
        if cache_dir is not None:
            saveBitArrayToCache(edge_lut, cache_dir, cache_str)
    return edge_lut


def buildMismatchLUT(idx_lut, genome_str, ref_genome_str, cache_dir=None):
    ''' Build a binary lookup table of mismatches between the genome and
    reference genome, indexed by the base index of the genome. These are
    single base-pair mismatches or highly-localized mismatches (i.e., a 
    run of 5 mismatches) that are useful for the purposes of identifying
    good disciminatory primers.
    '''
    mismatch_lut = None
    if cache_dir is not None:
        cache_str = str(genome_str + ref_genome_str) + 'mismatchLUT'
        mismatch_lut = getCachedBitArray(cache_dir, cache_str)
    if mismatch_lut is None:
        mismatch_lut = mauve.indexutils.findMismatches(idx_lut, genome_str, 
                                                       ref_genome_str)
        if cache_dir is not None:
            saveBitArrayToCache(mismatch_lut, cache_dir, cache_str)
    return mismatch_lut    
