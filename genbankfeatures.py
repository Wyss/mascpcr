'''
genbankfeatures | genbankfeatures.py
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Tools to extract feature footprints / indices from genbank files.

The general idea is to generate bitmasks or index arrays for the purposes of
quick lookups during processes like MASC PCR developement.

'''
from __future__ import print_function

import re
import warnings

import bitarray

from Bio import SeqIO


def filterFeatures(gb_data, feature_types=None, feature_regexs=None):
    ''' Filter genbank features based on regular expressions in
    `feature_regexs` arranged in the following manner:

        <feature type regex>: <value regex>

    For example:

        'gene': 'thr.*'

    If `feature_types` is provided, only the specified feature types will
    be returned (e.g., 'CDS').

    Returns a list of SeqFeature objects from the SeqRecord object `gb_data`
    '''
    feature_regexs = feature_regexs or {}
    features = gb_data.features

    def _featureFilter(feature):
        if feature_types is not None and feature.type not in feature_types:
            return False
        for feat_key, feat_value_re in feature_regexs.values():
            print(feature.get(feat_key, ''))
            if re.search(feat_value_re, feature.get(feat_key, '')) is None:
                return False
        return True

    return filter(_featureFilter, features)


def generateMask(gb_data, feature_types=None, feature_regexs=None,
                 strand_specific=False):
    ''' Generate a bitarray mask using the given feature qualifiers (see
    `filterFeatures` docstring).

    If `strand_specific` is True, two separate bitarrays will be returned,
    one for the forward strand features, and one for the reverse strand
    features. Ambiguous features (i.e., strand == 0 or None) will be mapped
    to the forward strand. Sub features are not currently supported.

    Note: location.end is inclusive (so we need to add 1 for indexing purposes)
    '''
    filtered_features = filterFeatures(gb_data, feature_types, feature_regexs)
    if strand_specific:
        fwd_bt_arr = bitarray.bitarray(len(gb_data.seq))
        fwd_bt_arr.setall(0)
        rev_bt_arr = bitarray.bitarray(len(gb_data.seq))
        fwd_bt_arr.setall(0)
    else:
        bt_arr = bitarray.bitarray(len(gb_data.seq))
        bt_arr.setall(0)
    for feat in filtered_features:
        if feat.sub_features is None:
            warnings.warn('Feature {} not processed b/c sub features are '
                          'not currently supported'.format(repr(feat)))
            continue
        if strand_specific:
            bt_arr = rev_bt_arr if feat.location.strand == -1 else fwd_bt_arr
        bt_arr[feat.location.start:feat.location.end+1] = 1
    if strand_specific:
        return fwd_bt_arr, rev_bt_arr
    else:
        return bt_arr


def findBorders(gb_data, feature_types=None, feature_regexs=None,
                strand_specific=False):
    ''' Return a list of feature edges using the given feature qualifiers
    (see `filterFeatures` docstring). The list will be comprised of tuples
    with the following format:
        (index, 1)          where the second value is 1 if it's a 5' ('rising')
                            edge or 0 if it is a 3' ('falling') edge
    If `strand_specific` is True, two separate lists will be returned

    Note: location.end is inclusive (so we need to add 1 for indexing purposes)
    '''
    filtered_features = filterFeatures(gb_data, feature_types, feature_regexs)
    if strand_specific:
        fwd_feat_list = []
        rev_feat_list = []
    else:
        feat_list = []
    for feat in filtered_features:
        if strand_specific:
            if feat.location.strand == -1:
                feat_list = rev_feat_list
            else:
                feat_list = fwd_feat_list
        feat_list.append((feat.location.start, 1))
        feat_list.append((feat.location.end, 0))
    if strand_specific:
        return fwd_feat_list, rev_feat_list
    else:
        return feat_list

def genBorderArray(gb_data, feature_types=None, feature_regexs=None,
                   strand_specific=False):
    if strand_specific:
        fwd_feat_list, rev_feat_list = findBorders(gb_data, feature_types,
                                                   feature_regexs,
                                                   strand_specific)
        fwd_arr = bitarray.bitarray(len(gb_data.seq))
        fwd_arr.setall(0)
        rev_arr = bitarray.bitarray(len(gb_data.seq))
        rev_arr.setall(0)
        for rec in fwd_feat_list:
            fwd_arr[rec[0]] = 1
        for rec in rev_feat_list:
            rev_arr[rec[0]] = 1
        return fwd_arr, rev_arr
    else:
        feat_list = findBorders(gb_data, feature_types, feature_regexs,
                                strand_specific)
        feat_arr = bitarray.bitarray(len(gb_data.seq))
        feat_arr.setall(0)
        for rec in feat_list:
            try:
                feat_arr[rec[0]] = 1
            except IndexError:
                print('IndexError while generating border array {}'.format(rec[0]))
        return feat_arr


if __name__ == '__main__':
    gb_data = SeqIO.read(
        'test_data/2014_02_18_gen9_54_failed_seg_fixed_with_gc_fixes.gb', 'gb')
    # print(filterFeatures(gb_data, ['gen9_fragment']))
    # gen9_frags = filterFeatures(gb_data, ['gen9_fragment'])
    # generateMask(gb_data, ['gen9_fragment'])
    # generateMask(gb_data, ['gen9_fragment'], strand_specific=True)
    edge_list = findBorders(gb_data, ['gen9_fragment'])
    print('gen9_fragment edges')
    print(sum(genBorderArray(gb_data, ['gen9_fragment'])))
    # for edge in edge_list:
    #     print('\tidx: {:<8} | {}'.format(
    #             edge[0], 'rising' if edge[1] else 'falling'))

