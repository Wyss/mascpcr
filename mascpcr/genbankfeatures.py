'''
genbankfeatures | genbankfeatures.py
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Tools to extract feature footprints / indices from genbank files.

'''
from __future__ import print_function

import re

import bitarray


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
        for feat_key, feat_value_re in feature_regexs.items():
            q_values = feature.qualifiers.get(feat_key, [])
            for v in q_values:
                if re.search(feat_value_re, v) is None:
                    return False
        return True

    return filter(_featureFilter, features)


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


def buildBorderLUT(gb_data, feature_types=None, feature_regexs=None,
                   strand_specific=False):
    if strand_specific:
        fwd_feat_list, rev_feat_list = findBorders(gb_data, feature_types,
                                                   feature_regexs,
                                                   strand_specific)
        fwd_lut = bitarray.bitarray(len(gb_data.seq))
        fwd_lut.setall(0)
        rev_lut = bitarray.bitarray(len(gb_data.seq))
        rev_lut.setall(0)
        for rec in fwd_feat_list:
            fwd_lut[rec[0]] = 1
        for rec in rev_feat_list:
            rev_lut[rec[0]] = 1
        return fwd_lut, rev_lut
    else:
        feat_list = findBorders(gb_data, feature_types, feature_regexs,
                                strand_specific)
        feat_lut = bitarray.bitarray(len(gb_data.seq))
        feat_lut.setall(0)
        for rec in feat_list:
            try:
                feat_lut[rec[0]] = 1
            except IndexError:
                print('IndexError while generating border array {}'.format(
                      rec[0]))
        return feat_lut


def findAggregateBoundaries(gb_data, feature_types=None, feature_regexs=None):
    ''' Find and return the min and max indices of a group of features matching
    `feature_types` and/or `feature_regexs`. For example, let's say your 
    genbank file has the following features:

    synth_seg   1001..2000
                    /label="seg17_000"
    synth_seg   2001..3000
                    /label="seg17_001"
    synth_seg   3001..4000
                    /label="seg17_002"
    synth_seg   4001..5000
                    /label="seg18_000"

    Then the following call will produce this output:

    >>>findAggregateBoundaries(gb_data, ['synth_seg'], {'label': r'seg17.*'})
    (1001, 4000)

    '''
    filtered_features = list(filterFeatures(gb_data, feature_types, 
                                            feature_regexs))
    if len(filtered_features) == 0:
        return None, None
    min_idx = len(gb_data.seq) + 1
    max_idx = -1
    for ff in filtered_features:
        min_idx = min(int(ff.location.start), min_idx)
        max_idx = max(int(ff.location.end), max_idx)
    return min_idx, max_idx
