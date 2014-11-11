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
mascpcr.genbankfeatures
~~~~~~~~~~~~~~~~~~~~~~~

:copyright: (c) 2014 Ben Pruitt & Nick Conway.
:license: GPLv2, see LICENSE for more details.

Methods to extract feature footprints / indices from instantiated 
``SeqRecord`` objects (``Biopython`` objects that abstract the contents of a
genbank file).

"""
from __future__ import print_function

import re

import bitarray


def filterFeatures(sr_obj, feature_types=None, qualifier_regexs=None):
    """Filter a `SeqRecord` object's `SeqFeature` list by type and qualifiers.


    Args:
        ``sr_obj``: instantiated Biopython ``SeqRecord`` object

    Kwargs:
        ``feature_types``: list of feature types (e.g., ['gene', 'CDS'])
        ``qualifier_regexs``: dict of <field name>: <value regex> entries

    Returns:
        Filtered list of `SeqRecord` objects

    Raises:
        None


    Examples:

        Return a list of all ``SeqFeature`` objects from ``gb_rec`` 
        that are of type 'mRNA' or 'CDS'::
    
            >>>filterFeatures(gb_rec, ['mRNA', 'CDS'])

        Return a list of all ``SeqFeature`` objects from ``gb_rec`` 
        that are of type 'mRNA' or 'CDS' and that additionally have the
        qualifier field 'gene' with a value that matches the regular expression
        'ubc.*'::

            >>>filterFeatures(gb_rec, ['gene', 'CDS'], {'gene': 'ubc.*'})

        The latter example would match the following genbank records::

            CDS             join(54..567,789..1254)
                            /gene="ubc42"
                            /product="ubiquitin conjugating enzyme"
                            /function="cell division control"

            CDS             join(54..567,789..1254)
                            /gene="ubc51"
                            /product="ubiquitin conjugating enzyme"
                            /function="cell division control"

    """
    qualifier_regexs = qualifier_regexs or {}
    features = sr_obj.features

    def _featureFilter(feature):
        if feature_types is not None and feature.type not in feature_types:
            return False
        for feat_key, feat_value_re in qualifier_regexs.items():
            q_values = feature.qualifiers.get(feat_key, [])
            for v in q_values:
                if re.search(feat_value_re, v) is None:
                    return False
        return True

    return filter(_featureFilter, features)


def findBorders(sr_obj, feature_types=None, qualifier_regexs=None,
                strand_specific=False):
    """Filter a ``SeqFeature`` list and find the border indices of its members.

    See :func:`filterFeatures` for explanation of filtering functionality.


    Args:
        ``sr_obj``: instantiated Biopython ``SeqRecord`` object

    Kwargs:
        ``feature_types``: list of feature types (e.g., ['gene', 'CDS'])
        ``qualifier_regexs``: dict of <field name>: <value regex> entries
        ``strand_specific``: boolean determining whether separate lists 
                             should be returned for each strand (fwd / rev)
    Returns:
        List(s) of (<idx>, <1 if rising edge, 0 if falling edge>) tuples.

    Raises:
        None

    """
    filtered_features = filterFeatures(sr_obj, feature_types, qualifier_regexs)
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


def buildBorderLUT(sr_obj, feature_types=None, qualifier_regexs=None,
                   strand_specific=False):
    """Filter a ``SeqRecord``'s features and build a binary LUT of border edges.

    See :func:`filterFeatures` for explanation of filtering functionality.

    Args:
        ``sr_obj``: instantiated Biopython ``SeqRecord`` object

    Kwargs:
        ``feature_types``: list of feature types (e.g., ['gene', 'CDS'])
        ``qualifier_regexs``: dict of <field name>: <value regex> entries
        ``strand_specific``: boolean determining whether separate lists 
                             should be returned for each strand (fwd / rev)
    Returns:
        Binary bitarray(s) (``bitarray.bitarray``) indicating the indices of
        feature borders (border indices have a value of 1). Strand-specific
        bitarrays are returned if ``strand_specific`` is ``True``.

    Raises:
        None

    """
    if strand_specific:
        fwd_feat_list, rev_feat_list = findBorders(sr_obj, feature_types,
                                                   qualifier_regexs,
                                                   strand_specific)
        fwd_lut = bitarray.bitarray(len(sr_obj.seq))
        fwd_lut.setall(0)
        rev_lut = bitarray.bitarray(len(sr_obj.seq))
        rev_lut.setall(0)
        for rec in fwd_feat_list:
            fwd_lut[rec[0]] = 1
        for rec in rev_feat_list:
            rev_lut[rec[0]] = 1
        return fwd_lut, rev_lut
    else:
        feat_list = findBorders(sr_obj, feature_types, qualifier_regexs,
                                strand_specific)
        feat_lut = bitarray.bitarray(len(sr_obj.seq))
        feat_lut.setall(0)
        for rec in feat_list:
            try:
                feat_lut[rec[0]] = 1
            except IndexError:
                print('IndexError while generating border array {}'.format(
                      rec[0]))
        return feat_lut


def findAggregateBoundaries(sr_obj, feature_types=None, qualifier_regexs=None):
    """Determine the outermost border indices of a group of filtered features.

    See :func:`filterFeatures` for explanation of filtering functionality.


    Args:
        ``sr_obj``: instantiated Biopython ``SeqRecord`` object

    Kwargs:
        ``feature_types``: list of feature types (e.g., ['gene', 'CDS'])
        ``qualifier_regexs``: dict of <field name>: <value regex> entries

    Returns:
        Tuple of (<min index>, <max index>) of the filtered features

    Raises:
        None


    For example, let's say your genbank file has the following features:

        synth_seg   1001..2000
                        /label="seg17_000"
        synth_seg   2001..3000
                        /label="seg17_001"
        synth_seg   3001..4000
                        /label="seg17_002"
        synth_seg   4001..5000
                        /label="seg18_000"

        Then the following call will produce this output::

        >>>findAggregateBoundaries(sr, ['synth_seg'], {'label': r'seg17.*'})
        (1001, 4000)

    """
    filtered_features = list(filterFeatures(sr_obj, feature_types, 
                                            qualifier_regexs))
    if len(filtered_features) == 0:
        return None, None
    min_idx = len(sr_obj.seq) + 1
    max_idx = -1
    for ff in filtered_features:
        min_idx = min(int(ff.location.start), min_idx)
        max_idx = max(int(ff.location.end), max_idx)
    return min_idx, max_idx
