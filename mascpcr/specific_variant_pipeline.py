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
mascpcr.specific_variant_pipeline
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Entry point for MASC-PCR primer design where user specifically provides a list
of positions and alternate alleles to assay.

NOTE: Work in progress. Not exposed via CLI yet.

"""

import copy
from itertools import cycle
import tempfile

from Bio import SeqIO
from Bio.Alphabet import generic_dna
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import pandas as pd

from mascpcr.ioutil import getPrimerPairProductSize
from mascpcr.pipeline import DEFAULT_PARAMS
from mascpcr.pipeline import generateLUTs
from mascpcr.primercandidate import findDiscriminatoryPrimer
from mascpcr.primerpartition import findBestCommonPrimerIn3pEndRange


EXPECTED_VARIANT_CSV_COLUMNS = ['POSITION', 'REF', 'ALT']

PRODUCT_SIZES = (100, 150, 200, 250, 300, 400, 500, 600, 700, 850)

OUTPUT_COLUMNS_ORDERED = [
    'POSITION', 'REF', 'ALT',
    'P_WT', 'P_WT_TM', 'P_WT_SCORE',
    'P_MUT', 'P_MUT_TM', 'P_MUT_SCORE',
    'P_COMMON', 'P_COMMON_TM', 'P_COMMON_SCORE',
    'AMPLICON_SIZE', 'P_FWD_STRAND']

# Defaults optimized for MASC-PCR according to Marc J. Lajoie.
THERMO_PARAMS = {
    'mv_conc': 50,
    'dv_conc': 1.5,
    'dntp_conc': 0.4,
    'dna_conc': 250,
}

SPURIOUS_TM_CLIP = 50

# More pressure on hitting the right value.
TM_RANGE = (61, 63)


def find_masc_pcr_primers_given_variants(
        variant_list_csv, ref_genome_src_path, output_file_csv):
    """Find a set of MASC-PCR primers given a list of variants as input.

    Args:
        variant_list_csv: csv file with at least columns
            EXPECTED_VARIANT_CSV_COLUMNS
        ref_genome_src_path: Full path to ref genome file. Genbank or fasta.
        output_file_csv: Path where output is written.
    """
    # Strategy: Leverage the genome-diffing machinery previously
    # implemented. We generate a modified version of the reference genome
    # and then request discriminatory primers at the specific positions we
    # are interested in.

    # Parse the input.
    variants_df = pd.read_csv(variant_list_csv)

    # Make sure required columns are present.
    assert not set(EXPECTED_VARIANT_CSV_COLUMNS) - set(variants_df.columns), (
            'Missing required columns ' +
            str(set(EXPECTED_VARIANT_CSV_COLUMNS) - set(variants_df.columns)))

    # Read in the ref genome as a SeqRecord, either genbank or fasta.
    ref_genome_seq_record = None
    for ref_genome_type in ['genbank', 'fasta']:
        try:
            with open(ref_genome_src_path) as fh:
                ref_genome_seq_record = SeqIO.read(fh, ref_genome_type)
        except ValueError:
            continue
    assert ref_genome_seq_record is not None, (
            "Ref genome not recognized as genbank or fasta.")
    ref_genome_str = str(ref_genome_seq_record.seq)

    # Generate a modified version of the reference genome which has the
    # desired changes incorporated. We do this so we can leverage the code
    # already implemented which takes as input the reference genome and modified
    # genome.
    modified_genome_str = _create_modified_genome_str(
            ref_genome_str, variants_df)

    # Write genomes to temp locations for rest of pipeline.
    with tempfile.NamedTemporaryFile(suffix='.gb', delete=False) as genome_fh:
        genome_fp = genome_fh.name
        modified_genome_seq_record = SeqRecord(
            Seq(modified_genome_str, generic_dna))
        SeqIO.write(modified_genome_seq_record, genome_fh, 'genbank')
    with tempfile.NamedTemporaryFile(suffix='.gb', delete=False) as ref_fh:
        ref_genome_fp = ref_fh.name
        ref_genome_seq_record = SeqRecord(
            Seq(ref_genome_str, generic_dna))
        SeqIO.write(ref_genome_seq_record, ref_fh, 'genbank')

    # DEBUG: TOGGLE WITH ABOVE
    # import os
    # from mascpcr.pipeline import CACHE_DIR
    # with open(os.path.join(CACHE_DIR, 'genome.gb'), 'w') as genome_fh:
    #     genome_fp = genome_fh.name
    #     modified_genome_seq_record = SeqRecord(
    #         Seq(modified_genome_str, generic_dna))
    #     SeqIO.write(modified_genome_seq_record, genome_fh, 'genbank')
    # with open(os.path.join(CACHE_DIR, 'ref_genome.gb'), 'w') as ref_fh:
    #     ref_genome_fp = ref_fh.name
    #     ref_genome_seq_record = SeqRecord(
    #         Seq(ref_genome_str, generic_dna))
    #     SeqIO.write(ref_genome_seq_record, ref_fh, 'genbank')
    # genome_fp = os.path.join(CACHE_DIR, 'genome.gb')
    # ref_genome_fp = os.path.join(CACHE_DIR, 'ref_genome.gb')

    # Generate the lookup tables (LUTs) used for the primer design pipeline.
    # Note that we use the reference genome as the 'FROM' genome and the
    # modified genome as the 'TO' genome.
    generate_LUTs_result = generateLUTs(
            ref_genome_fp,
            genome_fp,
            0, # start_idx
            len(modified_genome_str)) # end_idx
    (validate_from_genome_str, validate_to_genome_str, idx_lut, edge_lut,
            mismatch_lut, border_lut) = generate_LUTs_result
    assert validate_from_genome_str == ref_genome_str
    assert validate_to_genome_str == modified_genome_str

    # Next, for each variant, identify the best forward pair and reverse pair.
    # We'll choose whether to use the forward or reverse when we design the
    # common primer.

    # Modify params to be passed into primer finding functions.
    primer_finder_params = copy.deepcopy(DEFAULT_PARAMS)
    primer_finder_params.update({
        'min_num_mismatches': 1,
        'thermo_params': THERMO_PARAMS,
        'spurious_tm_clip': SPURIOUS_TM_CLIP,
        'tm_range': TM_RANGE
    })

    # We'll store the results as a map from variant position to discriminatory
    # primers.
    variant_position_to_candidate_primer_map = {}
    for _, v_obj in variants_df.iterrows():
        # Variants provided are 1-indexed, while pythonic bioinformatics stuff
        # is pretty much always 0-indexed like python.
        zero_indexed_variant_position = v_obj.POSITION - 1

        # First identify whether the allele position can be discriminated.
        # It's possible it can't, for example at a homopolymer run.
        # For now we set no result and continue to next variant.
        if mismatch_lut[zero_indexed_variant_position] != 1:
            variant_position_to_candidate_primer_map[v_obj.POSITION] = None
            continue

        # Next, we find the best discriminatory primer pairs for each variant.
        # The helper method _find_primers_strict_then_lenient() below allows us
        # to first search in strict mode and then lenient mode. This should
        # guarantee that a primer is found, even if not ideal.
        def _find_primers_strict_then_lenient(strand):
            def _find_disc_primer(maybe_mod_params):
                return findDiscriminatoryPrimer(
                        zero_indexed_variant_position, strand, idx_lut,
                        ref_genome_str, modified_genome_str, edge_lut,
                        mismatch_lut, maybe_mod_params)
            # First try in strict mode. Then lenient mode.
            result = _find_disc_primer(primer_finder_params)
            assert len(result) == 2
            if result[0] is None:
                # Repeat in lenient_mode.
                mod_params = copy.deepcopy(primer_finder_params)
                mod_params['lenient_mode'] = True
                result = _find_disc_primer(mod_params)
            return result

        # Compute and save all primer candidates.
        fwd_wt, fwd_mut = _find_primers_strict_then_lenient(1)
        rev_wt, rev_mut = _find_primers_strict_then_lenient(-1)
        variant_position_to_candidate_primer_map[v_obj.POSITION] = {
            'fwd_mut': fwd_mut,
            'fwd_wt': fwd_wt,
            'rev_mut': rev_mut,
            'rev_wt': rev_wt
        }

    # Now that we have primer candidates for each variant position, we want to
    # identify the best common primer to match.
    variant_position_to_best_primers_map = {}

    # For varying amplicon sizes, we just step through them in the order in
    # which the input is provided.
    amplicon_size_iter = cycle(PRODUCT_SIZES)

    for _, v_obj in variants_df.iterrows():
        candidate_dict = variant_position_to_candidate_primer_map[
                v_obj.POSITION]

        if candidate_dict is None:
            variant_position_to_best_primers_map[v_obj.POSITION] = None
            continue

        amplicon_size = next(amplicon_size_iter)

        # Find the best common primer.
        variant_position_to_best_primers_map[v_obj.POSITION] = (
                _find_best_common_primer(
                        candidate_dict, amplicon_size, idx_lut, ref_genome_str,
                        modified_genome_str, edge_lut, mismatch_lut,
                        primer_finder_params))

    # Finally write output report.
    # TODO: This should probably go in ioutil.

    result_list = []
    for _, v_obj in variants_df.iterrows():
        best_primer_set = variant_position_to_best_primers_map[v_obj.POSITION]
        if best_primer_set is not None:
            assert best_primer_set['wt'].strand == best_primer_set['mut'].strand
            assert (best_primer_set['wt'].strand !=
                    best_primer_set['common'].strand)
            v_obj['P_FWD_STRAND'] = best_primer_set['wt'].strand
            v_obj['AMPLICON_SIZE'] = getPrimerPairProductSize(
                    best_primer_set['wt'], best_primer_set['common'])
            report_prefix_and_candidate_key_pairs = [
                ('P_WT', 'wt'),
                ('P_MUT', 'mut'),
                ('P_COMMON', 'common')
            ]
            for prefix, key in report_prefix_and_candidate_key_pairs:
                v_obj[prefix] = best_primer_set[key].seq
                v_obj[prefix + '_TM'] = best_primer_set[key].tm
                v_obj[prefix + '_SCORE'] = best_primer_set[key].score
        result_list.append(v_obj)

    result_df = pd.DataFrame(result_list)

    # Explicitly set column order.
    result_df = result_df[OUTPUT_COLUMNS_ORDERED]

    result_df.to_csv(output_file_csv, index=False)


def _create_modified_genome_str(ref_genome_str, variants_df):
    """Generates a modified version of the reference genome which has the
    changes in variant_object_list incorporated.

    Args:
        ref_genome_str: String original genome seq.
        variants_df: pandas.DataFrame containing data about variants to
            incorporate.

    Returns:
        String genome seq with variants applied to ref_genome_str.
    """
    # Generate a modified version of the reference genome which has the
    # desired changes incorporated.
    # Strategy: Partition genome at variant positions. Make changes. Reassemble.
    # TODO: If we stick with using the two-genome strategy, we shoud make this
    # code more robust, perhaps using
    # https://github.com/churchlab/reference_genome_maker

    # Build this list of tuples (start_offset, partition). Any mutations
    # are inserted at the end of a partition.
    ordered_partition_obj_list = [(0, ref_genome_str[:])]
    variant_object_list = sorted(
            [pair[1] for pair in variants_df.iterrows()],
            key=lambda x: x.POSITION)
    # Track net change in genome length for assertion at the end.
    genome_len_delta_tracker = 0
    for v_obj in variant_object_list:
        # Pop the last partition object, split it at the breakpoint
        # corresponding to this variant, and add these two to the ordered list.
        last_partition_obj = ordered_partition_obj_list.pop()
        last_start_offset = last_partition_obj[0]
        last_seq_str = last_partition_obj[1]

        # Break just after the mutation.
        break_idx = v_obj.POSITION + len(v_obj.REF) - 1 - last_start_offset

        # Modify the pre-break sequence str and build the new partition object.
        is_ref_valid = bool(v_obj.REF ==
                last_seq_str[break_idx - len(v_obj.REF):break_idx])
        assert is_ref_valid, (
                'Variant at {pos} has incorrect REF'.format(pos=v_obj.POSITION))
        pre_break_seq_str = (
                last_seq_str[:break_idx - len(v_obj.REF)] + v_obj.ALT)
        pre_break_partition_obj = (last_start_offset, pre_break_seq_str)
        genome_len_delta_tracker += (len(v_obj.ALT) - len(v_obj.REF))

        # The post-break is what would normally follow.
        post_break_partition_obj = (
                break_idx + last_start_offset, last_seq_str[break_idx:])

        # Add the pre- and post-break partition objects back to the list.
        ordered_partition_obj_list.append(pre_break_partition_obj)
        ordered_partition_obj_list.append(post_break_partition_obj)

    # Assemble the modified genome.
    modified_genome_str = ''.join([p[1] for p in ordered_partition_obj_list])
    assert len(modified_genome_str) == (
            len(ref_genome_str) + genome_len_delta_tracker)

    return modified_genome_str


def _find_best_common_primer(
        candidate_dict, amplicon_size, idx_lut, from_genome_str,
        to_genome_str, edge_lut, mismatch_lut, primer_finder_params):
    """Identify the best set of mut, fwd, common primer given pairs of forward
    and reverse primer pairs.

    Slides across the possible common primer positions for this amplicon_size
    and returns the best scoring one.

    Args:
        candidate_dict: Dictionary with keys:
            * fwt_mut
            * fwd_wt
            * rev_mut
            * rev_wt

    Returns:
        Dictionary with keys:
            * wt
            * mut
            * common
        Each key points to value of type primercandidate.CandidatePrimer.
    """
    # Strategy: Identify common primer for each of the forward and reverse
    # candidates. Compare overall scores. Return the best mut/wt/common set.

    # Validate input according to assumptions for rest of code.
    for required_key in ['fwd_mut', 'fwd_wt', 'rev_mut', 'rev_wt']:
        assert required_key in candidate_dict
        assert candidate_dict[required_key] is not None

    # Forward.
    fwd_common_primer_end3p_range = _find_fwd_common_primer_end3p_range(
            candidate_dict, amplicon_size, edge_lut, primer_finder_params)
    best_common_primer_for_fwd_primers = findBestCommonPrimerIn3pEndRange(
        fwd_common_primer_end3p_range, -1, candidate_dict['fwd_mut'],
        from_genome_str, to_genome_str, idx_lut, edge_lut, mismatch_lut,
        primer_finder_params, amplicon_size=amplicon_size)

    # Reverse.
    rev_common_primer_end3p_range = _find_rev_common_primer_end3p_range(
            candidate_dict, amplicon_size, edge_lut, primer_finder_params)
    best_common_primer_for_rev_primers = findBestCommonPrimerIn3pEndRange(
        rev_common_primer_end3p_range, 1, candidate_dict['rev_mut'],
        from_genome_str, to_genome_str, idx_lut, edge_lut, mismatch_lut,
        primer_finder_params, amplicon_size=amplicon_size)

    # TODO: It should be okay if one of these is None, but need to explore
    # failure cases more carefully.
    assert (best_common_primer_for_fwd_primers is not None or
            best_common_primer_for_rev_primers is not None)

    def _build_fwd_result():
        return {
            'wt': candidate_dict['fwd_wt'],
            'mut': candidate_dict['fwd_mut'],
            'common': best_common_primer_for_fwd_primers
        }

    def _build_rev_result():
        return {
            'wt': candidate_dict['rev_wt'],
            'mut': candidate_dict['rev_mut'],
            'common': best_common_primer_for_rev_primers
        }

    # If one is None, return the other.
    if best_common_primer_for_fwd_primers is None:
        assert best_common_primer_for_rev_primers is not None
        return _build_rev_result()
    elif best_common_primer_for_rev_primers is None:
        assert best_common_primer_for_fwd_primers is not None
        return _build_fwd_result()

    # Otherwise compares scores, calculated as sum of the three primer scores.
    # TODO: Is there a better way to come up with the total score?
    fwd_score = (
            best_common_primer_for_fwd_primers.score +
            candidate_dict['fwd_mut'].score +
            candidate_dict['fwd_wt'].score)
    rev_score = (
            best_common_primer_for_rev_primers.score +
            candidate_dict['rev_mut'].score +
            candidate_dict['rev_wt'].score)

    if fwd_score >= rev_score:
        return _build_fwd_result()
    else:
        return _build_rev_result()


def _find_fwd_common_primer_end3p_range(
        candidate_dict, amplicon_size, edge_lut, primer_finder_params):
    """Helper function for finding the end3p range.
    """
    # Identify the range of possible 3-prime ends for the common primer using
    # the 5-prime end of discriminatory_primer and amplicon_size.
    discriminatory_primer_end5p = candidate_dict['fwd_mut'].idx

    # e.g. 0 + 100 - 30 - 10 = 60
    common_primer_end3p_min = (
            discriminatory_primer_end5p +
            amplicon_size -
            primer_finder_params['size_range'][1] -
            primer_finder_params['product_size_tolerance'])

    # e.g. 0 + 100 - 18 + 10 = 82
    common_primer_end3p_upper_bound = (
            discriminatory_primer_end5p +
            amplicon_size -
            primer_finder_params['size_range'][0] +
            primer_finder_params['product_size_tolerance'])

    # Define the search space of common primer 3-prime ends to avoid edges.
    common_primer_end3p_max = common_primer_end3p_min
    for safe_upper_bound in range(common_primer_end3p_min + 1,
            common_primer_end3p_upper_bound):
        if edge_lut[safe_upper_bound] != 0:
            break
        else:
            common_primer_end3p_max = safe_upper_bound

    return range(common_primer_end3p_min, common_primer_end3p_max)


def _find_rev_common_primer_end3p_range(
        candidate_dict, amplicon_size, edge_lut, primer_finder_params):
    """Helper function for finding the end3p range.
    """
    # Identify the range of possible 3-prime ends for the common primer using
    # the 5-prime end of discriminatory_primer and amplicon_size.
    discriminatory_primer_end5p = (candidate_dict['rev_mut'].idx +
            candidate_dict['rev_mut'].length)

    # e.g. 1000 - 100 + 30 + 10
    common_primer_end3p_max = (
            discriminatory_primer_end5p -
            amplicon_size +
            primer_finder_params['size_range'][1] +
            primer_finder_params['product_size_tolerance'])

    # e.g. 0 + 100 - 18 + 10 = 82
    common_primer_end3p_lower_bound = (
            discriminatory_primer_end5p -
            amplicon_size +
            primer_finder_params['size_range'][0] -
            primer_finder_params['product_size_tolerance'])

    # Define the search space of common primer 3-prime ends to avoid edges.
    min3p_no_edge = common_primer_end3p_max
    for min3p_no_edge in range(
            common_primer_end3p_max - 1, common_primer_end3p_lower_bound, -1):
        if edge_lut[min3p_no_edge] != 0:
            min3p_no_edge += 1
            break

    return range(common_primer_end3p_max, min3p_no_edge - 1, -1)
