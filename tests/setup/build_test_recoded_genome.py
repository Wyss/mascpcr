
import os
import random
import sys
import time

from Bio import SeqIO

from Bio.Alphabet import IUPAC
from Bio.Seq import MutableSeq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation

LOCAL_DIR = os.path.dirname(os.path.realpath(__file__))
TEST_DIR = os.path.dirname(LOCAL_DIR)
REF_GENOME_FP = os.path.join(TEST_DIR, 'test_input', 'reference_genome.gb')

ref_genome_obj = SeqIO.read(REF_GENOME_FP, 'gb')

rec_genome_seq = MutableSeq(str(ref_genome_obj.seq), 
                            IUPAC.IUPACUnambiguousDNA())

MUTATION_RATE = 0.03
PROXIMAL_MUTATION_RATE = 0.25

prev_base_mutated = False
num_mutated_bases = 0

genome_length = len(rec_genome_seq)

# minor recoding simulation
for idx in range(genome_length):
    if prev_base_mutated:
        mutate_base = random.random() < PROXIMAL_MUTATION_RATE
    else:
        mutate_base = random.random() < MUTATION_RATE
    if mutate_base:
        rec_genome_seq[idx] = random.choice(
          'ATGC'.replace(rec_genome_seq[idx], ''))
        num_mutated_bases += 1
        prev_base_mutated = True
    else:
        prev_base_mutated = False
    sys.stdout.write('Base: %d/%d, Mutation rate: %f\r' % (idx, genome_length, 
                   float(num_mutated_bases)/(idx+1) * 100))
    sys.stdout.flush()

# major swap simulation
for x in range(100):
    swap_size = random.randint(20, 4000)
    swap_idx_1 = random.randint(0, genome_length-swap_size-1)
    swap_idx_2 = random.randint(0, genome_length-swap_size-1)
    while swap_idx_1 < swap_idx_2 < swap_idx_1 + swap_size:
        swap_idx_2 = random.randint(0, genome_length-swap_size-1)
    piece_1 = rec_genome_seq[swap_idx_1:swap_idx_1+swap_size]
    piece_2 = rec_genome_seq[swap_idx_2:swap_idx_2+swap_size]
    rec_genome_seq[swap_idx_1:swap_idx_1+swap_size] = piece_2
    rec_genome_seq[swap_idx_2:swap_idx_2+swap_size] = piece_1


rec_genome = SeqRecord(rec_genome_seq, id='001', name='recoded_genome')

SYNTH_SEG_SIZE = 49000
SYNTH_FRAG_SIZE = 2500

for seg_num in range(genome_length / SYNTH_SEG_SIZE):
    start_idx = seg_num * SYNTH_SEG_SIZE
    end_idx = start_idx + SYNTH_SEG_SIZE
    if seg_num == (genome_length / SYNTH_SEG_SIZE - 1):
        end_idx += genome_length % SYNTH_SEG_SIZE
    s_feat = SeqFeature(FeatureLocation(start_idx, end_idx), 
                        type='synth_segment', qualifiers={'label':
                        'seg%02d' % (seg_num)})
    rec_genome.features.append(s_feat)
    for frag_num in range(SYNTH_SEG_SIZE / SYNTH_FRAG_SIZE):
        f_start_idx = start_idx + SYNTH_FRAG_SIZE * frag_num
        f_end_idx = f_start_idx + SYNTH_FRAG_SIZE
        if frag_num == (SYNTH_SEG_SIZE / SYNTH_FRAG_SIZE - 1):
            f_end_idx += SYNTH_SEG_SIZE % SYNTH_FRAG_SIZE
        else: 
            f_end_idx += 100  # simulated overlap for assembly
        f_feat = SeqFeature(FeatureLocation(f_start_idx,f_end_idx), 
                            type='synth_fragment', qualifiers={'label':
                            'seg%02d_%03d' % (seg_num, frag_num)})
        rec_genome.features.append(f_feat)

with open('recoded_genome_%d.gb' % time.time(), 'w') as fd:
    SeqIO.write(rec_genome, fd, 'gb')

