
import atexit
import shutil
import os

from Bio import SeqIO


LOCAL_DIR = os.path.dirname(os.path.realpath(__file__))
TEST_OUTPUT_DIR = os.path.join(LOCAL_DIR, 'test_output')
TEST_CACHE_DIR = os.path.join(LOCAL_DIR, 'test_output', '.cache')

try:
    os.makedirs(TEST_CACHE_DIR)
except OSError as e:
    if e.errno != 17:
        raise

def rmCache():
    shutil.rmtree(TEST_CACHE_DIR)

atexit.register(rmCache)


RECODED_GENOME_FP = os.path.join(
    LOCAL_DIR, 'test_input', 'recoded_genome_1417977689.gb')


REFERENCE_GENOME_FP = os.path.join(
    LOCAL_DIR, 'test_input', 'reference_genome.gb')

RECODED_GB = SeqIO.read(RECODED_GENOME_FP, 'gb')
REFERENCE_GB = SeqIO.read(REFERENCE_GENOME_FP, 'gb')
RECODED_GB_STR = str(RECODED_GB.seq).encode('utf-8')
REFERENCE_GB_STR = str(REFERENCE_GB.seq).encode('utf-8')
