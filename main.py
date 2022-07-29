from bloom_filter2 import BloomFilter
from collections import defaultdict
from pathlib import Path
from tqdm import tqdm
from math import floor, e
import json
import os


INPUT_DIR = os.path.join(".", "input")                  # Directory for RNG output
OUTPUT_DIR = os.path.join(".", "output")                # Directory for BitReps JSON
RESULTS_DIR = os.path.join(".", "results")              # Directory for BitReps analysis results
MODEL_DIR = os.path.join(".", "model")                  # Directory for baseline chi-square distribution
POSSIBLE_BLKS = [8, 16, 32, 64, 128, 256, 512]          # Supported blocksizes for BitReps


def dir_setup():
    """
    Create the directories necessary for BitReps operation, if necessary
    :return: None
    """
    if not os.path.isdir(OUTPUT_DIR):
        os.mkdir(OUTPUT_DIR)
    if not os.path.isdir(RESULTS_DIR):
        os.mkdir(RESULTS_DIR)


def calc_current_fpr(k, m, n):
    """
    Calculate the probability of a false positive when inserting element number n into a bloom filter composed of k
    hashes and m bits.
    :param k: Number of hashes
    :param n: Nth element
    :param m: Number of bits
    :return: False positive rate
    """
    return (1-e**(-k*n/m))**k


def generate_slides(bin_x, bin_y):
    """
    Generate all windows across the concatenation of two bitstrings.
    :param bin_x: First bitstring.
    :param bin_y: Second bitstring.
    :return: A list of all windows across the given bitstrings.
    """
    slides = []
    combined = bin_x + bin_y
    for i in range(len(bin_x)):
        slides.append(combined[i:i+len(bin_x)])
    print(slides)
    return slides


def slide_blocks(blocks, blocksize):
    """
    Convert list of non-overlapping blocks into overlapping blocks.
    :param blocks: List of non-sliding blocks built from input data.
    :param blocksize: User-specified blocksize.
    :return: A list of overlapping blocks based on user-specified input data.
    """
    completed_blocks = []
    aux = iter(blocks)
    for x, y in zip(aux, aux):
        print(x, y)
        bin_x = "{0:0{blocksize}b}".format(x, blocksize=blocksize)
        bin_y = "{0:0{blocksize}b}".format(y, blocksize=blocksize)
        completed_blocks.extend(generate_slides(bin_x, bin_y))
        print()
    return completed_blocks


def get_blocks(input_data, blocksize):
    """
    Split input data into blocks.
    :param input_data: User-specified input data.
    :param blocksize: User-specified blocksize.
    :return:
    """
    blocks = []
    with open(input_data, "rb") as f:
        while True:
            block = f.read(int(blocksize/8))
            if not block:
                break
            blocks.append(int.from_bytes(block, "big"))
    return blocks


def get_num_blocks(input_file, blocksize):
    """
    Determine the number of blocks in input data based on specified blocksize.
    :param input_file: User-specified input data.
    :param blocksize: User-specified blocksize.
    :return:
    """
    size = os.path.getsize(input_file)
    blocksize_bytes = blocksize / 8
    return size / blocksize_bytes


def tracker_dict():
    """
    Returns dictionary to be used as default dictionary.
    :return: Template for output dictionary.
    """
    return {
        "num_reps": 0,
        "bin_rep": "",
        "indices": []
    }


def bitreps_measure(input_file, blocksize, progress_bar, sliding, err_rate):
    """
    Orchestrate the BitReps test for a given input
    :param input_file: Path of input data
    :param blocksize: Desired blocksize for BitReps test
    :param progress_bar: Representation of progress bar (passed from GUI)
    :param sliding: No sliding window (0) or sliding window (1)
    :param err_rate: Desired error rate for the underlying bloom filter
    :return: None
    """
    # Ensure chosen blocksize is valid
    try:
        assert blocksize in POSSIBLE_BLKS
    except AssertionError:
        print("Invalid blocksize! Must be 16, 32, 64, 128, 256 or 512.")
        exit(1)

    # Ensure chosen error rate is valid
    try:
        assert 0 <= err_rate <= 1
    except AssertionError:
        print("Invalid error rate! Must be a float between 0 and 1 inclusive.")
        exit(1)

    blocks = get_blocks(input_file, blocksize)              # Split input data into blocks

    if sliding:                                             # If the user specifies a sliding window
        blocks = slide_blocks(blocks, blocksize)            # Convert blocks into sliding blocks

    num_blocks = len(blocks)                                        # Number of blocks in the input data
    bf = BloomFilter(max_elements=num_blocks, error_rate=err_rate)  # Instantiate bloom filter data structure
    percent = len(blocks) / 100                                     # Determine percentage increment requirement for GUI

    outer_hits = {
        "hits": defaultdict(tracker_dict),
        "blocksize": blocksize,
        "sliding": sliding,
        "err_rate": err_rate,
        "num_blocks": num_blocks,
        "avg_err_rate": 0
    }

    fprs = []

    for i in tqdm(range(len(blocks))):                      # For each block
        if blocks[i] in bf:                                 # If the block is in the bloom filter
            outer_hits["hits"][blocks[i]]["num_reps"] += 1  # Increase the number of observed repetitions for this block
            outer_hits["hits"][blocks[i]]["indices"].append(i)  # Associate current block index with nth repetition
            bin_rep = "{0:b}".format(blocks[i])
            while len(bin_rep) < blocksize:                 # Pad binary representation to blocksize if required
                bin_rep = "0" + bin_rep

            outer_hits["hits"][blocks[i]]["bin_rep"] = bin_rep
        else:                                               # If the block was not in the bloom filter
            bf.add(blocks[i])                               # Add it to the bloom filter, but not as a repetition

        # Record FPR for current insertion
        fprs.append(calc_current_fpr(bf.num_probes_k, bf.num_bits_m, i+1))

        # Increment the progress bar (for the GUI)
        if i % floor(percent) == 0:
            progress_bar.setValue(i / percent)

    # Determine average FPR across insertion time and add this to metadata
    afpr = sum(fprs) / len(fprs)
    outer_hits["avg_err_rate"] = afpr

    # Obtain input file name in preparation for output file
    output_name = Path(input_file).stem

    # Write the output JSON
    with open(os.path.join(OUTPUT_DIR, "%s-%s-%s-%s.json" % (output_name, blocksize, str(err_rate).replace(".", "_"), sliding)), "w+") as of:
        json.dump(outer_hits, of, indent=2)

    # Complete the progress bar
    progress_bar.setValue(100)
