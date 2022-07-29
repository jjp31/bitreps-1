from collections import Counter
from main import MODEL_DIR
from decimal import *
import json
import os


EXPECTED_PATH = os.path.join(MODEL_DIR, "urandom100M-32-1e-05-False.json")


def custom_chi(obs, exp):
    """
    Calculate the chi-square value of the provided observed and expected distributions
    :param obs: Observed distribution
    :param exp: Expected distribution
    :return: Chi-square value
    """
    chisq = 0
    for i in range(len(obs)):
        chisq += ((obs[i] - exp[i])**2) / exp[i]
    return chisq


def get_ratio(num_obs, num_exp):
    """
    Calculate the ratio of observed repetitions compared to expected repetitions, as a basic indicator of the amount
    of repetition in an input file.
    :param num_obs: The number of observed repetitions
    :param num_exp: The number of expected repetitions (due to false positives and genuine repetitions)
    :return: A value representing the level of repetition within an input, tending to 1.0 for random data
    """
    return round((num_exp / num_obs), 2)


def get_exp_fps(num_blocks, avg_err_rate):
    """
    Determine the number of expected false positives given the bloom filter's target error rate
    :param num_blocks: Number of blocks in the input data
    :param avg_err_rate: Average false positive rate of bloom filter across insertion time
    :return: Expected number of false positives
    """
    return round(num_blocks * avg_err_rate)


def get_exp_dupes(num_blocks, blocksize):
    """
    Determine the expected number of non-unique blocks assuming a uniform distribution of blocks.
    Note the use of Decimal to ensure small numbers are not rounded to zero. Without this, in the case of large
    block sizes, the result of the division (1 / x) is rounded to zero, thus invalidating output. Unfortunately, in the
    event of much larger blocksizes (e.g. 256 bits), the division results in such a small number that it rounds to zero
    even using Decimal. The implication here is that all blocks are expected to be unique - as such, zero is
    automatically returned in these cases.

    :param num_blocks: The number of blocks in the output
    :param blocksize: The size of the blocks, in bits
    :return: The expected number of non-unique blocks in the output
    """
    n = num_blocks
    x = 2 ** blocksize

    if 1 - Decimal(1 / x) == 1:     # All blocks are expected to be unique, therefore return 0
        return 0
    else:
        return round(num_blocks - (x * (1 - (1 - Decimal(1 / x)) ** n)))


def get_highest_rep(inputfile):
    """
    Determine the highest individually-repeating block within the output
    :param inputfile: File path of JSON output
    :return: The number of times that the maximally-repeating block within the output occurs
    """
    with open(inputfile) as f:
        data = json.load(f)

    current_max = 0
    block = 0
    for k, v in data["hits"].items():
        if v["num_reps"] > current_max:
            current_max = v["num_reps"]
            block = v["bin_rep"]
    return "%s (%s)" % (current_max, block)


def get_meta_data(inputfile):
    """
    Read BitReps JSON output and obtain metadata relating to blocksize, window type and error rate
    :param inputfile: BitReps JSON output supplied by user
    :return: A tuple of (blocksize, sliding and err_rate)
    """
    with open(inputfile) as f:
        data = json.load(f)
    blocksize = data["blocksize"]
    sliding = data["sliding"]
    err_rate = data["err_rate"]
    num_blocks = data["num_blocks"]
    obs_hits = len(data["hits"])
    avg_err_rate = data["avg_err_rate"]
    return blocksize, sliding, err_rate, num_blocks, obs_hits, avg_err_rate


def get_distri(inputfile):
    """
    Obtain the distribution of repetitions from a BitReps JSON file to be used for later processing
    :param inputfile: Path of BitReps JSON file
    :return: A list of values representing the number of repetitions in a file
    """
    distri = []
    with open(inputfile) as f:
        data = json.load(f)
    for k, v in data["hits"].items():
        distri.append(v["num_reps"])
    return distri


def trim_expected(distri):
    """
    Ensure that each histogram bucket has, at minimum, 5 elements, otherwise remove that bucket from the expected
    distribution
    :param distri: Expected distribution
    :return: Expected distribution where each element occurs at a minimum of 5 times
    """
    counts = Counter(distri)
    trimmed_distri = []
    unwanted = []
    for k, v in counts.items():
        if v < 5:
            unwanted.append(k)
    for num in distri:
        if num not in unwanted:
            trimmed_distri.append(num)
    return trimmed_distri


def calc_chi(data, exp_path):
    """
    Using the model file at exp_path as an expected distribution, calculate the chi-square value for a given repetition
    distribution
    :param data: Path of a JSON BitReps file to analyse
    :param exp_path: Path of a JSON BitReps file representing the expected chi-square distribution
    :return: (Chi-square value, observed distribution, expected distribution)
    """
    exp_vals = get_distri(exp_path)
    obs_vals = get_distri(data)
    trimmed_expected = trim_expected(exp_vals)

    fixed_obs = []                          # Only consider observed values which are in the expected distribution
    for val in obs_vals:
        if val in trimmed_expected:
            fixed_obs.append(val)

    count_obs = Counter(fixed_obs)
    count_exp = Counter(trimmed_expected)

    for k, v in count_exp.items():          # If a value in the expected distribution doesn't occur in the observed,
        if k not in count_obs:              # append a zero to the observed to ensure the same number of buckets
            count_obs[k] = 0

    chi_in_obs = []
    chi_in_exp = []

    for k, v in sorted(count_obs.items()):  # Convert distributions to lists of their frequencies for the chi-square
        chi_in_obs.append(v)                # function
    for k, v in sorted(count_exp.items()):
        chi_in_exp.append(v)

    chi = custom_chi(chi_in_obs, chi_in_exp)

    return chi, chi_in_obs, chi_in_exp
