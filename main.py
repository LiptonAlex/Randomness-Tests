import math
import numpy
from scipy import special as spc
from random import getrandbits


def monobit():
    rng_output = [getrandbits(1) for i in range(20000)]
    bin_data = "".join(str(x) for x in rng_output[:20000])
    print("Generated binary for Monobit: " + bin_data)
    count = 0

    for char in bin_data:
        if char == "0":
            count -= 1
        else:
            count += 1

    sobs = count / math.sqrt(len(bin_data))
    p_val = spc.erfc(math.fabs(sobs) / math.sqrt(2))

    if p_val > 0.01:
        print(">> Monobit Test: PASSED (P val > 0.01) -->", p_val)
    else:
        print(">> Monobit Test: NOT PASSED (P val < 0.01)  -->  ", p_val)
    return p_val


def max_serial():
    rng_output = [getrandbits(1) for i in range(20000)]
    bin_data = "".join(str(x) for x in rng_output[:20000])
    print("Generated binary for Max Frequency Test: " + bin_data)
    current_run = 0
    longest_run = 0
    threshold = 36

    for i in bin_data:
        if i == "1":
            current_run += 1
            if current_run > longest_run:
                longest_run = current_run
        else:
            current_run = 0

    return longest_run <= threshold


def longest_run():
    rng_output = [getrandbits(1) for i in range(20000)]
    bin_data = "".join(str(x) for x in rng_output[:20000])
    print("Generated binary for Longest Run: " + bin_data)
    if len(bin_data) < 128:
        print("\t", "Longest Run >> Not enough data to run test!")
        return -1.0
    elif len(bin_data) < 6272:
        k, m = 3, 8
        v_values = [1, 2, 3, 4]
        pik_values = [0.21484375, 0.3671875, 0.23046875, 0.1875]
    elif len(bin_data) < 75000:
        k, m = 5, 128
        v_values = [4, 5, 6, 7, 8, 9]
        pik_values = [0.1174035788, 0.242955959, 0.249363483, 0.17517706, 0.102701071, 0.112398847]
    else:
        k, m = 6, 10000
        v_values = [10, 11, 12, 13, 14, 15, 16]
        pik_values = [0.0882, 0.2092, 0.2483, 0.1933, 0.1208, 0.0675, 0.0727]

    num_blocks = math.floor(len(bin_data) / m)
    frequencies = numpy.zeros(k + 1)
    block_start, block_end = 0, m
    for i in range(num_blocks):
        # Slice the binary string into a block
        block_data = bin_data[block_start:block_end]
        # Keep track of the number of ones
        max_run_count, run_count = 0, 0
        for j in range(0, m):
            if block_data[j] == '1':
                run_count += 1
                max_run_count = max(max_run_count, run_count)
            else:
                max_run_count = max(max_run_count, run_count)
                run_count = 0
        max_run_count = max(max_run_count, run_count)
        if max_run_count < v_values[0]:
            frequencies[0] += 1
        for j in range(k):
            if max_run_count == v_values[j]:
                frequencies[j] += 1
        if max_run_count > v_values[k - 1]:
            frequencies[k] += 1
        block_start += m
        block_end += m
    chi_squared = 0
    for i in range(len(frequencies)):
        chi_squared += (pow(frequencies[i] - (num_blocks * pik_values[i]), 2.0)) / (num_blocks * pik_values[i])
    p_val = spc.gammaincc(float(k / 2), float(chi_squared / 2))
    print(">> Longest Run Test: " + str(p_val))
    return p_val


monobit()
longest_run()
test_result = max_serial()
if test_result:
    print(">> Max Serial Test passed.")
else:
    print(">> Max Serial Test failed.")