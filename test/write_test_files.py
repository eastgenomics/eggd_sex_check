#!/usr/bin/env python
"""
Script to generate test data files for testing get_mapped_reads function.
Each file represents different scenarios:
- Complete data
- Without chromosome 1
- Without chromosome Y
- Without chromosome 1 but with chr11
"""

import os

# Test data for correct idxstat output
CORRECT_DATA = """1\t100\t100\t0
2\t90\t80\t10
X\t80\t70\t10
Y\t70\t60\t10
Z\t50\t0\t50
"""

# Test data without chromosome 1
WITHOUT_CHR1_DATA = """2\t90\t80\t10
X\t80\t70\t10
Y\t70\t60\t10
Z\t50\t0\t50
M\t40\t0\t40
"""

# Test data without chromosome Y
WITHOUT_CHRY_DATA = """1\t100\t100\t0
2\t90\t80\t10
X\t80\t70\t10
Z\t50\t0\t50
M\t40\t0\t40
"""

# Test data without chromosome 1 but with chr11
WITH_CHR11_DATA = """11\t100\t100\t0
2\t90\t80\t10
X\t80\t70\t10
Y\t70\t60\t10
Z\t50\t0\t50
"""


def write_test_data(directory, filename, data):
    """
    Write test data to a file.

    Args:
        directory (str): Directory to write the file.
        filename (str): Name of the file to write.
        data (str): Content to write to the file.
    """
    os.makedirs(directory, exist_ok=True)
    filepath = os.path.join(directory, filename)
    with open(filepath, "w", encoding="utf-8") as file:
        file.write(data)


TEST_CASES = [
    ("correct_data.tsv", CORRECT_DATA),
    ("without_chr1.tsv", WITHOUT_CHR1_DATA),
    ("without_chry.tsv", WITHOUT_CHRY_DATA),
    ("with_chr11.tsv", WITH_CHR11_DATA),
]

TEST_DATA_DIR = "test_data"

if __name__ == "__main__":
    for filename, data in TEST_CASES:
        write_test_data(TEST_DATA_DIR, filename, data)
