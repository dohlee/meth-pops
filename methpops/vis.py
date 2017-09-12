import os
import argparse
import subprocess
import numpy as np
import matplotlib

import samtools_wrapper as samtools
import lollipop
from read import generate_read_from_bam_line
from utils import parse_genomic_region
matplotlib.use('agg')


def get_reads_aligned_to_region(inputFile, genomicRegion):
    output = samtools.view(inputFile, genomicRegion)
    reads = [generate_read_from_bam_line(line) for line in output.splitlines()]

    return reads


def get_sorted_list_of_CpG_positions(reads):
    positionSet = set(pos for read in reads for pos, meth in read.get_CpGs())
    return sorted(list(positionSet))


def pretty_print(reads):
    if len(reads) > 75:
        print("There are many reads to show. You may think of creating a methylation lollipop plot instead.")
        return

    if len(reads) == 0:
        print("There aren't any reads to show.")
        return

    positions = get_sorted_list_of_CpG_positions(reads)

    for read in reads:
        CpGs = read.get_CpGs()
        if len(CpGs) == 0:
            continue

        i = positions.index(CpGs[0][0])

        print('-' * i + read.get_CpG_colored_string() +
              '-' * (len(positions) - len(CpGs) - i))


def pretty_print_region(sortedBamFile, genomicRegion):
    reads = get_reads_aligned_to_region(sortedBamFile, genomicRegion)
    pretty_print(reads)


def print_region(sortedBamFile, genomicRegion):
    reads = get_reads_aligned_to_region(sortedBamFile, genomicRegion)
    if len(reads) == 0:
        print("There aren't any reads to show.")
        return

    startIndex = reads[0].position

    for read in reads:
        print(' ' * (read.position - startIndex) +
              read.get_full_CpG_colored_string())


def CpG_matrix(reads, positions):
    matrix = np.full([len(reads), len(positions)], -1, dtype=int)

    for i, read in enumerate(reads):
        CpGs = read.get_CpGs()
        for coordinate, meth in CpGs:
            matrix[i, positions.index(coordinate)] = meth

    return matrix


def save_lollipop_plot(sortedBamFile, genomicRegion, output=None):
    reads = get_reads_aligned_to_region(sortedBamFile, genomicRegion)
    positions = get_sorted_list_of_CpG_positions(reads)

    matrix = CpG_matrix(reads, positions)
    lollipop.save(matrix, genomicRegion, output)
