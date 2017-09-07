#!/packages/Anaconda2-4.1.1/bin/python2.7

import os
import argparse
import subprocess

from read import generate_read_from_bam_line
from utils import parse_genomic_region


def parse_arguments():
    """TODO: Docstring for parse_arguments.

    :returns: TODO
     """

    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", required=True, help="Input file")
    parser.add_argument("-o", "--output", required=True,
                        help="Output directory")
    parser.add_argument("-p", "--threads", type=int,
                        default=1, help="Number of threads to be used")
    parser.add_argument("-r", "--region", type=str,
                        default=None, help="Genomic region to visualize")

    return parser.parse_args()


def get_sorted_bam_file_name(inputFile):
    fileName = os.path.basename(inputFile)
    # Assumes fileName ends with .bam
    return '.'.join(fileName.split('.')[:-1]) + '.sorted.bam'


def get_bam_index_file_name(inputFile):
    fileName = os.path.basename(inputFile)
    return fileName + '.bai'


def check_if_sorted_bam_exists(inputFile, outputDirectory):
    sortedFileName = get_sorted_bam_file_name(inputFile)

    return os.path.exists(os.path.join(outputDirectory, sortedFileName))


def check_if_bam_index_exists(inputFile, outputDirectory):
    sortedFileName = get_sorted_bam_file_name(inputFile)
    indexName = get_bam_index_file_name(sortedFileName)

    return os.path.exists(os.path.join(outputDirectory, indexName))


def try_samtools_sort(inputFile, outputDirectory, threads):
    if check_if_sorted_bam_exists(inputFile, outputDirectory):
        print('Sorted bam file already exists. Continuing...')
        return

    command = lambda i, o, t: 'samtools sort %s -o %s --threads %s' % (i, o, t)

    sortedFileName = get_sorted_bam_file_name(inputFile)
    outputFile = os.path.join(outputDirectory, sortedFileName)

    subprocess.call(command(inputFile, outputFile, str(threads)).split())

    if not check_if_sorted_bam_exists(inputFile, outputDirectory):
        raise Exception('Samtools sort did not successfully finish.')


def try_samtools_index(inputFile, outputDirectory, threads):
    if check_if_bam_index_exists(inputFile, outputDirectory):
        print('Bam index file already . Continuing...')
        return

    command = lambda i: 'samtools index %s' % i

    try_samtools_sort(inputFile, outputDirectory, threads)
    sortedFileName = get_sorted_bam_file_name(inputFile)
    inputFile = os.path.join(outputDirectory, sortedFileName)

    subprocess.call(command(inputFile).split())


def get_reads_aligned_to_region(inputFile, genomicRegion):
    # gr : string indicating genomic region
    command = lambda i, gr: 'samtools view %s %s' % (i, gr)

    output = subprocess.check_output(command(inputFile, genomicRegion).split())
    reads = [generate_read_from_bam_line(line) for line in output.splitlines()]

    return reads


def get_sorted_list_of_CpG_positions(reads):
    positionSet = set(pos for read in reads for pos, meth in read.get_CpGs())
    return sorted(list(positionSet))


def pretty_print(reads):
    if len(reads) > 75:
        print("There are many reads to show. You may think of creating a methylation lollipop plot instead.")
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
    startIndex = reads[0].position

    for read in reads:
        print(' ' * (read.position - startIndex) +
              read.get_full_CpG_colored_string())

if __name__ == '__main__':
    args = parse_arguments()
    try_samtools_index(args.input, args.output, args.threads)
    sortedFileName = get_sorted_bam_file_name(args.input)
    sortedBamFile = os.path.join(args.output, sortedFileName)

    if args.region is None:
        print_region(sortedBamFile, 8, 97643825, 97644200)
    else:
        chromosome, start, end = parse_genomic_region(args.region)
        print_region(sortedBamFile, chromosome, start, end)
