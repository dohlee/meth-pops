#!/packages/Anaconda2-4.1.1/bin/python2.7

import os
import argparse
import subprocess
import matplotlib
import pandas as pd
matplotlib.use('agg')
import matplotlib.pyplot as plt

from read import generate_read_from_bam_line
from utils import colored, parse_genomic_region

def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", required=True, help="Input file")
    parser.add_argument("-o", "--output", required=True, help="Output directory")
    parser.add_argument("-p", "--threads", type=int, default=1, help="Number of threads to be used")
    parser.add_argument("-r", "--region", type=str, default=None, help="Genomic region to visualize")

    return parser.parse_args()

def get_sorted_bam_file_name(inputFile):
    fileName = os.path.basename(inputFile)
    # assumes fileName ends with .bam
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

    sortedFileName = get_sorted_bam_file_name(inputFile)
    subprocess.call(['samtools', 'sort', inputFile, '-o', os.path.join(outputDirectory, sortedFileName), '--threads', str(threads)])

    if not check_if_sorted_bam_exists(inputFile, outputDirectory):
        raise Exception('Samtools did not successfully finish.')

def try_samtools_index(inputFile, outputDirectory, threads):
    if check_if_bam_index_exists(inputFile, outputDirectory):
        print('Bam index file already . Continuing...')
        return

    try_samtools_sort(inputFile, outputDirectory, threads)
    sortedFileName = get_sorted_bam_file_name(inputFile)

    subprocess.call(['samtools', 'index', os.path.join(outputDirectory, sortedFileName)])

def get_reads_aligned_to_region(inputFile, chromosome, start, end):
    output = subprocess.check_output(['samtools', 'view', inputFile, '%s:%d-%d' % (str(chromosome), start, end)])
    reads = [generate_read_from_bam_line(line) for line in output.splitlines()]

    return reads

def pretty_print(reads):
    if len(reads) > 75:
        print("There are many reads to show. You may think of creating a metylation lollipop plot instead.")
        return

    positions = sorted(list(set(pos for read in reads for pos, meth in read.get_CpGs())))

    for read in reads:
        CpGs = read.get_CpGs()
        if len(CpGs) == 0:
            continue

        i = positions.index(CpGs[0][0])
        
        print('-' * i + read.get_CpG_colored_string() + '-' * (len(positions) - len(CpGs) - i))

def pretty_print_region(sortedBamFile, chromosome, start, end):
    reads = get_reads_aligned_to_region(sortedBamFile, chromosome, start, end)
    pretty_print(reads)

def print_region(sortedBamFile, chromosome, start, end):
    reads = get_reads_aligned_to_region(sortedBamFile, chromosome, start, end)
    startIndex = reads[0].position

    for read in reads:
        print(' ' * (read.position - startIndex) + read.get_full_CpG_colored_string())

if __name__ == '__main__':
    args = parse_arguments()
    try_samtools_index(args.input, args.output, args.threads)
    sortedFileName = get_sorted_bam_file_name(args.input)
    sortedBamFile = os.path.join(args.output, sortedFileName)

    # df = pd.read_csv('cpgIslandExt.txt', sep='\t', header=None)
    # for i, line in df.iterrows():
    #     start, end = line[2], line[3]
    #     reads = get_reads_aligned_to_region(os.path.join(args.output, sortedFileName), 1, start, end)
        # print(start, end)
        # if len(reads) > 5:
        #     startPosition = min([r.position for r in reads])
        #     print(start, end)
        #     for read in reads:
        #         print(read.get_CpG_colored_string())
        #         print(read.reverse, read.flag)
        #         print(read.get_CpGs())

        # if len(reads) > 10:
        #     print(start, end)
        #     pretty_print(reads)

    if args.region == None:
        print_region(sortedBamFile, 8, 97643825, 97644200)
    else:
        chromosome, start, end = parse_genomic_region(args.region)
        print_region(sortedBamFile, chromosome, start, end)

    # reads = get_reads_aligned_to_region(os.path.join(args.output, sortedFileName), 1, 3188392, 3188689)
    # pretty_print(reads)
    # for read in reads:
    #     print(read.get_CpG_colored_string())
    #     print(read.reverse, read.flag)
    #     print(read.get_CpGs())



            
