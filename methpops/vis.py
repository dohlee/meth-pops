#!/packages/Anaconda2-4.1.1/bin/python2.7

import os
import argparse
import subprocess
import matplotlib
import pandas as pd
matplotlib.use('agg')
import matplotlib.pyplot as plt

from read import generate_read_from_bam_line

def parse_arguments():
	parser = argparse.ArgumentParser()
	parser.add_argument("-i", "--input", required=True, help="Input file")
	parser.add_argument("-o", "--output", required=True, help="Output directory")
	parser.add_argument("-p", "--threads", type=int, default=1, help="Number of threads to be used")

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
		print('Sorted bam file already exits. Continuing...')
		return

	sortedFileName = get_sorted_bam_file_name(inputFile)
	subprocess.call(['samtools', 'sort', inputFile, '-o', os.path.join(outputDirectory, sortedFileName), '--threads', str(threads)])

	if not check_if_sorted_bam_exists(inputFile, outputDirectory):
		raise Exception('Samtools did not successfully finish.')

def try_samtools_index(inputFile, outputDirectory, threads):
	if check_if_bam_index_exists(inputFile, outputDirectory):
		print('Bam index file already exits. Continuing...')
		return

	try_samtools_sort(inputFile, outputDirectory, threads)
	sortedFileName = get_sorted_bam_file_name(inputFile)

	subprocess.call(['samtools', 'index', os.path.join(outputDirectory, sortedFileName)])

def get_reads_aligned_to_region(inputFile, chromosome, start, end):
	output = subprocess.check_output(['samtools', 'view', inputFile, '%d:%d-%d' % (chromosome, start, end)])
	reads = [generate_read_from_bam_line(line) for line in output.splitlines()]

	return reads

if __name__ == '__main__':
	args = parse_arguments()
	try_samtools_index(args.input, args.output, args.threads)
	sortedFileName = get_sorted_bam_file_name(args.input)

	df = pd.read_csv('cpgIslandExt.txt', sep='\t', header=None)
	for i, line in df.iterrows():
		start, end = line[2], line[3]
		print(start, end)
		reads = get_reads_aligned_to_region(os.path.join(args.output, sortedFileName), 1, start, end)
		
		for read in reads:
			print(read.get_CpGs())

		if i == 300:
			break

	
