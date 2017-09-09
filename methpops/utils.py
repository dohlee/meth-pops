import os


def colored(string, color='green'):
    if color == 'green':
        return '\033[1;42m' + string + '\033[1;m'
    elif color == 'red':
        return '\033[1;41m' + string + '\033[1;m'
    elif color in ['gray', 'grey']:
        return '\033[1;47m' + string + '\033[1;m'


def parse_genomic_region(string):
    chromosome, coordinates = string.split(':')
    start, end = list(map(int, coordinates.split('-')))
    return chromosome, start, end


def get_sorted_bam_file_name(bam):
    fileName = os.path.basename(bam)
    # Assumes fileName ends with .bam
    return '.'.join(fileName.split('.')[:-1]) + '.sorted.bam'


def get_bam_index_file_name(bam):
    fileName = os.path.basename(bam)
    return fileName + '.bai'


def check_if_sorted_bam_exists(bam, outputDirectory):
    sortedBamName = get_sorted_bam_file_name(bam)

    return os.path.exists(os.path.join(outputDirectory, sortedBamName))


def check_if_bam_index_exists(bam, outputDirectory):
    sortedBamName = get_sorted_bam_file_name(bam)
    indexName = get_bam_index_file_name(sortedBamName)

    return os.path.exists(os.path.join(outputDirectory, indexName))