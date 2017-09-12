import subprocess
import utils
import os

def view(bam, genomicRegion=None):
    if genomicRegion is not None:

        def command(bam, gr):
            return ('samtools view %s %s' % (bam, gr)).split()

        return subprocess.check_output(command(bam, genomicRegion))

    else:

        def command(bam):
            return ('samtools view %s' % bam).split()

        return subprocess.check_output(command(bam))


def index(bam):
    def command(bam):
        return ('samtools index %s' % bam).split()

    subprocess.call(command)


def sort(bam, out, threads):
    def command(bam, out, threads):
        return ('samtools sort %s -o %s --threads %d' % (bam, out, threads)).split()

    subprocess.call(command)


def try_samtools_sort(bam, outputDirectory, threads):
    if utils.check_if_sorted_bam_exists(bam, outputDirectory):
        print('Sorted bam file already exists. Continuing...')
        return

    sortedBamName = utils.get_sorted_bam_file_name(bam)
    outputFile = os.path.join(outputDirectory, sortedBamName)

    samtools.sort(bam, outputFile, threads)

    if not utils.check_if_sorted_bam_exists(bam, outputDirectory):
        raise Exception('Samtools sort did not successfully finish.')


def try_samtools_index(bam, outputDirectory, threads):
    if utils.check_if_bam_index_exists(bam, outputDirectory):
        print('Bam index file already exists. Continuing...')

        if utils.check_if_sorted_bam_exists(bam, outputDirectory):
            print('Sorted bam file already exists. Continuing...')
            return
        else:
            try_samtools_sort(bam, outputDirectory, threads)

    else:
        try_samtools_sort(bam, outputDirectory, threads)
        sortedBamName = utils.get_sorted_bam_file_name(bam)
        sortedBam = os.path.join(outputDirectory, sortedBamName)

        samtools.index(sortedBam)
