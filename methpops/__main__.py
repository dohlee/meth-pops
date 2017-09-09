import os
import vis
import argparse
import utils
import samtools_wrapper as samtools
EXAMPLE_REGION = '1:6466114-6466361'


def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", required=True, help="Input file.")
    parser.add_argument("-o", "--output", help="Output lollipop plot.")
    parser.add_argument("-d", "--tempdir",
                        help="Temporary directory in which the .sorted.bam file and bam index file will be saved.")
    parser.add_argument("-t", "--threads", type=int,
                        default=1, help="Number of threads to be used.")
    parser.add_argument("-c", "--commandline", action='store_true',
                        default=False, help="Show methylation status at commandline.")
    parser.add_argument("--pretty", action='store_true',
                        default=False, help="Show methylation status at commandline in prettier way.")
    parser.add_argument("-l", "--lollipop", action='store_true',
                        default=True, help="Save lollipop plot. Default: True")
    parser.add_argument("-r", "--region", type=str,
                        default=None, help="Genomic region to visualize.")

    args = parser.parse_args()
    if args.tempdir is None:
        args.tempdir = os.path.dirname(args.input)

    return args


def main():
    args = parse_arguments()
    # At first try to index bam file.
    samtools.try_samtools_index(args.input, args.tempdir, args.threads)
    sortedFileName = utils.get_sorted_bam_file_name(args.input)
    sortedBamFile = os.path.join(args.tempdir, sortedFileName)

    if args.region is None:
        args.region = EXAMPLE_REGION

    # Generate lollipop plot.
    if args.lollipop:
        vis.save_lollipop_plot(sortedBamFile, args.region, args.output)

    # Show commandline output.
    if args.pretty:
        vis.pretty_print_region(sortedBamFile, args.region)
    elif args.commandline:
        vis.print_region(sortedBamFile, args.region)


if __name__ == '__main__':
    main()
