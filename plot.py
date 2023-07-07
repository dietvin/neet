import argparse
from helper_functions import float_between_zero_and_one

from plot_coverage_track import CoverageTrackPlotter
from plot_compositions import CompositionPlotter

if __name__=="__main__":
    parser = argparse.ArgumentParser(prog="CoverageTrackPlotter",
                                     description="Provides different plotting functionalities.")
    subparsers = parser.add_subparsers(dest="subcommand")
    coverage_track_parser = subparsers.add_parser("coverage_track")
    coverage_track_parser.add_argument("-i", "--input", type=str, required=True,
                                       help="Path to the input file")
    coverage_track_parser.add_argument("-o", "--output", type=str, required=False,
                                       help="Path to the output file")
    coverage_track_parser.add_argument("-r", "--reference", type=str, required=True,
                                       help="Path to the reference fasta file")
    coverage_track_parser.add_argument('-t', '--error_threshold', type=float_between_zero_and_one, required=True,
                        help="Percentage of error threshold for coloring positions [0,1]")
    coverage_track_parser.add_argument('-c', '--chromosome', type=str, required=True,
                        help="Chromosome to be plotted")
    
    args = parser.parse_args()

    if args.subcommand == "coverage_track":
        plotter = CoverageTrackPlotter(in_path=args.input,
                          ref_path=args.reference,
                          threshold=args.error_threshold,
                          chromosome=args.chromosome,
                          out_path=args.output)
        plotter.create_plot()
    else:
        parser.print_help()