from pyfiglet import Figlet
import argparse
from . import helper_functions as hs
from .__version__ import __version__
from .pileup_extractor import FeatureExtractor
from .summary import SummaryCreator
from .poi_analyzer import POIAnalyzer
from .two_sample_extractor import PositionExtractor
from .position_summary import PositionSummary
from .filter import Filter
from .bed_ops import tsv_to_bed, intersect_beds, add_bed_info, merge, difference

def print_figlet(text: str) -> None:
    print(Figlet(font="slant").renderText(text))

def setup_parsers() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="NEET", 
                                     description="""
                                        Nanopore Error pattern Exploration Toolkit (NEET) -- 
                                        Framework for the analysis and exploration of (direct RNA) nanopore sequencing
                                        data. Github: https://github.com/dietvin/neet
                                        """)
    
    parser.add_argument("-v", "--version", action="version", version=f"NEET v{__version__}", help="show version number and exit")

    subparsers = parser.add_subparsers(title="Modules", dest="subcommand")

    # add parser for pileup extractor
    extractor_parser = subparsers.add_parser("extractor", 
                                             help="""
                                                Extract mapped sequencing data from pileup format into an accessible feature table. 
                                                """)

    extractor_parser.add_argument("-i", "--input", type=str, required=True,
                                  help="""
                                  Path to the input pileup file(s). If multiple are available, 
                                  specify paths comma-separated (<pileup1>,<pileup2>,...).
                                  """)
    extractor_parser.add_argument("-o", "--output", type=str, required=True,
                                  help="""
                                  Path to output file(s) (comma-separated if multiple) or directory. 
                                  If filename(s) are given the order corresponds to the input files. 
                                  If a directory is given, the output files are created using the 
                                  basename from an input file with the suffix "_extracted.tsv".
                                  Accordingly the summary file will be created with the suffix
                                  "_extracted_summary.html".
                                  """)
    extractor_parser.add_argument("-r", "--ref", type=str, required=True, 
                                  help="Path to the reference file in fasta format")
    extractor_parser.add_argument("-n", "--num_reads", type=hs.positive_int, required=False,
                                  help="""
                                  Filter by minimum number of reads at a position.  A position will 
                                  be dropped of if the coverage is lower than the given value.
                                  (num_reads = #A + #C + #G + #U + #DEL + #REFSKIP) 
                                  """)
    extractor_parser.add_argument("-m", "--mismatch_rate", type=hs.float_between_zero_and_one, required=False,
                                  help="""
                                  Filter by mismatch rate. A position will be dropped of if mismatch
                                  rate is smaller than the given value. 
                                  (mismatch_rate = #mismatches / (#A + #C + #G + #U + #DEL + #REFSKIP))
                                  """)
    extractor_parser.add_argument("-ma", "--mismatch_rate_alt", type=hs.float_between_zero_and_one, required=False,
                                  help="""
                                  Filter by minimum fraction of mismatched. Relative values measured
                                  on only the number of matched/mismatched reads without deletion and
                                  reference skip rate. 
                                  (mismatch_rate_alt = #mismatches / (#A + #C + #G + #U))
                                  """)
    extractor_parser.add_argument("-d", "--perc_deletion", type=hs.float_between_zero_and_one, required=False,
                                  help="""
                                  Filter by deletion rate. A position will be dropped of if the 
                                  deletion rate is lower than the given value.
                                  """)
    extractor_parser.add_argument("-q", "--mean_quality", type=hs.positive_float, required=False,
                                  help="""
                                  Filter by mean read quality scores. A position will be dropped of if 
                                  the mean quality is lower than the given value.
                                  """)
    extractor_parser.add_argument("-g", "--genomic_region", type=str, required=False,
                                  help="""
                                  Genomic region in "CHR:START-END" format or "CHR" for whole chromosome. 
                                  Specify to only extract information from a specific region.
                                  """)
    extractor_parser.add_argument("--use_multiprocessing", action="store_true", 
                                  help="""
                                  Specify whether to use multiprocessing for processing. Recommended for shorter, 
                                  deeply sequenced data. For low coverage, whole genome data set to false for 
                                  faster runtime.
                                  """)
    extractor_parser.add_argument("-t", "--num_processes", type=int, required=False, default=4,
                                  help="Number of parallel processes to use. Default: 4")
    extractor_parser.add_argument("-nw", "--window_size", type=int, required=False, default=2,
                                  help="""
                                  Number up-/downstream positions to consider for neighbouring error search.
                                  A value w corresponds to a sliding window of site 2*w+1.  Default: w=2
                                  """)
    extractor_parser.add_argument("-nt", "--neighbour_thresh", type=hs.float_between_zero_and_one, required=False, default=0.5,
                                  help="""
                                  Error threshold for neighbourhood search. During the search a neighbour
                                  with a mismatch rate above the given value will be counted as an error.
                                  Default: 0.5
                                  """)
    extractor_parser.add_argument("-b", "--n_bins", type=int, required=False, default=5000,
                                  help="""Number of bins to split the data into when creating the summary report.
                                  This only affects the plots, not the underlying feature table. Used to improve 
                                  performance and clarity of the created plots with many outliers. Note that 
                                  setting the value to a low number can lead to misleading plots. Set to "-1" 
                                  to disable binning. 
                                  Default: 5000
                                  """)
    extractor_parser.add_argument("--plot_alt", action="store_true", 
                                  help="""
                                  Specify whether to use the mismatch_rate or mismatch_rate_alt values for plotting
                                  in the summary report.
                                  """)
    extractor_parser.add_argument("--export_svg", action="store_true", 
                                  help="""
                                  Specify to export the created plots in SVG format. Files will be created in the
                                  output directory.
                                  """)

    # add parser for summary
    summary_parser = subparsers.add_parser("summary", 
                                           help="""
                                            Create a visual overview of a feature table. 
                                            """)
    summary_parser.add_argument("-i", "--input", type=str, required=True,
                                help="""
                                Path to the input TSV file. 
                                """)
    summary_parser.add_argument("-o", "--output", type=str, required=True,
                                help="""
                                Path to output a output directory or file. If a directory is given, the output 
                                files are created using the basename from an input file with the suffix 
                                "_summary.html". If plots --export_svg is set the exported 
                                plots will be stored in the same directory as the HTML report.                                
                                """)
    summary_parser.add_argument("-b", "--n_bins", type=int, required=False, default=5000,
                                help="""Number of bins to split the data into when creating the summary report.
                                This only affects the plots, not the underlying feature table. Used to improve 
                                performance and clarity of the created plots with many outliers. Note that 
                                setting the value to a low number can lead to misleading plots. Set to "-1" 
                                to disable binning. 
                                Default: 5000
                                """)
    summary_parser.add_argument("--plot_alt", action="store_true", 
                                help="""
                                Specify whether to use the mismatch_rate or mismatch_rate_alt values in the 
                                summary report.
                                """)
    summary_parser.add_argument("--export_svg", action="store_true", 
                                help="""
                                Specify to export the created plots in SVG format. Files will be created in the
                                output directory.
                                """)


    # add parser for position-of-interest analyzer
    poi_analyzer_parser = subparsers.add_parser("analyzepoi", 
                                                help="""
                                                    Analyze and explore specific positions of interest
                                                    """)
    
    poi_analyzer_parser.add_argument("-i", "--input", type=str, required=True,
                                     help="""
                                     Path to the input file. Must be of type tsv, as returned by the PileupExtractor.
                                     """)
    poi_analyzer_parser.add_argument("-o", "--output", type=str, required=True,
                                     help="""
                                     Path to the output files or directory. If files are given, one path must
                                     be provided for each bed-category that is given. If a directory is given
                                     a summary file will be created for each given category named
                                     "<category>_<input-basename>_summary.html".
                                     """)
    poi_analyzer_parser.add_argument("-b", "--bed", type=str, required=True,
                                     help="""
                                     Path to the bed file containing information in the fourth column.
                                     """)
    poi_analyzer_parser.add_argument("-r", "--ref", type=str, required=True,
                                     help="""
                                     Path to the reference file. Must be of type fasta.
                                     """)
    poi_analyzer_parser.add_argument("-c", "--bed_categories", type=str, required=True, 
                                     help="""
                                     One or more categories from the bed file to aggregate the data by. 
                                     Must be in the format: "cat1" or "cat1,cat2,cat3"
                                     """)
    poi_analyzer_parser.add_argument("-cc", "--counterparts", type=str, required=False, 
                                     help="""
                                     Canonical base corresponding to each category specified in --bed_categories. 
                                     Same format as --bed_categories flag. Enables comparisons between categories
                                     and counterparts.
                                     """)
    poi_analyzer_parser.add_argument("--update_tsv", action="store_true", 
                                     help="""
                                     If specified, creates an updated feature table (TSV) containing the 
                                     information from the name column in the bed file. Suffix "_w_bed_info.tsv" 
                                     will be added to newly created file.
                                     """)
    poi_analyzer_parser.add_argument("--use_mismatch_rate_alt", action="store_true", 
                                     help="""
                                     If specified use the mismatch_rate_alt values for plotting in the summary report.
                                     """)
    poi_analyzer_parser.add_argument("--export_svg", action="store_true", 
                                     help="""
                                     Specify to export the created plots in SVG format. Files will be created in the
                                     output directory.
                                     """)

    # add parser for two sample extractor
    two_extractor_parser = subparsers.add_parser("twosample", 
                                                 help="""
                                                    Extract positions of differential error rates between two samples.
                                                    """)
    
    two_extractor_parser.add_argument("-i", "--sample1", type=str, required=True,
                                      help="""
                                      Path to the input file(s). If replicates are available, specify paths 
                                      comma-separated (<repl1.tsv>,<repl2.tsv>,...). Must be of type tsv,
                                      as returned by the PileupExtractor.
                                      """)
    two_extractor_parser.add_argument("-bn", "--basename1", type=str, default="sample1",
                                      help="""
                                      Basename of the given sample. Used to create the pileup and extracted 
                                      features files. 
                                      Default: "sample1"
                                      """)
    two_extractor_parser.add_argument("-i2", "--sample2", type=str, required=True,
                                      help="""
                                      Path to the input file(s) from a second sample. If replicates are 
                                      available, specify paths comma-separated (<repl1.tsv>,<repl2.tsv>,...).
                                      Must be of type tsv, as returned by the PileupExtractor. 
                                      """)
    two_extractor_parser.add_argument("-bn2", "--basename2", type=str, default="sample2",
                                      help="""
                                      Basename of the second sample. Used to create the pileup and extracted 
                                      features files. 
                                      Default: "sample2"
                                      """)
    two_extractor_parser.add_argument("-o", "--output", type=str, required=True,
                                      help="""
                                      Path to output a output directory, in which all output files will be stored.
                                      The following files will be created: "<basename>_<basename2>_summary.html", 
                                      "<basename>_<basename2>.bed", "<basename>_excl.bed", "<basename2>_excl.bed"
                                      """)
    two_extractor_parser.add_argument("-r", "--ref", type=str, required=True, 
                                      help="""
                                      Path to the reference file in fasta format
                                      """)
    two_extractor_parser.add_argument("-f", "--error_feature", type=str, default="mismatch_rate_alt", 
                                      help="""
                                      Error feature to use during extraction. Can be one of the following: 
                                      deletion_rate, insertion_rate, mismatch_rate, mismatch_rate_alt. 
                                      Default: "mismatch_rate_alt"
                                      """)
    two_extractor_parser.add_argument("-e", "--error_threshold", type=hs.float_between_zero_and_one, default=0.5, 
                                      help="""
                                      Threshold to identify positions of iterest. Uses the mismatch_rate_alt feature. 
                                      Default: 0.5
                                      """)
    two_extractor_parser.add_argument("-c", "--coverage_threshold", type=hs.positive_int, default=40,
                                      help="""
                                      Minimum coverage of a position to be regarded in the extraction. 
                                      Default: 40
                                      """)
    two_extractor_parser.add_argument("--export_svg", action="store_true", 
                                      help="""
                                      Specify to export the created plots in SVG format. Files will be created in the
                                      output directory.
                                      """)

    # add parser for position summary
    pos_summary_parser = subparsers.add_parser("possummary", 
                                               help="""
                                                Visualize base compositions and error rates at individual positions
                                                between two conditions. 
                                                """)
    
    pos_summary_parser.add_argument("-i", "--sample1", type=str, required=True,
                                    help="""
                                    Path to the input file(s). If replicates are available, 
                                    specify paths comma-separated (<repl1.tsv>,<repl2.tsv>,...).
                                    Must be of type tsv, as returned by the PileupExtractor.
                                    """)
    pos_summary_parser.add_argument("-bn", "--basename1", type=str, default="sample1",
                                    help="""
                                    Basename of the given sample.
                                    Default: "sample1"
                                    """)
    pos_summary_parser.add_argument("-i2", "--sample2", type=str, required=True,
                                    help="""
                                    Path to the input file(s) from a second sample. If replicates are available, 
                                    specify paths comma-separated (<repl1.tsv>,<repl2.tsv>,...).
                                    Must be of type tsv, as returned by the PileupExtractor. 
                                    """)
    pos_summary_parser.add_argument("-bn2", "--basename2", type=str, default="sample2",
                                    help="""
                                    Basename of the second sample.
                                    Default: "sample2"
                                    """)
    pos_summary_parser.add_argument("-b", "--bed", type=str, required=True,
                                    help="""
                                    Path to the input BED file.
                                    """)
    pos_summary_parser.add_argument("-o", "--output_dir", type=str, required=True,
                                    help="""
                                    Path to output directory. The output file is created using the basenames from the input samples
                                    as follows: <basename1>_<basename2>_<basename-bed-file>_summary.html
                                    """)
    pos_summary_parser.add_argument("-n", "--n_surrounding", type=hs.positive_int, required=False, default=2,
                                    help="""
                                    Number of neighbouring positions to consider when plotting a position X. 
                                    In the created plots sites X-n, ..., X, ..., X+n are shown. 
                                    Default: 2
                                    """)
    pos_summary_parser.add_argument("--export_svg", action="store_true", 
                                    help="""
                                    Specify to export the created plots in SVG format. Files will be created in a directory of the 
                                    following format: <output_dir>/<basename1>_<basename2>_<basename-bed-file>
                                    For each position three plots will be created. The files will be named as follows:
                                    <chromosome>_<coordinate>_<plot-type>.svg
                                    """)

    
    # add parser for filtering
    filter_parser = subparsers.add_parser("filter", 
                                          help="""
                                            Filter a feature table on one or more features.
                                            """)
    
    filter_parser.add_argument("-i", "--input", type=str, required=True,
                               help="Path to input TSV file.")
    filter_parser.add_argument("-o", "--output", type=str, required=True,
                               help="""
                               Path to output TSV file or directory. If a directory is given
                               writes the output file with the basename of the input file
                               + "_filtered".
                               """)
    filter_parser.add_argument("-c", "--chromosome", type=str, required=False,
                               help="Filter by given chromosome.")
    filter_parser.add_argument("-s", "--site", type=str, required=False,
                               help="""
                               Filter by given site(s) or range. For single site: "x"; 
                               for multiple sites: "x,y,z,..."; for range: "x-y"
                               """)
    filter_parser.add_argument("-n", "--n_reads", type=str, required=False,
                               help="""
                               Filter by coverage. To filter coverage >= x: "x"; 
                               coverage <= x: "<=x"; coverage == x: "==x"
                               """)
    filter_parser.add_argument("-b", "--base", type=str, required=False,
                               help="""
                               Filter by reference base(s). To filter single base (e.g. A):
                               "A"; multiple bases (e.g. A, C & U): "A,C,U"
                               """)
    filter_parser.add_argument("-m", "--mismatched", action="store_true", required=False,
                               help="Filter mismatched positions.")
    filter_parser.add_argument("-mt", "--mismatch_types", type=str, required=False,
                               help="""
                               Filter one or more specific types of mismatches. 
                               E.g: filter A-to-T mismatches --> "A-T"; filter A-to-T and 
                               C-to-T mismatches --> "A-T,C-T"
                               """)
    filter_parser.add_argument("-p", "--mismatch_rate", type=str, required=False,
                               help="""
                               Filter by mismatch rate. To filter mismatch_rate >= x:
                               "x"; mismatch_rate <= x: "<=x"
                               """)
    filter_parser.add_argument("-pa", "--mismatch_rate_alt", type=str, required=False,
                               help="""
                               Filter by percent of mismatched reads using the alternative measure. 
                               To filter mismatch_rate_alt >= x: "x"; mismatch_rate_alt <= x: "<=x"
                               """)
    filter_parser.add_argument("-d", "--deletion_rate", type=str, required=False,
                               help="""
                               Filter by percent of deleted reads. To filter deletion_rate >= x:
                               "x"; deletion_rate <= x: "<=x"
                               """)
    filter_parser.add_argument("-pi", "--insertion_rate", type=str, required=False,
                               help="""
                               Filter by percent of inserted reads. To filter insertion_rate >= x: 
                               "x"; insertion_rate <= x: "<=x"
                               """)
    filter_parser.add_argument("-pr", "--refskip_rate", type=str, required=False,
                               help="""
                               Filter by percent of reads with reference skip. To filter refskip_rate >= x:
                               "x"; refskip_rate <= x: "<=x"
                               """)
    filter_parser.add_argument("-f", "--motif", type=str, required=False,
                               help="Filter by motif around position.")
    filter_parser.add_argument("-q", "--q_score", type=str, required=False,
                               help="""
                               Filter by mean quality. To filter q_mean >= x: "x"; q_mean <= x: "<=x"
                               """)
    filter_parser.add_argument("-bi", "--filter_bed", type=str, required=False,
                               help="""
                               Path to a bed file. The TSV file will be filtered by the positions from the bed file,
                               keeping only positions that are found in the bed file.
                               """)
    filter_parser.add_argument("-be", "--exclude_bed", type=str, required=False,
                               help="""
                               Path to a bed file. The TSV file will be filtered by the positions from the bed file,
                               keeping only positions that are NOT found in the bed file.
                               """)
    
    # add parser for bed ops
    bedops_parser = subparsers.add_parser("bedops", 
                                          help="""
                                            Perform minor tasks related to bed files. 
                                            """)
    subsubparsers = bedops_parser.add_subparsers(title="Bed-Ops commands", dest="subsubcommand",
                                                 help="""
                                                    Perform minor tasks related to bed files.
                                                    """)

    tsv2bed_parser = subsubparsers.add_parser("tsv2bed", help="Transform feature table into BED format")
    tsv2bed_parser.add_argument("-i", "--input", type=str, required=True,
                                help="Path to the input file")
    tsv2bed_parser.add_argument("-o", "--output", type=str, required=True,
                                help="Path to the output file.")
    
    intersect_parser = subsubparsers.add_parser("intersect", help="Extracts shared and exclusive positions from two bed files.")
    intersect_parser.add_argument("-a", "--file_a", type=str, required=True,
                                  help="First BED file")
    intersect_parser.add_argument("-b", "--file_b", type=str, required=True,
                                  help="Second BED file")
    intersect_parser.add_argument("-o", "--output", type=str, required=True,
                                  help="Path to the output file.")
    intersect_parser.add_argument("--label_a", type=str, required=True,
                                  help="Label given to file a")
    intersect_parser.add_argument("--label_b", type=str, required=True,
                                  help="Label given to file b")
    
    add_bed_parser = subsubparsers.add_parser("addinfo", help="Add the name column from a BED file to the respective line in the TSV file.")
    add_bed_parser.add_argument("-i", "--input", type=str, required=True,
                                help="Feature table output from PileupExtractor/NeighbourhoodSearcher")
    add_bed_parser.add_argument("-o", "--output", type=str, required=True,
                                help="Path to the output file")
    add_bed_parser.add_argument("-b", "--bed", type=str, required=True,
                                help="BED file containing additional information in the 'name' column")
    
    merge_parser = subsubparsers.add_parser("merge", help="Merges the position contained in multiple bed files.")
    merge_parser.add_argument("-i", "--input", type=str, required=True,
                              help="Paths to multiple bed files, comma separated. For example: '/path/to/file1.bed,/path/to/file2.bed,...'")
    merge_parser.add_argument("-o", "--output", type=str, required=True,
                              help="Path to the output file")

    difference_parser = subsubparsers.add_parser("difference", help="Extract the positons present in file 1 and not present in file 2.")
    difference_parser.add_argument("-i1", "--input1", type=str, required=True,
                                   help="Path to bed file 1")
    difference_parser.add_argument("-i2", "--input2", type=str, required=True,
                                   help="Path to bed file 2")
    difference_parser.add_argument("-o", "--output", type=str, required=True,
                                   help="Path to the output file")


    return parser, bedops_parser

def main():

    parser, bedops_parser = setup_parsers()
    args = parser.parse_args()

    if args.subcommand == "extractor":
        print_figlet("NEET - Pileup Extractor")
        feature_extractor = FeatureExtractor(in_paths=args.input, 
                                             out_paths=args.output, 
                                             ref_path=args.ref,
                                             num_reads=args.num_reads, 
                                             mismatch_rate=args.mismatch_rate,
                                             mismatch_rate_alt=args.mismatch_rate_alt,
                                             perc_deletion=args.perc_deletion,
                                             mean_quality=args.mean_quality,
                                             genomic_region=args.genomic_region,
                                             use_multiprocessing=args.use_multiprocessing,
                                             num_processes=args.num_processes,
                                             window_size=args.window_size,
                                             neighbour_error_threshold=args.neighbour_thresh,
                                             n_bins_summary = args.n_bins,
                                             use_alt_summary=args.plot_alt,
                                             export_svg_summary=args.export_svg)
        feature_extractor.main()

    elif args.subcommand == "summary":
        print_figlet("NEET - Summary")
        summary = SummaryCreator(in_path=args.input, 
                                 out_path=args.output, 
                                 n_bins=args.n_bins, 
                                 use_perc_mismatch_alt=args.plot_alt, 
                                 export_svg=args.export_svg,)
        summary.main()

    elif args.subcommand == "analyzepoi":
        print_figlet("NEET - Position-of-Interest Analyzer")
        poi_analyzer = POIAnalyzer(in_path=args.input,
                                   out_path=args.output,
                                   bed_path=args.bed,
                                   ref_path=args.ref,
                                   categories=args.bed_categories,
                                   canonical_counterpart=args.counterparts,
                                   output_tsv=args.update_tsv,
                                   use_perc_mismatch_alt=args.use_mismatch_rate_alt,
                                   export_svg=args.export_svg)
        poi_analyzer.main()

    elif args.subcommand == "twosample":
        print_figlet("NEET - Two-Sample Extractor")
        posextr = PositionExtractor(in_paths_a=args.sample1, 
                                    in_paths_b=args.sample2, 
                                    out_dir=args.output, 
                                    ref_path=args.ref, 
                                    label_a=args.basename1, 
                                    label_b=args.basename2, 
                                    error_feature=args.error_feature,
                                    error_threshold=args.error_threshold, 
                                    coverage_threshold=args.coverage_threshold,
                                    export_svg=args.export_svg)
        posextr.main()

    elif args.subcommand == "possummary":
        print_figlet("NEET - Position Summary")
        pos_summary = PositionSummary(paths_a = args.sample1,
                                      paths_b = args.sample2,
                                      basename_a = args.basename1,
                                      basename_b = args.basename2,
                                      bed_path = args.bed,
                                      out_path = args.output_dir,
                                      nb_size = args.n_surrounding,
                                      export_svg=args.export_svg)
        pos_summary.main() 
    
    elif args.subcommand == "filter":
        print_figlet("NEET - Filter")
        filter = Filter(input_path=args.input,
                        output_path=args.output,
                        chrom=args.chromosome,
                        site=args.site,
                        n_reads=args.n_reads,
                        base=args.base,
                        mismatched=args.mismatched,
                        mismatch_types=args.mismatch_types,
                        perc_mismatched=args.percent_mismatched,
                        perc_mismatched_alt=args.mismatch_rate_alt,
                        perc_deletion=args.percent_deletion,
                        perc_insertion=args.percent_insertion,
                        perc_refskip=args.percent_refskip,
                        motif=args.motif,
                        q_score=args.q_score,
                        bed_include=args.filter_bed,
                        bed_exclude=args.exclude_bed)
        filter.main()
    
    elif args.subcommand == "bedops":
        print_figlet("NEET - Bed Ops")
        if args.subsubcommand == "tsv2bed":
            tsv_to_bed(args.input, args.output)
        elif args.subsubcommand == "intersect":
            intersect_beds(args.file_a, args.file_b, args.output, args.label_a, args.label_b)
        elif args.subsubcommand == "addinfo":
            add_bed_info(args.input, args.bed, args.output)
        elif args.subsubcommand == "merge":
            merge(args.input, args.output)
        elif args.subsubcommand == "difference":
            difference(bed1=args.input1, bed2=args.input1, out=args.output)
        else:
            bedops_parser.print_help()
    else:
        print_figlet("NEET")
        parser.print_help()

if __name__=="__main__":
    main()


