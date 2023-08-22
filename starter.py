import subprocess, argparse, os, warnings, datetime
from typing import List, Dict
from pileup_extractor import FeatureExtractor
from stat_compare import StatComparer
import helper_functions as hs


class Processor:
    in_path1: List[str]
    in_path2: List[str]|None
    basename1: str
    basename2: str
    two_samples: bool
    out_path: str
    ref_path: str
    optional_args: Dict[str, int|float|bool|None]
    stat_comp_args: Dict[str, float|bool]

    def __init__(self, in_path1, basename1, out_path, ref_path, in_path2, basename2,
                 num_reads,
                 perc_mismatched,
                 perc_deletion,
                 mean_quality,
                 genomic_region,
                 use_multiprocessing,
                 num_processes,
                 coverage_alt,
                 window_size,
                 neighbour_thresh,
                 n_bins,
                 alpha,
                 no_tsv) -> None:
        hs.print_update("Start initialization.")
        self.process_paths(in_path1, in_path2, out_path, ref_path)
        self.basename1 = basename1
        self.basename2 = basename2

        self.optional_args = {"num_reads": num_reads,
                              "perc_mismatch": perc_mismatched,
                              "perc_deletion": perc_deletion,
                              "mean_quality": mean_quality,
                              "genomic_region": genomic_region,
                              "use_multiprocessing": use_multiprocessing,
                              "num_processes": num_processes,
                              "use_alt_coverage": coverage_alt,
                              "window_size": window_size,
                              "neighbour_error_threshold": neighbour_thresh,
                              "n_bins_summary": n_bins}
        self.stat_comp_args = {"alpha": alpha, "write_tsv": (not no_tsv)}

        if self.two_samples:
            message = f"Found {len(self.in_path1)} bam files for the first sample and {len(self.in_path2)} for the second sample."
        else:
            message = f"Found {len(self.in_path1)} bam files from one sample."
        hs.print_update("Initialization sucessful. " + message)


    def process_paths(self, in1: str, in2: str|None, out: str, ref_path: str) -> None:
        self.in_path1 = self.process_in(in1)
        if in2:
            self.in_path2 = self.process_in(in2)
            self.two_samples = True
        else:
            self.two_samples = False

        hs.check_path(ref_path, extensions=[".fa", ".fasta", ".fn"])
        self.ref_path = ref_path

        if not out.endswith("/"): out+="/"
        self.check_out(out)
        self.out_path = out

    def process_in(self, in_paths: str) -> List[str]:
        in_list = in_paths.split(",")
        for path in in_list:
            self.check_path(path)            
        return in_list

    def check_out(self, out: str): 
        if not os.path.isdir(out):
            warnings.warn(f"Directory '{out}' does not exist. Creating directory.")
            try:
                os.makedirs(out)
            except:
                raise Exception(f"Directory '{out}' was not found and could not be created.")

    def check_path(self, path: str) -> None:
        """
        Check if the specified file path exists and has the expected file extension.

        This function verifies whether the file specified by the given path exists and has a valid extension.
        If the file does not exist, it raises a FileNotFoundError with a detailed error message.
        If the file extension does not match any of the expected extensions, it raises an Exception.

        Parameters:
            path (str): The file path to be checked.

        Raises:
            FileNotFoundError: If the specified file path does not exist.
            Exception: If the file format is not correct.
        """
        if not os.path.exists(path): # does file exist?
            raise FileNotFoundError(f"Input file not found. File '{path}' does not exist.")
        file_type = os.path.splitext(path)[1]
        if file_type != ".bam":
            raise Exception(f"Found file extension {file_type}. File extension to must be .bam.")


    def main(self):
        pileup_path = self.out_path + self.basename1 + ".pileup"
        feature_extract_path = self.out_path + self.basename1 + "_extracted.tsv"
        
        self.create_pileup(self.in_path1, pileup_path)
        if self.two_samples:
            pileup_path2 = self.out_path + self.basename2 + ".pileup"
            self.create_pileup(self.in_path2, pileup_path2)

            pileup_path += ","+pileup_path2 
            feature_extract_path += ","+ self.out_path + self.basename2 + "_extracted.tsv"
        
        self.run_feature_extractor(pileup_path, feature_extract_path)

        if self.two_samples:
            samples = feature_extract_path.split(",")
            self.run_stat_comparer(samples[0], samples[1])

    def create_pileup(self, paths: List[str], out: str) -> None:
        hs.print_update(f"Creating pileup file at {out} from file(s): {paths}")
        subprocess.run(["samtools", "mpileup", "-ABQ0", "-d", "0", "-f", self.ref_path, "-o", out] + paths)

    def run_feature_extractor(self, pileup_path: str, out_path: str) -> None:
        hs.print_update("Starting feature extractor.")
        args = self.optional_args
        extractor = FeatureExtractor(in_paths = pileup_path, 
                                     out_paths = out_path, 
                                     ref_path = self.ref_path,
                                     **args)
        extractor.process_files()

    def run_stat_comparer(self, sample1: str, sample2: str) -> None:
        out = f"{self.out_path}{self.basename1}_{self.basename2}_comp.tsv" 
        hs.print_update(f"Starting statistical comparison of files {sample1} and {sample2}. Writing to {out}")
        stat_comp = StatComparer(sample1, sample2, out, **self.stat_comp_args)
        stat_comp.main()

def setup_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="Neet", description="Extract different characteristics from given bam file(s).")
    ### Input / output arguments ###
    parser.add_argument('-i', '--sample1', type=str, required=True,
                        help="""
                            Path to the input file(s). If replicates are available, specify paths comma-separated (<repl.1>,<repl.2>,...).
                            Can be of type .bam or .pileup. If bam files are given, samtools mpileup is executed first.
                            """)
    parser.add_argument('-bn', '--basename', type=str, required=False, default="sample1",
                    help="""
                        Basename of the given sample. Used to create the pileup and extracted features files. Default: 'sample1'
                        """)
    parser.add_argument('-i2', '--sample2', type=str, required=False,
                        help="""
                            Path to the input file(s) from a second sample. If replicates are available, specify paths comma-separated (<repl.1>,<repl.2>,...).
                            Can be of type .bam or .pileup. If provided the extracted features will be statistically compared to the ones from the first sample. 
                            """)
    parser.add_argument('-bn2', '--basename2', type=str, required=False, default="sample2",
                    help="""
                        Basename of the second sample. Used to create the pileup and extracted features files. Default: 'sample2'
                        """)
    parser.add_argument('-o', '--output', type=str, required=True,
                        help="""
                            Path to output a output directory, in which all output files will be stored.
                            """)
    parser.add_argument('-r', '--reference', type=str, required=True, help='Path to the reference file')

    ### Pileup extractor filters ###
    parser.add_argument('-n', '--num_reads', type=hs.positive_int, required=False,
                        help='Filter by minimum number of reads at a position')
    parser.add_argument('-p', '--perc_mismatched', type=hs.float_between_zero_and_one, required=False, default=0,
                        help='Filter by minimum fraction of mismatched/deleted/inserted bases')
    parser.add_argument('-d', '--perc_deletion', type=hs.float_between_zero_and_one, required=False,
                        help='Filter by minimum percentage of deleted bases')
    parser.add_argument('-q', '--mean_quality', type=hs.positive_float, required=False,
                        help='Filter by mean read quality scores')
    parser.add_argument('-g', '--genomic_region', type=str, required=False,
                        help='Genomic region in "CHR:START-END" format or "CHR" for whole chromosome. Specify to only extract information from a specific region.')
    parser.add_argument('--use_multiprocessing', action="store_true", 
                        help="""
                            Specify whether to use multiprocessing for processing. Recommended for shorter, deeply sequenced data.
                            For low coverage, whole genome data set to false for faster runtime.
                            """)
    parser.add_argument('-t', '--num_processes', type=int, required=False, default=4,
                        help='Number of threads to use for processing.')
    parser.add_argument('--coverage_alt', action="store_true", 
                        help="""
                            Specify which approach should be used to calculate the number of reads a given position.
                            Default: use coverage as calculated by samtools mpileup (i.e. #A+#C+#G+#T+#del).
                            Alternative: calculates coverage considering only matched and mismatched reads, not considering deletions
                            """)
    
    ### Neighbour search arguments ###
    parser.add_argument('-nw', '--window_size', type=int, required=False, default=2,
                        help='Size of the sliding window = 2*w+1. Required when -s flag is set')
    parser.add_argument("-nt", "--neighbour_thresh", type=hs.float_between_zero_and_one, required=False, default=0.5,
                        help="""
                            Threshold of error percentage (--perc_mismatched / --perc_deletion), from which a neighbouring position
                            should be considered as an error.
                            """)
    parser.add_argument('-b', '--n_bins', type=int, required=False, default=5000,
                        help="""Number of bins to split the data into when creating the summary plots. This does not affect the extracted data.
                            Used only to improve performance and clarity of the created plots. Note that setting the value to a low number 
                            can lead to misleading results. Set to '-1' to disable binning. Default: 5000
                            """)
    
    ### Statistical comparison arguments ###
    parser.add_argument('-a', '--alpha', type=float, required=False, default=0.01,
                        help="""
                            Significane level used as a threshold for classifying significant results during the statistical comparison.
                            """)
    parser.add_argument('--no_tsv', action="store_true", 
                        help="""
                            When stated does not store the calculated p-values from the statistical comparison in a TSV file and only as HTML.
                            """)

    return parser

if __name__=="__main__":

    parser = setup_parser()
    args = parser.parse_args()

    processor = Processor(args.sample1, 
                          args.basename,
                          args.output, 
                          args.reference,
                          args.sample2, 
                          args.basename2,
                          num_reads=args.num_reads, 
                          perc_mismatched=args.perc_mismatched,
                          perc_deletion=args.perc_deletion,
                          mean_quality=args.mean_quality,
                          genomic_region=args.genomic_region,
                          use_multiprocessing=args.use_multiprocessing,
                          num_processes=args.num_processes,
                          coverage_alt=args.coverage_alt,
                          window_size=args.window_size,
                          neighbour_thresh=args.neighbour_thresh,
                          n_bins = args.n_bins,
                          alpha = args.alpha,
                          no_tsv = args.no_tsv)
    
    processor.main()