import re, os, warnings, argparse, sys, io
from typing import Dict, List, Tuple, Union, Any
from tqdm import tqdm
from pyfiglet import Figlet
import numpy as np
from multiprocessing import Pool
from itertools import takewhile, repeat
from helper_functions import positive_int, positive_float, float_between_zero_and_one

class FeatureExtractor:
    """
    From a pileup alignment file, extract the following features
    at each genomic position:
    - chromosome
    - position on chromosome
    - number of reads
    - reference base
    - majority base
    - number of A, C, G, T (absolute & relative)
    - number of insertions and deletions (abs. & rel.)
    - motif surrounding a position
    - percentage of reads with an error at a given position
    - mean and standard deviation of the Qscore a a given position

    Also allows for filtering positions by the number of reads, the error percentage,
    the mean quality. Additionally, a specific genomic region can be extracted exlusively.

    The pileup file should be created using Samtool's mpileup method as follows:

    `samtools mpileup -f REF.fa -A -d 0 -Q 0 MAPPING.bam > OUTFILE.pileup` 

    Attributes
    ----------
    input_path : str
        path to a pileup file
    output_path : str
        path to newly created output file
    ref_path : str
        path to the reference fasta file
    ref_sequences : Dict[str, str]
        dictionary containing all sequences (values) and sequence names (keys) 
        extracted from the reference fasta file
    filter_num_reads : int
        filtering option to keep only positions with at least filter_num_reads
        number of reads 
    filter_perc_mismatch : float
        filtering option to keep only position with an error percentage of at 
        least filter_perc_mismatch
    filter_mean_quality : float
        filtering option to keep only position with a mean Qscore of at least
        filter_mean_quality 
    filter_genomic_region : str
        string of format CHR:START-END to specify a genomic region that should
        be processed exclusively
    num_processes
        number of parallel processes during featrue extraction. Decreases 
        computational time

    Methods
    -------
    - get_references(path: str) -> Dict[str, str]
    - extract_positional_info(data_string: str) -> Tuple[str, int, int]
        - region_is_valid(self, chr, start, end) -> bool
    - process_file(self) -> None
        - get_num_lines(self, path: str) -> int
        - process_position(self, line: List[str]) -> str
            - remove_indels(self, pileup_string: str) -> str
            - parse_pileup_string(self, pileup_string: str, ref_base: str) -> Dict[str, Union[str, int]]
            - get_relative_count(self, count_dict: Dict[str, Union[str, int]], n_reads: int) -> Dict[str, Union[str, int, float]]
            - get_majority_base(self, count_dict: Dict[str, Union[str, int, float]]) -> str
            - get_motif(self, chr: str, site: int, ref: str, k: int) -> str
            - get_mismatch_perc(self, count_dict: Dict[str, Union[str, int, float]], ref_base: str) -> int
            - get_read_quality(self, read_qualities: str) -> Tuple[float, float]
    """


    input_path : str
    output_path : str
    ref_path : str
    ref_sequences : Dict[str, str]
    filter_num_reads: int
    filter_perc_mismatch: float
    filter_mean_quality: float
    filter_genomic_region: Tuple[str, int, int]
    num_processes: int
    use_alt_coverage: bool

    def __init__(self, ref_path: str,
                 in_paths: str = None,
                 out_paths: str = None, 
                 num_reads: int = None, 
                 perc_mismatch: float = None,
                 perc_deletion: float = None,
                 mean_quality: float = None,
                 genomic_region: str = None,
                 num_processes: int = None,
                 use_alt_coverage: bool = False) -> None:
        
        self.process_paths(ref_path, in_paths, out_paths)

        self.ref_path = self.check_get_ref_path(ref_path)
        self.input_path = self.check_get_in_path(in_path) if in_path else None
        self.output_path = self.check_get_out_path(out_path, in_path) if out_path else None
    
        self.ref_sequences = self.get_references(ref_path)

        # if no argument is given (e.g. num_reads=None) the minimum number of reads is set to the value that includes all pos.
        self.filter_num_reads = num_reads if num_reads is not None else 1
        self.filter_perc_mismatch = perc_mismatch if perc_mismatch else 0
        self.filter_perc_deletion = perc_deletion if perc_deletion else 0
        self.filter_mean_quality = mean_quality if mean_quality else 0
        self.filter_genomic_region = self.extract_positional_info(genomic_region) if genomic_region else None

        self.num_processes = num_processes

        self.use_alt_coverage = use_alt_coverage

    def __str__(self) -> str:
        return "Pileup extractor"

    #################################################################################################################
    #                                   Functions called during initialization                                      #
    #################################################################################################################

    def process_paths(self, ref_path: str, in_paths: str, out_paths: str) -> None:
        # process path to reference fasta
        self.check_path(ref_path, ["fasta", "fna", "ffn", "faa", "frn", "fa"])

        # process input path(s)
        in_paths = in_paths.split(",")
        for path in in_paths: self.check_path(path, ["msf", "pup", "pileup"]) 
        # process output path(s)  

    def check_path(self, path: str, extensions: List[str]) -> None:
        """
        Check if the specified file path exists and has the expected file extension.

        This function verifies whether the file specified by the given path exists and has a valid extension.
        If the file does not exist, it raises a FileNotFoundError with a detailed error message.
        If the file extension does not match any of the expected extensions, it raises a Warning.

        Parameters:
            path (str): The file path to be checked.
            extensions (List[str]): A list of expected file extensions (e.g., ['txt', 'csv']).

        Raises:
            FileNotFoundError: If the specified file path does not exist.
            Warning: If the file extension is not among the expected extensions.
        """
        if not os.path.exists(path): # does file exist?
            raise FileNotFoundError(f"Input file not found. File '{path}' does not exist.")
        file_type = os.path.splitext(path)[1]
        if not file_type in extensions:
            warnings.warn(f"Found file extension {file_type}. Expected file extension to be one of: {extensions}. If this is deliberate, ignore warning.", Warning)
        


    def check_get_in_path(self, in_path: str, 
                          exp_extensions: List[str] = [".msf", ".pup", ".pileup"],
                          warn_expected_text: str = "Expected .pileup file or similar.") -> str:
        """
        Check if the given input path is valid and of the expected file type.

        Parameters
        ----------
        in_path : str
            Path to the input file given by the user.
        exp_extensions : List[str]
            File extensions of the given input format
        warn_expected_text : str
            Text to be displayed in the warning if another format is given
        Returns
        -------
        str
            Valid input file path.

        Raises
        ------
        FileNotFoundError
            If the input file does not exist.

        Warns
        -----
        UserWarning
            If the input file is not of the expected file type (.msf, .pup, .pileup).
            Warns the user to ensure it is a pileup file as produced by Samtools' mpileup function.
        """
        if not os.path.exists(in_path): # does file exist?
            raise FileNotFoundError(f"Input file not found. File '{in_path}' does not exist.")
        file_type = os.path.splitext(in_path)[1]
        if not file_type in exp_extensions: # is file likely in pileup format?
            warnings.warn(f"Input file of type {file_type}. {warn_expected_text}", Warning)
        
        return in_path

    def check_get_out_path(self, out_path: str, in_path: str, suffix: str = "_extracted.tsv") -> str:
        """
        Check if the given out_put path is valid. Can be either a filename or directory.
        If a directory is given, output path will be '[DIR-path]/[INPUT-FILE-BASENAME]_extracted.tsv'.

        Parameters
        ----------
        out_path : str
            Output path given by the user. Either path to a (non-existing) file or a directory
        in_path : str
            Path to input file given by the user

        Returns
        -------
        str
            Valid path to output file 

        Raises
        ------
        FileNotFoundError
            If the output directory or path to the output file does not exist.
        """
        if os.path.isdir(out_path):
            if not os.path.exists(out_path):
                raise FileNotFoundError(f"Output directory not found. '{out_path}' does not exist.")

            in_basename = os.path.splitext(os.path.basename(in_path))[0]
            if not out_path.endswith("/"):
                out_path += "/"

            return os.path.join(out_path, in_basename + suffix)
        
        else:
            dirname = os.path.dirname(out_path)
            if not os.path.exists(dirname):
                raise FileNotFoundError(f"Path to output file not found. '{dirname}' does not exist.")

            file_extension = os.path.splitext(out_path)[1]
            if file_extension != ".tsv":
                warnings.warn(f"Given output file has extension '{file_extension}'. Note that the output file will be of type '.tsv'.")

            return out_path

    def check_get_ref_path(self, in_path: str) -> str:
        """
        Check if the given path to reference file is valid 
        and of the expected file type.

        Parameters
        ----------
        in_path : str
            Path to the reference file given by the user.

        Returns
        -------
        str
            Valid input file path.

        Raises
        ------
        FileNotFoundError
            If the input file does not exist.

        Warns
        -----
        UserWarning
            If the input file is not of the expected file type (.fa, .fasta).
            Warns the user to ensure it is a fasta file.
        """
        if not os.path.exists(in_path): # does file exist?
            raise FileNotFoundError(f"Input file not found. File '{in_path}' does not exist.")
        file_type = os.path.splitext(in_path)[1]
        if not file_type in [".fa", ".fasta"]: # is file likely in pileup format?
            warnings.warn(f"Reference file of type {file_type}. Make sure that this is a fasta file", Warning)
        return in_path

    def get_references(self, path: str) -> Dict[str, str]:
        """
        Reads a fasta file and stores the sequences in a dictionary (values) with the 
        corresponding chromosome names (keys).

        Parameters
        ----------
        path : str
            filepath to a fasta file

        Returns
        -------
        dict[str]
            Dictionary where the key is the chromosome name and the value is the sequence
        """
        with open(path, "r") as ref:
            lines = ref.readlines()
            i = 0
            refs = {}
            for i in range(len(lines)):
                line = lines[i]
                if line.startswith(">"):
                    refs[line[1:].strip()] = lines[i+1].strip()
        
        return refs
    
    def extract_positional_info(self, data_string: str) -> Tuple[str, int, int]:
        """
        Extracts the chromosome name, start value, and end value from a string in the format "chromosome_name:start-end".

        Parameters
        ----------
        data_string : str
            The input string in the format "chromosome_name:start-end".

        Returns
        -------
        tuple
            A tuple containing the chromosome name (str), start value (int), and end value (int).
        """
        chromosome, positions = data_string.split(':')
        start_str, end_str = positions.split('-')
        start = int(start_str.replace(',', ''))
        end = int(end_str.replace(',', ''))

        if self.region_is_valid(chromosome, start, end):
            return chromosome, start, end
    
    def region_is_valid(self, chr, start, end) -> bool:
        """
        Checks if the chromosome is found in the reference sequences and if so, whether the given start and end
        coordinates are in range of the corresponding sequence.

        Parameters
        ----------
        pos_info : Tuple[str, int, int]
            Positional information extracted in self.extract_positional_info()

        Returns
        -------
        bool
            True, if all information is valid

        Raises
        ------
        Exception, if not all information is valid. 
        """
        # check if chromosome name is found in self.ref_sequences
        if chr not in list(self.ref_sequences.keys()):
            raise Exception(f"Chromosome region error: Name '{chr}' not found in reference sequences from file '{self.ref_path}'")
        # check if start < end
        if start >= end:
            raise Exception(f"Chromosome region error: End position {end} must be larger than start position {start}.")
        # check if start is in range
        chr_len = len(self.ref_sequences[chr])
        if start <= 0 or start > chr_len:
            raise Exception(f"Chromosome region error: Start position {start} not in range of corrdinates 1-{chr_len} (both incl.).")
        # check if end is in range
        if end <= 0 or end > chr_len:
            raise Exception(f"Chromosome region error: End position {end} not in range of corrdinates 1-{chr_len} (both incl.).")
        return True


    #################################################################################################################
    #                                  Functions called during feature extraction                                   #
    #################################################################################################################

    def process_file(self) -> None:
        """
        Reads .pileup file and processes it, writing the results to a new file
        using multiprocessing.

        Parameters
        ----------
        infile : str
            path to the input .pileup file
        outfile : str
            path to the output tsv file
        ref : str
            path to the reference fasta
            
        Returns
        -------
        None
        """
        from_stdin = False if self.input_path else True
        to_stdout = False if self.output_path else True
        def output_line(line: str, output: io.TextIOWrapper = None) -> None:
            if output:
                output.write(line)
            else:
                sys.stdout.write(line)
        def write(file_input):
            with Pool(processes=self.num_processes) as pool:
                o = None if to_stdout else open(self.output_path, "w")

                if not to_stdout:
                    progress_bar = tqdm() if from_stdin else tqdm(total=self.get_num_lines(self.input_path))

                header = f"chr\tsite\tn_reads\tref_base\tmajority_base\tn_a\tn_c\tn_g\tn_t\tn_del\tn_ins\tn_a_rel\tn_c_rel\tn_g_rel\tn_t_rel\tn_del_rel\tn_ins_rel\tperc_mismatch\tmotif\tq_mean\tq_std\n"
                output_line(header, o)
                results = []

                for line in file_input:
                    result = pool.apply_async(self.process_position, args=(line.split("\t"),))
                    results.append((line, result))

                for line, result in results:
                    outline = result.get()
                    if len(outline) > 0:
                        output_line(outline, o)
                    if not to_stdout:
                        progress_bar.update()

                if not to_stdout:
                    o.close()
                    progress_bar.close()

        if not to_stdout:
            f = Figlet(font="slant")
            print(f.renderText("Neet - pileup extractor"))
            print(str(self))

        if from_stdin:
            write(sys.stdin)
        else:
            with open(self.input_path, "r") as i:
                write(i)

    def get_num_lines(self, path: str) -> int:
        """
        Calculate the number of lines in a given file. Function taken from
        https://stackoverflow.com/questions/845058/how-to-get-line-count-of-a-large-file-cheaply-in-python
        
        Parameters
        ----------
        path : str
            Path to a file

        Returns
        -------
        int
            Number of lines in the given file
        """
        f = open(path, 'rb')
        bufgen = takewhile(lambda x: x, (f.raw.read(1024*1024) for _ in repeat(None)))
        return sum( buf.count(b'\n') for buf in bufgen )

    def process_position(self, line: List[str]) -> str:
        """
        Takes a line from a pileup file and processes it.

        Parameters
        ----------
        line : list[str]
            list containing each element from the pileup line.

        Returns
        -------
        str
            New line derived from the initial one. Can be written to a new file in consequent
            steps.
        """
        # extract elements from list
        chr, site, ref_base, read_bases, read_qualities = line[0], int(line[1]), line[2], line[4], line[5]
        
        # filter by genomic region
        region = self.filter_genomic_region
        if region is not None:
            if not(chr == region[0] and site >= region[1] and site <= region[2]): # both start and end inclusive
                return ""

        # extract coverage and filter by number of reads if the standard coverage option is used 
        if not self.use_alt_coverage:
            n_reads = int(line[3])
            if n_reads < self.filter_num_reads: return ""

        # get qualitiy measures
        quality_mean, quality_std = self.get_read_quality(read_qualities)

        # filter by mean read quality
        if quality_mean < self.filter_mean_quality: return ""

        # get reference sequence 
        ref = self.ref_sequences[chr]
        # get absolute number of A, C, G, T, ins, del
        count = self.parse_pileup_string(read_bases, ref_base)

        # in case the alternative way of calculating the coverage is specified
        # could use if else statement and get the other case down here, but then 
        # the count will be calculated each time, potentially wasting time in case the 
        # filter_num_reads is used
        if self.use_alt_coverage:
            n_reads = count["a"]+count["c"]+count["g"]+count["t"]
            if n_reads < self.filter_num_reads: return ""

        # get relative number of A, C, G and T counts
        count = self.get_relative_count(count, n_reads)

        # filter by percentage of deletions
        if count["del_rel"] < self.filter_perc_deletion: return ""

        # get allele fraction
        perc_mismatch = self.get_mismatch_perc(count, ref_base)

        # filter by perc_mismatch
        if perc_mismatch < self.filter_perc_mismatch:
            return ""

        # get majority base
        majority_base = self.get_majority_base(count)

        # get 11b motif
        motif = self.get_motif(chr, site, ref, k=5)

        out = f'{chr}\t{site}\t{n_reads}\t{ref_base}\t{majority_base}\t{count["a"]}\t{count["c"]}\t{count["g"]}\t{count["t"]}\t{count["del"]}\t{count["ins"]}\t{count["a_rel"]}\t{count["c_rel"]}\t{count["g_rel"]}\t{count["t_rel"]}\t{count["del_rel"]}\t{count["ins_rel"]}\t{perc_mismatch}\t{motif}\t{quality_mean}\t{quality_std}\n'
        return out

    def remove_indels(self, pileup_string: str) -> str:
        """
        Takes a pileup string and removes all occurences of the following patterns:
        '\+[0-9]+' for insertions
        '\-[0-9]+' for deletions
        In addition to the pattern itself, remove the following n characters,
        where n is the number specified after + or -.

        Parameters
        ----------
        pileup_string : str
            Pileup string extracted from the fifth column of a pileup file

        Returns
        -------
        str
            Pileup strings with all occurences of the patterns above removed
        """
        pattern = "(\+|\-)[0-9]+"
        
        # get the start and end indices of all found patterns 
        coords = []
        for m in re.finditer(pattern, pileup_string):
            str_len_as_str = pileup_string[m.start()+1:m.end()]
            num_digits = len(str_len_as_str)
            str_len = int(str_len_as_str)
            coords.append((m.start(), m.start()+1+num_digits+str_len))

        # remove the patterns by the indices
        for start, end in reversed(coords): # reverse list as to not shift the index downstream
            pileup_string = pileup_string[:start] + pileup_string[end:]

        return pileup_string

    def parse_pileup_string(self, pileup_string: str, ref_base: str) -> Dict[str, Union[str, int]]:
        """
        Extracts the number of each base called at a given position, as well as the number
        of insertions and deletions. Information is extracted from a pileup string (fifth
        column in a pileup file).

        Parameters
        ----------
        pileup_string : str
            Pileup string extracted from the fifth column of a pileup file
        ref_base : str
            reference base at the position corresponding to the pileup string

        Returns
        -------
        dict
            Dictionary containing the number of A, T, C, G, 
            insertions and deletions.
        """
        pileup_string = pileup_string.lower()
        # remove all occurences of a caret and the following letter (could contain a,c,g,t)
        pileup_string = re.sub(r'\^.', '', pileup_string)

        ref_base = ref_base.lower()
        count_dict = {"a": 0, "t": 0, "c": 0, "g": 0, "del": 0, "ins": 0}
        
        # get number of deletions
        count_dict["del"] = pileup_string.count("*")
        # get number of insertions
        count_dict["ins"] = len(re.findall(r'\+[0-9]+[ACGTNacgtn]+', pileup_string))

        # remove indel patterns to count the number of mismatches correctly
        pileup_string = self.remove_indels(pileup_string)

        # get number of mismatches (i.e. [ACGT])
        count_dict["a"] = pileup_string.count("a")
        count_dict["t"] = pileup_string.count("t")
        count_dict["c"] = pileup_string.count("c")
        count_dict["g"] = pileup_string.count("g")

        # get number of matches (determine where to count matches bases on ref_base)
        n_matches = pileup_string.count('.') + pileup_string.count(',')
        count_dict[ref_base] = n_matches

        # for now removing the count of insertions, as the pileup format cannot provide the information

        return count_dict

    def get_relative_count(self, count_dict: Dict[str, Union[str, int]], n_reads: int) -> Dict[str, Union[str, int, float]]:
        """
        Gets a dictionary containing the absolute counts for A, C, G and T 
        and calculates the relative proportions

        Parameters
        ----------
        count_dict : dict[int]
            Dictionary containing the absolute counts for A, C, G and T
        n_reads : int
            Number of reads at the given position

        Returns
        -------
        dict[float]
            Dictionary containing the relative counts for A, C, G and T
        """
        #n_reads = sum([count_dict["a"], count_dict["c"], count_dict["g"], count_dict["t"]])
        try:
            count_dict["a_rel"] = count_dict["a"] / n_reads
            count_dict["c_rel"] = count_dict["c"] / n_reads
            count_dict["g_rel"] = count_dict["g"] / n_reads
            count_dict["t_rel"] = count_dict["t"] / n_reads
            count_dict["del_rel"] = count_dict["del"] / n_reads
            count_dict["ins_rel"] = count_dict["ins"] / n_reads

        except ZeroDivisionError:
            count_dict["a_rel"] = 0
            count_dict["c_rel"] = 0
            count_dict["g_rel"] = 0
            count_dict["t_rel"] = 0
            count_dict["del_rel"] = 0
            count_dict["ins_rel"] = 0


        return count_dict

    def get_majority_base(self, count_dict: Dict[str, Union[str, int, float]]) -> str:
        """
        Gets a dictionary containing the absolute counts for A, C, G and T and returns the
        key of the one with the highest count.

        Parameters
        ----------
        count_dict : dict
            dictionary containing the absolute counts for A, C, G and T

        Returns
        -------
        str
            Key from the dictionary corresponding to the largest value
        """
        dict_subset = dict((k, count_dict[k]) for k in ("a", "c", "g", "t"))
        return max(dict_subset, key = dict_subset.get).upper()

    def get_motif(self, chr: str, site: int, ref: str, k: int) -> str:
        """
        Extracts the motif of k bases up- and downstream from a given chromosomal site.
        Around the start and end of a refernce sequence the missing bases are filled with
        Ns.

        Parameters
        ----------
        chr : str
            name of the chromosome
        site : int
            position on the chromosome (1-indexed)
        ref : str
            reference sequence for the given chromosome 
        k : int
            number of bases to be regarded in both up- and downstream direction 
            
        Returns
        -------
        str
            sequence of k bases around the center site
        """ 
        idx = site-1
        n_ref = len(ref)

        if idx >= 0 and idx < n_ref:
            idx_l = idx-k
            idx_r = idx+k+1
            # left overhang
            if idx_l < 0:
                len_overhang = abs(idx_l)
                overhang = "N" * len_overhang
                motif = overhang + ref[:idx_r]
            # right overhang
            elif idx_r > n_ref:
                len_overhang = idx_r - n_ref
                overhang = "N" * len_overhang
                motif = ref[idx_l:] + overhang
            # no overhang
            else:
                motif = ref[idx_l:idx_r]

            return motif
        
    def get_mismatch_perc(self, count_dict: Dict[str, Union[str, int, float]], ref_base: str) -> int:
        """
        Calculates the number of reads containing a mismatch, insertion or deletion 
        at a given position.

        Parameters
        ----------
        count_dict : dict
            Dictionary containing the number of occurences of A,C,G,T,ins,del for a given position
        ref_base : str
            reference base at the given position

        Returns
        -------
        int
            Number of mismatched reads a the given position
        """
        mismatch_perc_sum = 0
        for b in ["a", "c", "g", "t"]:
            if b != ref_base.lower():
                mismatch_perc_sum += count_dict[b+"_rel"]

        return mismatch_perc_sum

    def get_read_quality(self, read_qualities: str) -> Tuple[float, float]:
        """
        Calculates the mean and std from the read qualities given in the sixth row
        of a pileup file.

        Parameters
        ----------
        read_qualities : str
            Read quality string from pileup file

        Returns
        -------
        tuple[float, float]
            Mean and standard deviation of read qualities
        """
        # transform string to list of corresponding phred numeric values
        vals = [code - 33 for code in read_qualities.encode("ascii")]

        mean = sum(vals)/len(vals)
        std = np.std(vals)

        return mean, std 
    


if __name__ == "__main__":
    parser = argparse.ArgumentParser(prog="Pileup feature extractor",
                                        description="Extracs different characteristics from a\
                                        given pileup file.")
    parser.add_argument('-r', '--reference', type=str, required=True,
                        help='Path to the reference file')
    parser.add_argument('-i', '--input', type=str, required=True, default="-",
                        help="""
                            Path to the input file(s). If replicates are available, specify paths comma-separated (<repl.1>,<repl.2>,...). 
                            To read from stdin, give "-" for input.
                            """)
    parser.add_argument('-o', '--output', type=str, required=True, default="-",
                        help="""
                            Path to output file(s) (comma-separated if multiple) or directory. If filename(s) are given the order corresponds to the input files. 
                            If a directory is given, the output files are created using the basename from an input file with the suffix "_extracted.tsv".
                            To write to stdout, give "-" for output.
                            """)
    parser.add_argument('-n', '--num_reads', type=positive_int, required=False,
                        help='Filter by minimum number of reads at a position')
    parser.add_argument('-p', '--perc_mismatched', type=float_between_zero_and_one, required=False,
                        help='Filter by minimum fraction of mismatched/deleted/inserted bases')
    parser.add_argument('-d', '--perc_deletion', type=float_between_zero_and_one, required=False,
                        help='Filter by minimum percentage of deleted bases')
    parser.add_argument('-q', '--mean_quality', type=positive_float, required=False,
                        help='Filter by mean read quality scores')
    parser.add_argument('-g', '--genomic_region', type=str, required=False,
                        help='Genomic region in "CHR:START-END" format. Specify to only extract information from a specific region.')
    parser.add_argument('-t', '--num_processes', type=int, required=False,
                        help='Number of threads to use for processing.')
    parser.add_argument('--coverage_alt', action="store_true", 
                        help="""
                            Specify which approach should be used to calculate the number of reads a given position. Default: use coverage as calculated 
                            by samtools mpileup (i.e. #A+#C+#G+#T+#del). Alternative: calculates coverage considering only matched and mismatched reads, 
                            not considering deletions.
                            """)

    args = parser.parse_args()
    
    feature_extractor = FeatureExtractor(args.reference, args.input, args.output, 
                                         num_reads=args.num_reads, 
                                         perc_mismatch=args.perc_mismatched,
                                         perc_deletion=args.perc_deletion,
                                         mean_quality=args.mean_quality,
                                         genomic_region=args.genomic_region,
                                         num_processes=args.num_processes,
                                         use_alt_coverage=args.coverage_alt)
    
    feature_extractor.process_file()