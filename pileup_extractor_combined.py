import re, os, warnings, argparse, sys, io, tempfile
from typing import Dict, List, Tuple, Union, Any
from tqdm import tqdm
from pyfiglet import Figlet
import numpy as np
from multiprocessing import Pool
from helper_functions import positive_int, positive_float, float_between_zero_and_one, get_num_lines

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
    tmp_data: List[str]
    no_neighbour_search: bool
    window_size: int
    neighbour_error_threshold: float

    def __init__(self, ref_path: str,
                 in_path: str = None,
                 out_path: str = None, 
                 num_reads: int = None, 
                 perc_mismatch: float = None,
                 perc_deletion: float = None,
                 mean_quality: float = None,
                 genomic_region: str = None,
                 num_processes: int = None,
                 use_alt_coverage: bool = False,
                 no_neighbour_search: bool = False,
                 window_size: int = None,
                 neighbour_error_threshold: float = None) -> None:
        
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

        self.no_neighbour_search = no_neighbour_search
        self.window_size = window_size
        self.neighbour_error_threshold = neighbour_error_threshold

        self.tmp_data = []

    def __str__(self) -> str:
        o = f"FeatureExtractor instance information:\n\n"
        o += f" - input pileup file: {self.input_path}\n"
        o += f" - writing output file to {self.output_path}\n"
        o += f" - {len(self.ref_sequences)} reference sequence(s) found in file {self.ref_path}\n"
        o += f" - ignoring positions with  <{self.filter_num_reads} reads\n"
        o += f" - ignoring positions with  <{self.filter_perc_mismatch*100}% read errors\n" if self.filter_perc_mismatch>0 else " - no filtering by error percentage\n"
        o += f" - ignoring positions with mean read qualities  <{self.filter_mean_quality}\n" if self.filter_mean_quality>0 else " - no filtering by mean read quality\n"
        o += f" - extracting positions only on chromosome '{self.filter_genomic_region[0]}' from position {self.filter_genomic_region[1]} to {self.filter_genomic_region[2]}\n" if self.filter_genomic_region else " - extracting all genomic positions\n"
        o += f" - extracting number of reads from pileup file" if not self.use_alt_coverage else f" - calculating number of reads from matched & mismatched reads only"
        return o


    #################################################################################################################
    #                                   Functions called during initialization                                      #
    #################################################################################################################

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

        def store_output_line(line: str, output: io.TextIOWrapper = None) -> None:
            """
            Stores the output line based on the neighbour search settings.

            Args:
                line (str): The line to be stored.
                output (io.TextIOWrapper, optional): Output file. If not provided, the line is printed to the standard output. Default is None.

            Returns:
                None

            Raises:
                None

            Notes:
                If 'no_neighbour_search' is True and 'output' is provided, the line is written to the output file.
                If 'no_neighbour_search' is True and 'output' is not provided, the line is printed to the standard output.
                If 'no_neighbour_search' is False, the line is appended to the 'tmp_data' list.
            """
            if (self.no_neighbour_search) & (output is not None):
                output.write(line)
            elif self.no_neighbour_search:
                sys.stdout.write(line)
            else:
                self.tmp_data.append(line)
            
        def write(file_input):
            with Pool(processes=self.num_processes) as pool:
                o = None if (to_stdout | self.no_neighbour_search) else open(self.output_path, "w")

                if not to_stdout:
                    desc = "Processing pileup rows"
                    progress_bar = tqdm(desc=desc) if from_stdin else tqdm(desc=desc, total=get_num_lines(self.input_path))

                header = f"chr\tsite\tn_reads\tref_base\tmajority_base\tn_a\tn_c\tn_g\tn_t\tn_del\tn_ins\tn_a_rel\tn_c_rel\tn_g_rel\tn_t_rel\tn_del_rel\tn_ins_rel\tperc_mismatch\tmotif\tq_mean\tq_std\n"
                store_output_line(header, o)
                results = []

                for line in file_input:
                    result = pool.apply_async(self.process_position, args=(line.split("\t"),))
                    results.append((line, result))

                for line, result in results:
                    outline = result.get()
                    if len(outline) > 0:
                        store_output_line(outline, o)
                    if not to_stdout:
                        progress_bar.update()

                if not to_stdout:
                    o.close()
                    progress_bar.close()

        if not to_stdout:
            f = Figlet(font="slant")
            print(f.renderText("Neet - pileup extractor"))
            # print(str(self))

        if from_stdin:
            write(sys.stdin)
        else:
            with open(self.input_path, "r") as i:
                write(i)

        # start neighbourhood search
        if not self.no_neighbour_search:
            self.read_lines_sliding_window()

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
    
    #################################################################################################################
    #                                  Functions called during neighbour search                                     #
    #################################################################################################################

    def read_lines_sliding_window(self) -> None:
        """
        Reads lines from temporary data using a sliding window approach and performs processing on the lines.
        Writes the processed output to the specified output file or stdout.

        Returns:
            None
        """
        def output_line(line, output: io.TextIOWrapper = None) -> None:
            """
            Writes a line to the specified output file or stdout.

            Args:
                line: The line to be written.
                output: The output file to write the line to. If not provided, writes to stdout.

            Returns:
                None
            """
            if output:
                output.write(line)
            else:
                sys.stdout.write(line)

        def write_edge_lines(neighbourhood: List[str], outfile: io.TextIOWrapper, start: bool = True):
            """
            Writes the neighbourhood information for the first or last rows as returned by the process_edge method.

            Args:
                neighbourhood: A list of lines representing the current neighbourhood.
                outfile: The output file to write the edge lines to.
                start: A boolean indicating if it's the start or end of the neighbourhood.

            Returns:
                None
            """
            k = self.window_size
            r = range(k) if start else range(k+1, 2*k+1)
            for current_pos in r:
                outline = self.process_edge(current_pos, neighbourhood, start)
                output_line(outline, outfile)

        to_stdout = False if self.output_path else True

        window_size_full = 1 + 2 * self.window_size
        lines = []
        first = True

        o = None if to_stdout else open(self.output_path, "w")

        n_lines = len(self.tmp_data)
        header = self.tmp_data.pop(0)
        header = header.strip("\n")+"\thas_neighbour_error\tneighbour_error_pos\n" 
        output_line(header, o)

        if not to_stdout:
            desc = "Searching for nearby errors"
            progress_bar = tqdm(desc=desc, total=n_lines)

        if n_lines < window_size_full:
            lines = []
            for line in self.tmp_data:
                lines.append(line)
            for current_pos in range(n_lines):
                outline = self.process_small(current_pos, lines, n_lines)
                output_line(outline, o)
                if not to_stdout:
                    progress_bar.update()
        else:
            for line in self.tmp_data:
                lines.append(line)

                if len(lines) > window_size_full:  
                    lines.pop(0)
            
                if len(lines) == window_size_full:
                    # once lines hits the max. size for the first time, get and write the first k lines
                    if first:
                        first = False
                        write_edge_lines(lines, o, start=True)

                    outline = self.process_neighbourhood(lines)
                    output_line(outline, o)

                if not to_stdout:
                    progress_bar.update()

            # after the last center line was added to lines, process the last k lines
            write_edge_lines(lines, o, start=False)

        if not to_stdout:
            o.close()
            progress_bar.update()
            progress_bar.close()

    def process_small(self, current_pos: int, neighbourhood: List[str], n_lines: int) -> str:
        """
        Process a small neighbourhood in case the full window size is smaller than the number
        of lines.

        Args:
            current_pos: The current position within the neighbourhood.
            neighbourhood: A list of lines representing the neighbourhood.
            n_lines: The total number of lines in the input.

        Returns:
            The processed line as a string.
        """
        ref_str = neighbourhood[current_pos].strip("\n")
        nb = neighbourhood.copy()
        ref = nb[current_pos].strip("\n").split("\t")
        del(nb[current_pos])

        has_nb, nb_info = self.get_neighbour_info(ref, nb)
        ref_str += f"\t{has_nb}\t{nb_info}\n"
        return ref_str

    def process_edge(self, current_pos: int, neighbourhood: List[str], start: bool = True) -> str:
        """
        Process a neighbourhood at the beginning or the end of the input.

        Args:
            current_pos: The current position within the neighbourhood.
            neighbourhood: A list of lines representing the neighbourhood.
            start: A boolean indicating if it's the start or end of the neighbourhood.

        Returns:
            The processed line as a string.
        """

        k = self.window_size
        ref_str = neighbourhood[current_pos].strip("\n")
        nb = neighbourhood.copy()
        ref = nb[current_pos].strip("\n").split("\t")

        if start:
            nb = nb[:current_pos+k+1]
            del(nb[current_pos])
        else:       
            del(nb[current_pos])
            nb = nb[current_pos-k:]

        has_nb, nb_info = self.get_neighbour_info(ref, nb)
        ref_str += f"\t{has_nb}\t{nb_info}\n"
        return ref_str

    def process_neighbourhood(self, neighbourhood: List[str]) -> str:
        """
        Get 2*window_size+1 rows ([row i-k, ..., row i, ..., row i+k]) that 
        are next to each other in the tsv file and compare the row i to all 
        remaining ones. Check if the other rows are on the same chromosome 
        and if the relative distance btw them is smaller or equal to k.
        Create summary string that indicates the error position relative to 
        the center position.
        Add new information to row and return

        Parameters
        ----------
        neighbourhood : List[str]
            List of k number of rows extracted from a tsv file

        Returns
        -------
        str
            New line containing the neighbourhood information for a given
            center position
        """
        k = self.window_size

        ref_str = neighbourhood[k].strip("\n")
        nb = neighbourhood.copy()

        ref = nb[k].strip("\n").split("\t")
        del nb[k]

        has_nb, nb_info = self.get_neighbour_info(ref, nb)
        ref_str += f"\t{has_nb}\t{nb_info}\n"
        return ref_str

    def get_neighbour_info(self, ref: List[str], neighbourhood: List[str]) -> Tuple[bool, str]:
        """
        From a range of neighbouring lines in a file (neighbourhood), check if positions in these
        lines are neighbours on the reference genome (based on chromosome and site) and if close
        genomic positions can be regarded as errors based on the given error threshold.

        Parameters
        ----------
        ref : List[str]
            line of the central position for which neighbours are searched
        neighbourhood : List[str]
            list containing the lines surrounding the central position

        Returns
        -------
        Tuple[bool, str]
            - boolean indicating if any neighbouring error was found
            - string giving the relative distance to all neighbouring errors to the central position,
              if any were found 
        """
        k = self.window_size
        ref_chr = ref[0]
        ref_site = int(ref[1])

        has_nb = False
        nb_info = ""

        # for each neighbour check if they are 1.on the same chr and 2.neighbours
        for pos in neighbourhood:
            pos = pos.strip("\n").split("\t")
            chr = pos[0]
            perc_error = float(pos[17])

            if (chr == ref_chr) & (perc_error >= self.neighbour_error_threshold): # check if same chromosome & if pos is error
                site = int(pos[1])
                relative_pos = site - ref_site

                if (abs(relative_pos) <= k): # check if pos are close to each other
                    has_nb = True
                    nb_info += str(relative_pos)+","
        return has_nb, nb_info



if __name__ == "__main__":
    parser = argparse.ArgumentParser(prog="Pileup feature extractor",
                                        description="Extracs different characteristics from a\
                                        given pileup file.")
    parser.add_argument('-r', '--reference', type=str, required=True,
                        help='Path to the reference file')
    parser.add_argument('-i', '--input', type=str, required=False,
                        help='Path to the input file. If not specified, read from stdin')
    parser.add_argument('-o', '--output', type=str, required=False,
                        help='Path to the output file. If not specified, write to stdout')
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
                        help='Specify which approach should be used to calculate the number of reads a given position. \
                            Default: use coverage as calculated by samtools mpileup (i.e. #A+#C+#G+#T+#del).\
                            Alternative: calculates coverage considering only matched and mismatched reads, not considering deletions')
    parser.add_argument("--no_neighbour_search", action="store_true",
                        help="Specifies that the neighbourhood search should not be performed. If not specified, \
                            checks for each position, if and where errors appear in the genomic surroundings. \
                            When specified, make sure to add -nw and -nt flag.")
    parser.add_argument('-nw', '--window_size', type=int, required=False, default=2,
                        help='Size of the sliding window = 2*w+1. Required when -s flag is set')
    parser.add_argument("-nt", "--neighbour_thresh", type=float_between_zero_and_one, required=False, default=0.5,
                        help="Threshold of error percentage (--perc_mismatched / --perc_deletion), from which a neighbouring position \
                            should be considered as an error.")

    args = parser.parse_args()
    
    feature_extractor = FeatureExtractor(args.reference, args.input, args.output, 
                                         num_reads=args.num_reads, 
                                         perc_mismatch=args.perc_mismatched,
                                         perc_deletion=args.perc_deletion,
                                         mean_quality=args.mean_quality,
                                         genomic_region=args.genomic_region,
                                         num_processes=args.num_processes,
                                         use_alt_coverage=args.coverage_alt,
                                         no_neighbour_search=args.no_neighbour_search,
                                         window_size=args.window_size,
                                         neighbour_error_threshold=args.neighbour_thresh)
    
    feature_extractor.process_file()