import re, os, warnings, argparse, sys, io, random
from typing import Dict, List, Tuple, Union, Any
from tqdm import tqdm
import numpy as np
from multiprocessing import Pool
from . import helper_functions as hs
from .summary import SummaryCreator
class FeatureExtractor:

    input_paths : List[str]
    output_paths : List[str]
    ref_path : str
    ref_sequences : Dict[str, str]
    
    filter_num_reads: int
    filter_perc_mismatch: float
    filter_perc_mismatch_alt: float
    filter_mean_quality: float
    filter_genomic_region: Tuple[str, int, int] | None
    num_processes: int
    window_size: int
    neighbour_error_threshold: float
    n_bins_summary: int|None
    use_alt_summary: bool

    def __init__(self, 
                 in_paths: str,
                 out_paths: str, 
                 ref_path: str,
                 num_reads: int, 
                 mismatch_rate: float | None,
                 mismatch_rate_alt: float | None,
                 perc_deletion: float | None,
                 mean_quality: float | None,
                 genomic_region: str | None,
                 use_multiprocessing: bool,
                 num_processes: int,
                 window_size: int,
                 neighbour_error_threshold: float,
                 n_bins_summary: int,
                 use_alt_summary: bool,
                 export_svg_summary: bool) -> None:

        self.process_paths(ref_path, in_paths, out_paths)    

        # if one of the following arguments was not provided (i.e. arg=None), set variable to a value so nothing gets filtered out
        self.filter_num_reads = num_reads if num_reads is not None else 1
        self.filter_perc_mismatch = mismatch_rate if mismatch_rate else 0
        self.filter_perc_mismatch_alt = mismatch_rate_alt if mismatch_rate_alt else 0
        self.filter_perc_deletion = perc_deletion if perc_deletion else 0
        self.filter_mean_quality = mean_quality if mean_quality else 0
        self.filter_genomic_region = self.extract_positional_info(genomic_region) if genomic_region else None

        self.use_multiprocessing = use_multiprocessing
        self.num_processes = num_processes

        self.window_size = window_size
        self.neighbour_error_threshold = neighbour_error_threshold

        self.n_bins_summary = n_bins_summary if n_bins_summary > 0 else None 
        self.use_alt_summary = use_alt_summary
        self.export_svg_summary = export_svg_summary

    def __str__(self) -> str:
        return ""


    #################################################################################################################
    #                                   Functions called during initialization                                      #
    #################################################################################################################

    def process_paths(self, ref_path: str, in_paths_str: str, out_paths: str) -> None:
        """
        Process input and output file paths for the data processing.

        This method processes the file paths provided for reference sequence, input data,
        and output data. It checks the validity of the paths and stores the processed
        information in the respective instance variables.

        Parameters:
            ref_path (str): The file path to the reference fasta file.
            in_paths (str): Comma-separated list of input file paths.
            out_paths (str): Output file path or directory path.

        Returns:
            None

        Raises:
            FileNotFoundError: If the reference file or any of the input files do not exist.
            Warning: If the file extension of any input file is not among the expected extensions.

        Note:
            - The reference fasta file should have one or more sequence entries in the format:
            >chromosome_name description(optional)
            sequence
            - The input files are expected to have specific extensions (e.g., '.msf', '.pup', '.pileup').
            - The output files will be generated with the '.tsv' extension.
        """
        # process input path(s)
        in_paths = in_paths_str.split(",")
        for path in in_paths: hs.check_input_path(path, [".msf", ".pup", ".pileup"])
        self.input_paths = in_paths

        # process output path(s)
        self.output_paths = self.process_outpaths(out_paths)

        # process path to reference fasta
        hs.check_input_path(ref_path, [".fasta", ".fna", ".ffn", ".faa", ".frn", ".fa"])
        self.ref_sequences = hs.get_references(ref_path)

    def process_outpaths(self, out: str) -> List[str]:
        """
        Process the output file paths based on the input `out` argument.

        Args:
            out (str): The output file path or directory path. If set to "-", no output file paths
                    will be generated, and the function will return `None`. If `out` is a
                    directory path, output file paths will be generated using the input file
                    names with the `.tsv` extension appended. If `out` is a comma-separated list
                    of filenames, the function will verify if the specified files exist and if
                    their extensions are `.tsv`.

        Returns:
            List[str] | None: A list of output file paths if `out` is a directory path or a
                            comma-separated list of filenames. If `out` is set to "-", the
                            function returns `None`.

        Raises:
            FileNotFoundError: If `out` is a directory path, and the directory does not exist.
                            If `out` is a list of filenames, and the directory of any file
                            does not exist.
            UserWarning: If any filename in the list has an extension other than `.tsv`, a
                        user warning will be issued to inform the user that the output file
                        will be of type `.tsv`.
        """
        out_paths = out.split(",")
        if len(out_paths) > 1:
            for path in out_paths:
                hs.check_create_dir(os.path.dirname(path))
                file_extension = os.path.splitext(path)[1]
                if file_extension != ".tsv":
                    warnings.warn(f"Given output file has extension '{file_extension}'. Note that the output file will be of type '.tsv'.")
        else: 
            hs.check_create_dir(out)
            out_paths = []
            for in_path in self.input_paths:
                basename = os.path.splitext(os.path.basename(in_path))[0]
                out_paths.append(os.path.join(out, f"{basename}_extracted.tsv"))
        return out_paths

    def extract_positional_info(self, data_string: str) -> Tuple[str, int, int] | None:
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
        if (":" in data_string) & ("-" in data_string): 
            chromosome, positions = data_string.split(':')
            start_str, end_str = positions.split('-')
            start = int(start_str.replace(',', ''))
            end = int(end_str.replace(',', ''))
        else:
            chromosome = data_string
            start = 1
            end = len(self.ref_sequences[chromosome])

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
    def main(self) -> None:
        """
        Process multiple pairs of input and output files. Reads .pileup files and processes them in parallel using multiprocessing.

        Returns:
            None
        """ 
        for in_file, out_file in zip(self.input_paths, self.output_paths):
            hs.print_update(f"Processing file '{in_file}'. Writing to '{out_file}'.")
            self.process_file(in_file, out_file)
            self.create_summary_file(out_file)

    def process_file(self, in_file: str, out_file: str) -> None:
        """
        Reads a .pileup file, processes it in a sliding window to incorporate the neighbourhood search,
        and writes the results to a new file.

        Parameters
        ----------
        in_file : str
            Path to the input .pileup file.
        out_file : str
            Path to the output tsv file.

        Returns
        -------
        None
        """  

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
                outfile.write(outline)

        def write_center_line(neighbourhood: List[str], outfile: io.TextIOWrapper):
            """
            Writes the neighbourhood information for the center position of a full-sized neighbourhood.

            Args:
                neighbourhood: A list of lines representing the current neighbourhood.
                outfile: The output file to write the edge lines to.

            Returns:
                None
            """
            outline = self.process_neighbourhood(neighbourhood)
            outfile.write(outline)

        with open(in_file, "r") as i, open(out_file, "w") as o:
            desc = "Processing pileup rows"
            hs.print_update("Counting number of lines to process.")
            progress_bar = tqdm(desc=desc, total=hs.get_num_lines(in_file))

            header = f"chr\tsite\tn_reads\tref_base\tmajority_base\tn_a\tn_c\tn_g\tn_u\tn_del\tn_ins\tn_ref_skip\ta_rate\tc_rate\tg_rate\tu_rate\tdeletion_rate\tinsertion_rate\trefskip_rate\tmismatch_rate\tmismatch_rate_alt\tmotif\tq_mean\tq_std\tneighbour_error_pos\n"
            o.write(header)
            
            nb_size_full = 1 + 2 * self.window_size
            nb_lines = []
            nb_first = True

            if self.use_multiprocessing:
                with Pool(processes=self.num_processes) as pool:
                    for outline in pool.imap(self.process_position, i):
                        if len(outline) > 0: 
                            nb_lines.append(outline)

                            if len(nb_lines) > nb_size_full:  
                                nb_lines.pop(0)

                            if len(nb_lines) == nb_size_full:
                                if nb_first:
                                    nb_first = False
                                    write_edge_lines(nb_lines, o, start=True)
                                write_center_line(nb_lines, o)
                        progress_bar.update()
            else:
                for line in i:
                    outline = self.process_position(line) # extracting the features themselves
                    if len(outline) > 0: 
                        nb_lines.append(outline)

                        if len(nb_lines) > nb_size_full:  
                            nb_lines.pop(0)

                        if len(nb_lines) == nb_size_full:
                            if nb_first:
                                nb_first = False
                                write_edge_lines(nb_lines, o, start=True)
                            write_center_line(nb_lines, o)
                    progress_bar.update()

            if len(nb_lines) < nb_size_full:
                for current_pos in range(len(nb_lines)):
                    outline = self.process_small(current_pos, nb_lines)
                    o.write(outline)
            else:
                write_edge_lines(nb_lines, o, start=False)

            progress_bar.update()
            progress_bar.close()

    def create_summary_file(self, file_path: str) -> None:
        """
        Create a summary file from the newly created tsv file.

        Parameters:
            file_path (str): Path to the newly created tsv output.

        Returns:
            None
        """
        out_path = os.path.splitext(file_path)[0]+"_summary.html"
        summary_creator = SummaryCreator(file_path, out_path, n_bins=self.n_bins_summary, 
                                         use_perc_mismatch_alt=self.use_alt_summary,
                                         export_svg=self.export_svg_summary)
        summary_creator.main()

    def process_position(self, line_str: str) -> str:
        """
        Processes a single position from the pileup data.

        Parameters:
            pileup_data (str): A string containing tab-separated data for a single position from the pileup file.

        Returns:
            str: The processed line as a string.
        """
        line = line_str.split("\t")
        # extract elements from list
        try:
            chr, site, ref_base, read_bases, read_qualities = line[0], int(line[1]), line[2].replace("T", "U"), line[4], line[5]
        except:
            return ""
            
        # filter by genomic region
        region = self.filter_genomic_region
        if region is not None:
            if not(chr == region[0] and site >= region[1] and site <= region[2]): # both start and end inclusive
                return ""

        # extract coverage and filter by number of reads if the standard coverage option is used 
        n_reads = int(line[3])
        if n_reads < self.filter_num_reads: return ""

        # get reference sequence 
        try:
            ref = self.ref_sequences[chr]
        except:
            return ""
        # get absolute number of A, C, G, U, ins, del
        count, ref_skip_positions = self.parse_pileup_string(read_bases, ref_base)

        # get qualitiy measures
        quality_mean, quality_std = self.get_read_quality(read_qualities, ref_skip_positions)
        # filter by mean read quality
        if quality_mean < self.filter_mean_quality: return ""

        # in case the alternative way of calculating the coverage is specified
        # could use if else statement and get the other case down here, but then 
        # the count will be calculated each time, potentially wasting time in case the 
        # filter_num_reads is used
        n_reads_alt = count["a"]+count["c"]+count["g"]+count["u"]

        # get relative number of A, C, G and U counts
        count_rel = self.get_relative_count(count, n_reads)
        count_rel_alt = self.get_relative_count(count, n_reads_alt)

        # filter by percentage of deletions
        if count_rel["del"] < self.filter_perc_deletion: return ""

        # get allele fraction
        mismatch_rate = self.get_mismatch_perc(count_rel, ref_base)
        mismatch_rate_alt = self.get_mismatch_perc(count_rel_alt, ref_base)

        # filter by mismatch_rate
        if mismatch_rate < self.filter_perc_mismatch:
            return ""
        if mismatch_rate_alt < self.filter_perc_mismatch_alt:
            return ""

        # get majority base
        majority_base = self.get_majority_base(count)

        # get 11b motif
        motif = self.get_motif(chr, site, ref, k=2)

        out = f'{chr}\t{site}\t{n_reads}\t{ref_base}\t{majority_base}\t{count["a"]}\t{count["c"]}\t{count["g"]}\t{count["u"]}\t{count["del"]}\t{count["ins"]}\t{count["ref_skip"]}\t{count_rel["a"]}\t{count_rel["c"]}\t{count_rel["g"]}\t{count_rel["u"]}\t{count_rel["del"]}\t{count_rel["ins"]}\t{count_rel["ref_skip"]}\t{mismatch_rate}\t{mismatch_rate_alt}\t{motif}\t{quality_mean}\t{quality_std}\n'
        return out

    def remove_indels(self, pileup_string: str) -> str:
        """
        Takes a pileup string and removes all occurences of the following patterns:
        '\\+[0-9]+' for insertions
        '\\-[0-9]+' for deletions
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
        pattern = "(\\+|\\-)[0-9]+"
        
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

    def parse_pileup_string(self, pileup_string: str, ref_base: str) -> Tuple[Dict[str, int], List[int]]:
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
            Dictionary containing the number of A, U, C, G, 
            insertions and deletions.
        """
        pileup_string = pileup_string.lower()
        # remove all occurences of a caret and the following letter (could contain a,c,g,t)
        pileup_string = re.sub(r'\^.', '', pileup_string)

        ref_base = ref_base.lower()
        count_dict = {"a": 0, "u": 0, "c": 0, "g": 0, "del": 0, "ins": 0}
        
        # get number of deletions
        count_dict["del"] = pileup_string.count("*")
        # get number of insertions
        count_dict["ins"] = len(re.findall(r'\+[0-9]+[ACGTNacgtn]+', pileup_string))

        # remove indel patterns to count the number of mismatches correctly
        pileup_string = self.remove_indels(pileup_string)

        # get number of mismatches (i.e. [ACGT])
        count_dict["a"] = pileup_string.count("a")
        count_dict["u"] = pileup_string.count("t")
        count_dict["c"] = pileup_string.count("c")
        count_dict["g"] = pileup_string.count("g")

        # get number of matches (determine where to count matches bases on ref_base)
        n_matches = pileup_string.count('.') + pileup_string.count(',')
        count_dict[ref_base] = n_matches

        # get number of reference skips
        n_ref_skips = pileup_string.count("<") + pileup_string.count(">")
        count_dict["ref_skip"] = n_ref_skips

        # get the indices from the > & < positions, to filter out of the quality string later on
        # (the corresponding positions refer to the read quality, not the quality of the position on the reads)
        pileup_string = re.sub(r'\*', "", pileup_string)

        ref_skip_idc = [i for i, char in enumerate(pileup_string) if (char==">") | (char=="<")]

        return count_dict, ref_skip_idc

    def get_relative_count(self, count_dict: Dict[str, int], n_reads: int) -> Dict[str, float]:
        """
        Gets a dictionary containing the absolute counts for A, C, G and U
        and calculates the relative proportions

        Parameters
        ----------
        count_dict : dict[int]
            Dictionary containing the absolute counts for A, C, G and U
        n_reads : int
            Number of reads at the given position

        Returns
        -------
        dict[float]
            Dictionary containing the relative counts for A, C, G and U
        """
        #n_reads = sum([count_dict["a"], count_dict["c"], count_dict["g"], count_dict["u"]])
        rel_dict = {}
        for category in ["a", "c", "g", "u", "del", "ins", "ref_skip"]:
            try:
                rel_dict[category] = count_dict[category] / n_reads
            except:
                rel_dict[category] = 0

        return rel_dict

    def get_majority_base(self, count_dict: Dict[str, int]) -> str:
        """
        Gets a dictionary containing the absolute counts for A, C, G and U and returns the
        key of the one with the highest count.

        Parameters
        ----------
        count_dict : dict
            dictionary containing the absolute counts for A, C, G and U

        Returns
        -------
        str
            Key from the dictionary corresponding to the largest value
        """
        dict_subset = dict((k, count_dict[k]) for k in ("a", "c", "g", "u"))
        dict_subset["-"] = count_dict["del"]

        return max(dict_subset, key = lambda k: dict_subset[k]).upper()

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

            return motif.replace("T", "U")
        return ""
        
    def get_mismatch_perc(self, count_dict_rel: Dict[str, float], ref_base: str) -> float:
        """
        Calculates the number of reads containing a mismatch, insertion or deletion 
        at a given position.

        Parameters
        ----------
        count_dict : dict
            Dictionary containing the number of occurences of A,C,G,U,ins,del for a given position
        ref_base : str
            reference base at the given position

        Returns
        -------
        int
            Number of mismatched reads a the given position
        """
        mismatch_perc_sum = 0
        for b in ["a", "c", "g", "u"]:
            if b != ref_base.lower():
                mismatch_perc_sum += count_dict_rel[b]

        return mismatch_perc_sum

    def get_read_quality(self, read_qualities: str, ref_skip_positions: List[int]) -> Tuple[float, float]:
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
        # remove quality values corresponding to reference skips
        read_qualities_len = len(read_qualities)
        read_qualities_list = list(read_qualities)
        ref_skip_positions.reverse()
        
        for i in ref_skip_positions:
            if i < read_qualities_len: # when >/< comes at the end, the quality values don't seem to be added to the string
                read_qualities_list.pop(i)
        read_qualities = "".join(read_qualities_list)

        if len(read_qualities) > 0:
            # transform string to list of corresponding phred numeric values
            vals = [code - 33 for code in read_qualities.encode("ascii")]
            
            mean = np.mean(vals).astype(float)
            std = np.std(vals).astype(float)
        else:
            mean = -1
            std = -1

        return mean, std 
    
    #################################################################################################################
    #                                  Functions called during neighbour search                                     #
    #################################################################################################################
    def process_small(self, current_pos: int, neighbourhood: List[str]) -> str:
        """
        Process a small neighbourhood in case the full window size is smaller than the number of lines.

        Parameters:
            current_pos (int): The current position within the neighbourhood.
            neighbourhood (List[str]): A list of lines representing the neighbourhood.
            n_lines (int): The total number of lines in the input.

        Returns:
            str: The processed line as a string.
        """        
        ref_str = neighbourhood[current_pos].strip("\n")
        nb = neighbourhood.copy()
        ref = nb[current_pos].strip("\n").split("\t")
        del(nb[current_pos])

        nb_info = self.get_neighbour_info(ref, nb)
        ref_str += f"\t{nb_info}\n"
        return ref_str

    def process_edge(self, current_pos: int, neighbourhood: List[str], start: bool = True) -> str:
        """
        Process a neighbourhood at the beginning or the end of the input.

        Parameters:
            current_pos (int): The current position within the neighbourhood.
            neighbourhood (List[str]): A list of lines representing the neighbourhood.
            start (bool, optional): A boolean indicating if it's the start or end of the neighbourhood. Defaults to True.

        Returns:
            str: The processed line as a string.
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

        nb_info = self.get_neighbour_info(ref, nb)
        ref_str += f"\t{nb_info}\n"
        return ref_str

    def process_neighbourhood(self, neighbourhood: List[str]) -> str:
        """
        Get 2*window_size+1 rows ([row i-k, ..., row i, ..., row i+k]) that are next to each other in the tsv file and compare the row i to all 
        remaining ones. Check if the other rows are on the same chromosome and if the relative distance between them is smaller or equal to k.
        Create a summary string that indicates the error position relative to the center position.
        Add new information to the row and return.

        Parameters:
            neighbourhood (List[str]): List of k number of rows extracted from a tsv file.

        Returns:
            str: New line containing the neighbourhood information for a given center position.
        """        
        k = self.window_size

        ref_str = neighbourhood[k].strip("\n")
        nb = neighbourhood.copy()

        ref = nb[k].strip("\n").split("\t")
        del nb[k]

        nb_info = self.get_neighbour_info(ref, nb)
        ref_str += f"\t{nb_info}\n"
        return ref_str

    def get_neighbour_info(self, ref: List[str], neighbourhood: List[str]) -> str:
        """
        From a range of neighbouring lines in a file (neighbourhood), check if positions in these lines are neighbours on the reference genome 
        (based on chromosome and site) and if close genomic positions can be regarded as errors based on the given error threshold.

        Parameters:
            ref (List[str]): Line of the central position for which neighbours are searched.
            neighbourhood (List[str]): List containing the lines surrounding the central position.

        Returns:
            str: A string giving the relative distance to all neighbouring errors to the central position, if any were found.
        """        
        k = self.window_size
        ref_chr = ref[0]
        ref_site = int(ref[1])

        nb_info = ""

        # for each neighbour check if they are 1.on the same chr and 2.neighbours
        for pos in neighbourhood:
            pos = pos.strip("\n").split("\t")
            chr = pos[0]
            perc_error = float(pos[19])

            if (chr == ref_chr) & (perc_error >= self.neighbour_error_threshold): # check if same chromosome & if pos is error
                site = int(pos[1])
                relative_pos = site - ref_site

                if (abs(relative_pos) <= k): # check if pos are close to each other
                    nb_info += str(relative_pos)+","
        return nb_info