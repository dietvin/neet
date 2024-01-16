from typing import List, Dict, Tuple
import os, sys, warnings
import plotly.graph_objects as go
import plotly.express as px
from plotly.subplots import make_subplots
from plotly.io import to_html
from intervaltree import IntervalTree
import datetime, argparse
from collections import defaultdict
from statistics import median
from scipy.stats import normaltest, bartlett, ttest_ind, mannwhitneyu
from . import helper_functions as hs

class POIAnalyzer():

    in_path: str
    bed_path: str
    ref_path: str

    perc_mismatch_col: str
    output_tsv: bool
    export_svg: bool

    data: Dict[str, List[str]] | Dict[str, List[int]] | Dict[str, List[float]] | Dict[str, List[str | None]]
    bed_categories: List[str]
    bed_categories_counterparts: List[str] | List[None]
    output_paths: List[str]

    category_data: Dict[str, List[str]] | Dict[str, List[int]] | Dict[str, List[float]] | Dict[str, List[str | None]]
    counterpart_data: Dict[str, List[str]] | Dict[str, List[int]] | Dict[str, List[float]] | Dict[str, List[str | None]] | None
    current_category: str
    current_counterpart: str | None 
    
    def __init__(self, in_path: str, out_path: str, bed_path: str, ref_path: str,
                 categories: str,
                 canonical_counterpart: str | None,
                 output_tsv: bool = True, 
                 use_perc_mismatch_alt: bool = False,
                 export_svg: bool = False) -> None:
        """
        Initialize an instance of the GenomicDataProcessor class.

        Args:
            in_path (str): The path to the input TSV file containing genomic data.
            out_path (str): The path to the output directory where results will be saved.
            bed_path (str): The path to the BED file containing genomic intervals.
            ref_path (str): The path to the reference FASTA file.
            categories (str): A comma-separated string containing bed categories.
            canonical_counterpart (str): A comma-separated string containing bases corresponding to bed categories.
            output_tsv (bool, optional): Whether to output results as a TSV file. Default is True.
            use_perc_mismatch_alt (bool, optional): Whether to use 'mismatch_rate_alt' instead of 'mismatch_rate' column. Default is False.

        Returns:
            None

        Note:
            This constructor initializes an instance of the GenomicDataProcessor class. It sets various attributes
            based on the provided arguments and performs necessary data processing and validation steps.

        """
        self.process_path(in_path, bed_path, ref_path)
        self.load_data()
        self.add_bed_info()
        self.perc_mismatch_col = "mismatch_rate_alt" if use_perc_mismatch_alt else "mismatch_rate"

        self.process_categories(categories, canonical_counterpart, out_path)

        self.output_tsv = output_tsv
        self.export_svg = export_svg

    ##############################################################################################################
    #                                           Initialization methods                                           #
    ##############################################################################################################

    def process_path(self, in_path: str, bed_path: str, ref_path: str) -> None:
        """
        Process and validate input and output paths for data processing.

        Args:
            in_path (str): The path to the input TSV file containing genomic data.
            out_path (str): The path to the output directory where results will be saved.
            bed_path (str): The path to the BED file containing genomic intervals.
            ref_path (str): The path to the reference FASTA file.

        Returns:
            None

        Note:
            This method validates and sets the paths required for data processing. It checks the existence and
            file extensions of the provided paths, raising exceptions or issuing warnings as appropriate. After
            validation, it sets the 'in_path', 'bed_path', 'ref_path', and 'output_path' attributes of the class
            for use in subsequent processing steps.

        """
        hs.check_input_path(in_path, [".tsv"])
        self.in_path = in_path
        hs.check_input_path(bed_path, [".bed"])
        self.bed_path = bed_path
        hs.check_input_path(ref_path, [".fasta", ".fa", ".fn"])
        self.ref_path = ref_path

    def load_data(self) -> None:
        """
        Load data from the specified input file into the SummaryCreator instance.

        Reads the data from a tab-separated values (tsv) file as created by the PileupExtractor module
        and stores it in the 'data' attribute of the class.

        Returns:
        - None
        """
        cols = ["chr", "site", "n_reads", "ref_base", "majority_base", "a_rate", "c_rate", "g_rate", "u_rate", "deletion_rate", "insertion_rate", "refskip_rate", "mismatch_rate", "q_mean", "motif", "neighbour_error_pos"]

        col_idx = {'chr': 0, 'site': 1, 'n_reads': 2, 'ref_base': 3, 'majority_base': 4, 'n_a': 5, 'n_c': 6, 'n_g': 7, 'n_t': 8, 'n_del': 9, 'n_ins': 10, 'n_ref_skip': 11, 'a_rate': 12, 'c_rate': 13, 'g_rate': 14, 'u_rate': 15, 'deletion_rate': 16, 'insertion_rate': 17, 'refskip_rate': 18, 'mismatch_rate': 19, 'mismatch_rate_alt': 20, 'motif': 21, 'q_mean': 22, 'q_std': 23, 'neighbour_error_pos': 24}
        dtypes = {'chr': str, 'site': int, 'n_reads': int, 'ref_base': str, 'majority_base': str, 'n_a': int, 'n_c': int, 'n_g': int, 'n_t': int, 'n_del': int, 'n_ins': int, 'n_ref_skip': int, 'a_rate': float, 'c_rate': float, 'g_rate': float, 'u_rate': float, 'deletion_rate': float, 'insertion_rate': float, 'refskip_rate': float, 'mismatch_rate': float, 'mismatch_rate_alt': float, 'motif': str, 'q_mean': float, 'q_std': float, 'neighbour_error_pos': str}
        
        with open(self.in_path, "r") as file:
            next(file)
            data = dict(zip(cols, [[] for _ in cols]))
            for line in file:
                line = line.strip().split("\t")
                for col in cols[:-1]:
                    data[col].append(dtypes[col](line[col_idx[col]]))
                if len(line) < 25:
                    data["neighbour_error_pos"].append(None)
                else:
                    data["neighbour_error_pos"].append(line[24])

            self.data = data
    
    def add_bed_info(self) -> None:
        """
        Add BED information to a DataFrame based on a BED file.

        Returns:
            pd.DataFrame: The DataFrame with an additional 'bed_name' column containing BED information.

        Note:
            This method takes a DataFrame 'data' containing genomic data and adds an additional column 'bed_name'
            to it. The 'bed_name' column is populated by querying a previously built IntervalTree from the
            specified 'bed_path' for each row in the DataFrame. The method returns the enriched DataFrame.

        """
        tree = self.build_interval_tree(self.bed_path)
        self.data["bed_name"] = [self.get_name_for_position(chrom, site, tree) for chrom, site in zip(self.data["chr"], self.data["site"])]

    def get_name_for_position(self, chrom: str, site: int, interval_tree: IntervalTree) -> str|None:
        """
        Retrieve the name associated with a given genomic position from an interval tree.

        Args:
            position (tuple): A tuple containing the chromosome and position to query (e.g., ('chr1', 100)).
            interval_tree (IntervalTree): An IntervalTree data structure containing genomic intervals.

        Returns:
            str | None: The name associated with the provided position if found, or None if not found.

        Note:
            This method queries the provided 'interval_tree' for intervals that overlap with the specified
            genomic position. It iterates through the results and checks if the chromosome matches the
            provided chromosome in the 'position' tuple. If a matching interval is found, the associated
            name is returned. If no matching interval is found, None is returned.

        """
        results = interval_tree[site:site+1]  # Query the interval tree
        for interval in results:
            chromosome, name = interval.data
            if chromosome == chrom:
                return name
        return None

    def build_interval_tree(self, bed_file: str) -> IntervalTree:
        """
        Build an IntervalTree data structure from a BED file.

        Args:
            bed_file (str): The path to the BED file containing genomic intervals.

        Returns:
            IntervalTree: An IntervalTree data structure containing the genomic intervals from the BED file.

        Note:
            This method reads the specified BED file and extracts genomic intervals along with optional
            associated names. It constructs an IntervalTree data structure to efficiently query and store
            these intervals. The method accounts for the 0-based indexing convention in BED files by adding
            1 to the start and end positions. If a name is provided in the BED file, it is associated with
            the corresponding interval. The constructed IntervalTree is returned.

        """
        interval_tree = IntervalTree()
        with open(bed_file, 'r') as file:
            for line in file:
                parts = line.strip().split('\t')
                if len(parts) >= 3:
                    chromosome = parts[0]
                    start = int(parts[1])+1 # +1 to account for 0-index in BED files 
                    end = int(parts[2])+1 # +1 to account for 0-index in BED files
                    name = parts[3] if len(parts) >= 4 else None
                    interval_tree[start:end] = (chromosome, name)
        return interval_tree

    def process_categories(self, categories: str, counterparts: str|None, out_paths: str) -> None:
        """
        Process categories, counterparts, and output paths.

        Parameters:
        - categories (str): Comma-separated string of categories to be processed.
        - counterparts (str|None): Comma-separated string of counterparts or None if not provided.
        - out_paths (str): Comma-separated string of output paths or a single directory path.

        Raises:
        - Exception: If a specified category is not found in the bed file.
        - Exception: If the number of counterparts does not match the number of categories.
        - Exception: If the number of output paths does not match the number of categories.

        Note:
        - If only one output path is provided, it is considered a directory, and individual
        output paths for each category are generated based on this directory.
        """
        # process given categories
        cat_split = categories.split(",")
        unique_cat = list(sorted([i for i in set(self.data["bed_name"]) if i]))
        for category in cat_split:
            if category not in unique_cat:
                raise Exception(f"Given name '{category}' was not found in the bed file.")
            
        # if given, process counterparts
        if counterparts:
            cou_split = counterparts.upper().split(",")
            if len(cou_split) != len(cat_split):
                raise Exception(f"For the {len(cat_split)} categories {len(cou_split)} corresponding bases were given. Each category must have a base it corresponds to.")
        else:
            cou_split = [None for _ in self.bed_categories]

        # process ouput paths
        out_split = out_paths.split(",")
        if len(out_split) > 1:
            if len(out_split) != len(cat_split): # output for each category?
                raise Exception(f"For the {len(cat_split)} categories {len(out_split)} output paths were given. Each category must have an output path it corresponds to.")
            else:
                for out in out_split:
                    hs.check_create_dir(os.path.dirname(out))
        else: # output string must be a directory
            dirname = out_paths
            hs.check_create_dir(dirname)
            in_basename = os.path.splitext(os.path.basename(self.in_path))[0]
            out_split = [os.path.join(dirname, f"{c}_{in_basename}_summary.html") for c in cat_split]
    
        self.bed_categories = cat_split
        self.bed_categories_counterparts = cou_split
        self.output_paths = out_split

    ##############################################################################################################
    #                                           Main processing methods                                          #
    ##############################################################################################################

    def main(self):
        """
        Main entry point of the script.

        This method iterates through the bed categories and their corresponding bases, processing each category using
        the `process_category` method.

        Returns:
            None: The script performs the desired processing and saves the output files.
        """
        for category, corr_base, output_path in zip(self.bed_categories, self.bed_categories_counterparts, self.output_paths):
            self.process_category(category, corr_base, output_path)

    def process_category(self, category: str, corresponding_base: str|None, output_path: str) -> None:
        """
        Process and analyze the data for a specific category and (optionally) corresponding base.

        This method generates multiple plots and saves them along with a summary HTML template. Optionally, it also
        writes the processed data to a TSV file.

        Parameters:
            category (str): The category of interest.
            corresponding_base (str): The corresponding base for the category.

        Returns:
            None: The plots and summary are saved to files, and optionally, the data is saved in TSV format.
        """
        hs.print_update(f"Processing {category} sites.")

        hs.print_update("   - subsetting data... ", line_break=False)
        self.subset_category_data(category, corresponding_base)
        if corresponding_base:
            hs.print_update(f"Done. Found {len(self.category_data['chr'])} {category} and {len(self.counterpart_data['chr'])} {corresponding_base} sites.", with_time=False)
        else:
            hs.print_update(f"Done. Found {len(self.category_data['chr'])} {category} sites.", with_time=False)

        hs.print_update("   - creating position overview... ", line_break=False)
        plot_mod_map = self.create_map_plot()
        hs.print_update("Done", with_time=False)

        hs.print_update(f"   - creating overview of mismatch types at {category} sites... ", line_break=False)
        plot_mism_types = self.create_mism_types_plot()
        hs.print_update("Done", with_time=False)

        hs.print_update(f"   - creating base composition plots at {category} and {corresponding_base} sites... ", line_break=False)
        plot_comp = self.create_composition_plot()
        hs.print_update("Done", with_time=False)

        hs.print_update(f"   - creating overview of error rates at {category} and {corresponding_base} sites... ", line_break=False)
        plot_err_rate, p_val_table = self.create_error_rate_plot()
        hs.print_update("Done", with_time=False)
        
        hs.print_update(f"   - creating overview of neighbouring mismatches around {category} sites... ", line_break=False)
        plot_nb = self.create_nb_plot()
        hs.print_update("Done", with_time=False)

        plots = [plot_mod_map, plot_mism_types, plot_comp, plot_err_rate, plot_nb]
        tables = [p_val_table]

        hs.print_update(f"  - creating HTML summary file at {output_path}... ", line_break=False)
        self.write_template(plots, tables, output_path)
        hs.print_update("Done", with_time=False)

        if self.output_tsv:
            hs.print_update(f"  - creating updated feature table at {self.in_path}... ", line_break=False)
            self.write_tsv()
            hs.print_update("Done", with_time=False)
        hs.print_update(f"Finished processing {category}.")

    ##############################################################################################################
    #                                          Data preparation methods                                          #
    ##############################################################################################################
        
    def subset_category_data(self, category: str, counterpart: str|None) -> None:
        """
        Subset the data based on the given category and counterpart (if provided).

        Parameters:
        - category (str): The category to filter the data by.
        - counterpart (str|None): The counterpart to filter the data by. If None, counterpart filtering is skipped.

        Returns:
        None

        The method filters the internal data based on the provided category and counterpart values,
        creating two subsets: one for the specified category and another for the counterpart (if provided).

        It updates the following attributes:
        - self.category_data (dict): A dictionary containing subsets of data for the specified category.
        - self.counterpart_data (dict): A dictionary containing subsets of data for the specified counterpart.
        - self.current_category (str): The currently selected category.
        - self.current_counterpart (str|None): The currently selected counterpart.
        """
        bed_name_idx = list(self.data.keys()).index("bed_name")
        ref_base_idx = list(self.data.keys()).index("ref_base")
        data_keys = self.data.keys()
        key_idx = dict(zip([i for i in range(len(data_keys))], data_keys))

        subset_category = dict(zip(data_keys, [[] for _ in range(len(data_keys))]))
        subset_counterpart = dict(zip(data_keys, [[] for _ in range(len(data_keys))]))
        
        for elements in zip(*(self.data[key] for key in data_keys)):
            if elements[bed_name_idx] == category:
                for i, element in enumerate(elements):
                    subset_category[key_idx[i]].append(element)
            if counterpart:
                if elements[ref_base_idx] == counterpart:
                    for i, element in enumerate(elements):
                        subset_counterpart[key_idx[i]].append(element)
        
        self.category_data = subset_category
        self.counterpart_data = subset_counterpart

        self.current_category = category
        self.current_counterpart = counterpart

    def prepare_data_map(self) -> Tuple[List[str], List[int]]:
        """
        Prepare the data for the sequence map based on references and category data.

        Returns:
        Tuple[List[str], List[int]]:
            A tuple containing two lists:
            - List of chromosome names present in both references and category data.
            - List of corresponding chromosome lengths.

        The method reads a fasta file specified by 'ref_path' and extracts the sequences
        along with their chromosome names. It then compares the unique chromosome names in
        the category data with those in the references, ensuring only common ones are considered.

        The resulting tuple contains:
        - List of common chromosome names sorted numerically and alphabetically.
        - List of corresponding lengths of the chromosomes.
        """
        def get_references(path: str) -> Dict[str, str]:
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
                refs = {}
                line = next(ref)
                if not line.startswith(">"):
                    raise Exception(f"Fasta format error. The first line of fasta file '{path}' does not contain a header (starting with '>').")
                
                chr_name = line[1:].strip().split(" ")[0]
                seq = ""
                for line in ref:
                    if line.startswith(">"):
                        refs[chr_name] = seq
                        chr_name = line[1:].strip().split(" ")[0]
                        seq = ""
                    else:
                        seq += line.strip()
                        
                refs[chr_name] = seq # add the last dict entry
                sys.stdout.write("\n")
            return refs

        def custom_sort_key(item):
            if item.isdigit():  # Check if the item is a digit
                return (int(item),)  # Convert to integer and sort numerically
            else:
                return (float('inf'), item)  # Place non-digits at the end

        ref_dict = get_references(self.ref_path)

        all_keys = ref_dict.keys()
        present_chr = list(sorted(set(self.category_data["chr"]).intersection(all_keys), key=custom_sort_key))
        chr_lens = [len(ref_dict[x]) for x in present_chr]

        return present_chr, chr_lens

    def prepare_data_mism_types(self) -> Tuple[List[List[int]], List[List[str]], int]:
        """
        Prepare data for a confusion matrix based on mismatch counts.

        Parameters:
        - mis_count (Dict[str, int]): Dictionary containing mismatch counts for each base pair combination.

        Returns:
        - Tuple[List[List[int]], List[List[str]], int]: Tuple containing matrix data, matrix labels, and the maximum value in the matrix.
        """

        mis_types = [f"{f} - {t}" for f in ["A", "C", "G", "U"] for t in ["A", "C", "G", "U"]]
        mis_count = dict(zip(mis_types, [0]*len(mis_types)))

        for ref, maj in zip(self.category_data["ref_base"], self.category_data["majority_base"]):
            if (ref != "N") & (maj != "N"):
                mis_count[f"{ref} - {maj}"] += 1

        # prepare data for the confusion matrix
        matrix_data = [[None]*4 for _ in range(4)]
        matrix_labels = [[None]*4 for _ in range(4)]
        bases = ["A", "C", "G", "U"]
        # fill count matrix
        for i in range(4): 
            for j in range(4):
                count = mis_count[f"{bases[i]} - {bases[j]}"]
                matrix_data[i][j] = count
                if i != j:
                    matrix_labels[i][j] = count
        # fill the matrix containing corresponding labels
        vals_flat = [element for sublist in matrix_labels for element in sublist if (element is not None)]
        vals_sum = sum(vals_flat)
        for i in range(4): 
            for j in range(4):
                if matrix_labels[i][j]:
                    matrix_labels[i][j] = f"{round(matrix_labels[i][j] / vals_sum * 100, 2)}%"
                else: 
                    matrix_labels[i][j] = ""
        # get the max value
        val_max = max(vals_flat)

        return matrix_data, matrix_labels, val_max

    def prepare_data_composition(self) -> Tuple[Dict[str, Tuple[float,float]|Tuple[float,float,float,float]], Dict[str, Tuple[float,float]|Tuple[float,float,float,float]], Tuple[int,int]|Tuple[int,int,int,int]]:
        """
        Prepare data for composition analysis, including base rates, medians, and scaled medians.

        Returns:
        Tuple[Dict[str, Tuple[float, float] | Tuple[float, float, float, float]],
            Dict[str, Tuple[float, float] | Tuple[float, float, float, float]],
            Tuple[int, int] | Tuple[int, int, int, int]]:
            A tuple containing three elements:
            1. Dictionary with scaled median rates for each base.
            2. Dictionary with un-scaled median rates for each base.
            3. Tuple with counts of matches and mismatches.

        The method processes base rates in the category data, calculates medians,
        and scales the medians to represent proportions. If a counterpart is specified,
        the same process is applied, and the results are combined for further analysis.

        The returned tuple includes:
        1. Scaled median rates for each base.
        2. Un-scaled median rates for each base.
        3. Counts of matches and mismatches (category-match, counterpart-match, category-mismatch, counterpart-mismatch).
        """
        cols = ["ref_base", "majority_base", "a_rate", "c_rate", "g_rate", "u_rate"]

        def get_rates(data: Dict[str, List[str|int|float]]) -> Tuple[Dict[str, Tuple[List[float], List[float]]], Tuple[int]]:
            rate_dict = defaultdict(lambda: ([], [])) # first sublist -> match; second sublist -> mismatch
            counts = [0, 0]
            for ref, maj, *base_rates in zip(*(data[col] for col in cols)):
                i = 0 if ref == maj else 1
                counts[i] += 1
                for j, base in enumerate(["A", "C", "G", "U"]):
                    rate_dict[base][i].append(base_rates[j])

            return dict(rate_dict), tuple(counts)

        def get_median_rates(rate_dict: Dict[str, Tuple[List[float], List[float]]]) -> Dict[str, Tuple[float, float]|Tuple[float,float,float,float]]:
            median_rate_dict = {}
            for key in rate_dict.keys():
                median_rate_dict[key] = (median(rate_dict[key][0]), median(rate_dict[key][1]))
            return median_rate_dict

        def get_scaled_median_rates(median_rate_dict: Dict[str, Tuple[float, float]|Tuple[float,float,float,float]]) -> Dict[str, Tuple[float, float]|Tuple[float,float,float,float]]:
            scaled_rate_dict = {}
            # calculate the sum of the A-, C-, G- and U-rates for each x-value
            totals = [sum(median_rates) for median_rates in zip(*(median_rate_dict[key] for key in median_rate_dict.keys()))]
            for key, medians in median_rate_dict.items():
                # scale each value by the corresponding x-value (rates should add up to 1)
                scaled_rate_dict[key] = tuple(med/total for med, total in zip(medians, totals))
            return scaled_rate_dict

        rates, counts = get_rates(self.category_data)
        median_rates = get_median_rates(rates)

        if self.current_counterpart:
            rates_cou, counts_cou = get_rates(self.counterpart_data)
            median_rates_cou = get_median_rates(rates_cou)
            # combine count_dict and count_dict_counterpart to have the x-axis order:
            # category-match, counterpart-match, category-mismatch, counterpart-mismatch 
            rate_tmp = {}
            for (base, rate), (_, rate_cou) in zip(median_rates.items(), median_rates_cou.items()):
                rate_tmp[base] = (rate[0], rate_cou[0], rate[1], rate_cou[1])
            median_rates = rate_tmp

            counts = (counts[0], counts_cou[0], counts[1], counts_cou[1])

            x_vals = [f"<i>{self.current_category}</i> match<br>(n = {counts[0]})", 
                    f"{self.current_counterpart} match<br>(n = {counts[1]})", 
                    f"<i>{self.current_category}</i> mismatch<br>(n = {counts[2]})", 
                    f"{self.current_counterpart} mismatch<br>(n = {counts[3]})"]
        else:
            x_vals = [f"<i>{self.current_category}</i> match<br>(n = {counts[0]})", 
                    f"<i>{self.current_category}</i> mismatch<br>(n = {counts[1]})"]

        return get_scaled_median_rates(median_rates), median_rates, x_vals
    
    def prepare_data_errorrates(self) -> Tuple[Dict[str, List[float]], 
                                        Dict[str, List[float]], 
                                        Dict[str, List[float]], 
                                        Dict[str, List[float]], 
                                        Dict[str, List[float]],
                                        List[str]]:
        """
        Prepare error rates data for analysis.

        Returns:
        Tuple[Dict[str, List[float]], 
            Dict[str, List[float]], 
            Dict[str, List[float]], 
            Dict[str, List[float]], 
            Dict[str, List[float]],
            List[str]]:
            A tuple containing five dictionaries representing different error rates and a list of x-axis labels.

        The method processes error rates in the category data, including mismatch,
        deletion rate, insertion rate, reference skip rate, and mean quality score.
        If a counterpart is specified, the same process is applied, and the results are combined for further analysis.

        The returned tuple includes:
        1. Percentage mismatch rates.
        2. Relative deletion rates.
        3. Relative insertion rates.
        4. Relative reference skip rates.
        5. Mean quality scores.
        6. List of x-axis labels representing different conditions.
        """
        def get_rates(data: Dict[str, List[str|int|float]], x_mat: int, x_mis: int) -> Dict[str, Dict[str, List[float]]]:
            rate_dict = defaultdict(lambda: {"x": [], "y": []})
            cols = ["ref_base", "majority_base", self.perc_mismatch_col, "deletion_rate", "insertion_rate", "refskip_rate", "q_mean"]
            cols_idxs = dict(zip(cols[2:], range(len(cols))))
            for ref, maj, *rates in zip(*(data[col] for col in cols)):
                x = x_mat if ref == maj else x_mis
                for rate in [self.perc_mismatch_col, "deletion_rate", "insertion_rate", "refskip_rate", "q_mean"]:
                    rate_dict[rate]["x"].append(x)
                    rate_dict[rate]["y"].append(rates[cols_idxs[rate]])
            return dict(rate_dict)

        x_mis = 2 if self.current_counterpart else 1
        rates = get_rates(self.category_data, x_mat=0, x_mis=x_mis)

        if self.current_counterpart:
            rates_cou = get_rates(self.counterpart_data, x_mat=1, x_mis=3)
            for key in rates.keys():
                rates[key]["x"] += rates_cou[key]["x"]
                rates[key]["y"] += rates_cou[key]["y"]
            x_vals = [f"<i>{self.current_category}</i> match", 
                    f"{self.current_counterpart} match", 
                    f"<i>{self.current_category}</i> mismatch", 
                    f"{self.current_counterpart} mismatch"]
        else:
            x_vals = [f"<i>{self.current_category}</i> match", 
                    f"<i>{self.current_category}</i> mismatch"]

        return *rates.values(), x_vals

    def prepare_nb_counts(self) -> Tuple[List[int], List[int], List[int], List[int], List[str], List[int]]:
        """
        Prepare data for analyzing neighboring error positions in a specified bed category.

        Returns:
            Tuple[List[int], List[int], List[int], List[int], List[str], List[int]]:
            A tuple containing:
            - x_vals (List[int]): Positions of neighboring errors.
            - y_vals (List[int]): Counts of neighboring errors.
            - x_vals_mod (List[int]): Positions of neighboring errors due to modifications.
            - y_vals_mod (List[int]): Counts of neighboring errors due to modifications.
            - pie_labs (List[str]): Labels for a pie chart indicating the presence of surrounding errors.
            - pie_vals (List[int]): Counts for the pie chart categories.

        Note:
            This method prepares data for analyzing neighboring error positions in a specified bed category.
            It extracts relevant data from the DataFrame, including neighboring error positions, and counts
            of these positions. Additionally, it identifies neighboring errors that are due to modifications.
            The resulting data is organized into lists for plotting, including a pie chart indicating the presence
            of surrounding errors.

        """
        def to_numeric(pos_str: str) -> List[int]:
            if len(pos_str) > 0:
                if pos_str.endswith(","): 
                    pos_str = pos_str[:-1] 
                return list(map(int, pos_str.split(",")))
            return []

        sites = [(c, s) for c, s in zip(self.category_data["chr"], self.category_data["site"])]

        n_no_nb = 0
        n_has_nb = 0
        n_has_nb_mod = 0

        counts = defaultdict(lambda: 0)
        counts_mod = defaultdict(lambda: 0)

        for site, nb_info in zip(sites, self.category_data["neighbour_error_pos"]):
            if nb_info:
                nb_info = to_numeric(nb_info)
                for distance in nb_info:
                    # check if pos + distance corresponds to another mod site
                    if (site[0], site[1]+distance) in sites:
                        n_has_nb_mod += 1
                        counts_mod[distance] += 1
                    else:
                        n_has_nb += 1
                        counts[distance] += 1
            else:
                n_no_nb += 1

        pie_labs = ["No surrounding errors", "Has surrounding errors", "Has surrounding error due to mod."]
        pie_vals = [n_no_nb, n_has_nb, n_has_nb_mod]

        return list(counts.keys()), list(counts.values()), list(counts_mod.keys()), list(counts_mod.values()), pie_labs, pie_vals

    ##############################################################################################################
    #                                            p-val table methods                                             #
    ##############################################################################################################
    def perform_test(self, set1: List[float], set2: List[float], alpha: float = 0.05) -> str:
        """
        Perform a statistical test to compare two sets of data.

        Parameters:
        - set1 (List[float]): The first set of data.
        - set2 (List[float]): The second set of data.
        - alpha (float): The significance level for the tests (default is 0.05).

        Returns:
        str:
            A string representing the result of the statistical test.

        The method performs a series of tests to compare two sets of data. It checks for normality
        in both samples and equal variances using normaltest and bartlett tests. If the conditions
        are met, it performs an independent t-test (TT), otherwise, it performs a Mann-Whitney U test (MWU).

        The returned string includes the p-value of the test and the abbreviation of the test method.
        In case of an exception during the test, it returns "ERROR".
        """
        try:
            if (normaltest(set1)[1] >= alpha) & (normaltest(set2)[1] >= alpha) & (bartlett(set1, set2)[1] >= alpha): # normal distributions in both samples + equal variances
                return f"{ttest_ind(a=set1, b=set2)[1]:.5e} (TT)"
            else: # normal distribution in both samples and equal variances
                return f"{mannwhitneyu(x=set1, y=set2)[1]:.5e} (MWU)"
        except:
            return "ERROR"

    def get_table_header(self) -> str:
        """
        Generate an HTML table header based on the current category and counterpart.

        Returns:
        str:
            A string containing the HTML table header.

        The method generates an HTML table header with column labels based on the current
        category and counterpart (if applicable). It returns a string formatted as an HTML table row.
        """
        if self.current_counterpart:
            header = f"<th></th><th>{self.current_category} match vs. {self.current_counterpart} match</th><th>{self.current_category} mismatch vs. {self.current_counterpart} mismatch</th><th>{self.current_category} match vs. {self.current_category} mismatch</th><th>{self.current_counterpart} match vs. {self.current_counterpart} mismatch</th>"
        else: 
            header = f"<th></th><th>{self.current_category} match vs. {self.current_category} mismatch</th>"
        return "<tr>" + header + "</tr>"

    # TODO: add a second table for statistical evaluation of the base compositions

    def create_test_overview(self, datasets: List[Dict[str, List[float]]]) -> str:
        """
        Create an HTML table summarizing statistical tests for multiple datasets.

        Parameters:
        - datasets (List[Dict[str, List[float]]]): A list of datasets, each containing 'x' and 'y' lists.

        Returns:
        str:
            A string containing the HTML table with test results.

        The method generates an HTML table summarizing statistical tests for multiple datasets.
        Each dataset should include 'x' and 'y' lists. The table includes rows for different data types
        (Mismatch rate, Deletion rate, Insertion rate, Reference skip rate, Mean quality) and columns
        for various comparisons based on the current category and counterpart.

        The returned string contains the complete HTML table.
        """
        table = f"<table>" + self.get_table_header()

        for data_type, data in zip(["Mismatch rate", "Deletion rate", "Insertion rate", "Reference skip rate", "Mean quality"], datasets):
        # split data by groups ("x" entry in data dict)
            row_data = defaultdict(lambda: [])
            for x, y in zip(data["x"], data["y"]):
                row_data[x].append(y)
            row_data = dict(row_data)

            if self.current_counterpart:
                row = f"<tr><td>{data_type}</td><td>{self.perform_test(row_data[0], row_data[1])}</td><td>{self.perform_test(row_data[2], row_data[3])}</td><td>{self.perform_test(row_data[0], row_data[2])}</td><td>{self.perform_test(row_data[1], row_data[3])}</td></tr>"
            else:
                row = f"<tr><td>{data_type}</td><td>{self.perform_test(row_data[0], row_data[1])}</td></tr>"

            table += row

        return table + "</table>"

    ##############################################################################################################
    #                                              Plotting methods                                              #
    ##############################################################################################################

    def update_plot(self, fig, title: str|None = None, xlab: str|None = None, ylab: str|None = None, height: int = 500, width: int = 800):
        """
        Updates the layout of a Plotly figure.

        Args:
            fig (go.Figure): The Plotly figure to be updated.
            title (str|None): Title of the plot (optional).
            xlab (str|None): Label for the x-axis (optional).
            ylab (str|None): Label for the y-axis (optional).
            height (int): Height of the plot (default: 500).
            width (int): Width of the plot (default: 800).

        Returns:
            go.Figure: The updated Plotly figure.
        """
        fig.update_layout(template="seaborn",
                    title = title,
                    xaxis_title = xlab,
                    yaxis_title = ylab,
                    font=dict(family="Open sans, sans-serif", size=20),
                    plot_bgcolor="white",
                    margin=dict(l=50, r=50, t=50, b=50),
                    height=height,  # Set the height to a constant value
                    width=width)
        fig.update_xaxes(showline=True, linewidth=2, linecolor='black', mirror=True, showticklabels=True, ticks='outside', showgrid=False, tickwidth=2)
        fig.update_yaxes(showline=True, linewidth=2, linecolor='black', mirror=True, showticklabels=True, ticks='outside', showgrid=False, tickwidth=2)

        return fig

    def create_map_plot(self) -> go.Figure: 
        """
        Create a map plot visualizing chromosome lengths and category sites.

        Returns:
        go.Figure:
            A Plotly figure representing the map plot.

        The method generates a map plot that visualizes chromosome lengths as bars and
        category sites as markers. The x-axis represents chromosomes, and the y-axis represents
        coordinates. The plot is customized for clear visualization, including bar width,
        scatter marker size, and layout adjustments.
        """
        present_chr, chr_lens = self.prepare_data_map()

        width = 1200
        bargap = 0.9
        bar_width = 1-bargap+0.15
        n_bars = len(present_chr)
        scatter_size = width/n_bars*bar_width

        fig = self.update_plot(go.Figure(), height = 1000, width = width, ylab="Coordinate")
        fig.add_trace(go.Bar(x=list(present_chr), y=list(chr_lens), marker=dict(color="lightgrey", line=dict(color="black", width=2)), name="Chromosomes", showlegend=False))
        fig.update_layout(bargap=0.5, yaxis=dict(range=[0,max(chr_lens)+0.1*max(chr_lens)]))
        fig.add_trace(go.Scatter(x=self.category_data["chr"], y=self.category_data["site"], mode='markers', marker=dict(symbol='line-ew', color="#1f77b4", size=scatter_size, line=dict(width=1.1, color="#1f77b4")), name=f"{self.current_category} sites", hovertemplate="Chr%{x}:%{y}"))
        fig.update_xaxes(fixedrange=True)

        return fig

    def create_mism_types_plot(self) -> go.Figure:
        """
        Create a heatmap plot visualizing mismatch types.

        Returns:
        go.Figure:
            A Plotly figure representing the mismatch types heatmap.

        The method generates a heatmap plot that visualizes mismatch types between called and reference bases.
        The x-axis represents the called bases, the y-axis represents the reference bases, and the color
        represents the count of occurrences for each mismatch type. The plot is customized for clear visualization.
        """
        matrix_data, matrix_labels, matrix_max_mism = self.prepare_data_mism_types()

        fig = px.imshow(matrix_data, labels=dict(x="Called base", y="Reference base", color="Count"), zmin=0, zmax=1.2*matrix_max_mism, color_continuous_scale="portland")
        fig = self.update_plot(fig, None, "Called base", "Reference base", width=800)
        fig.update_traces(text=matrix_labels, texttemplate="%{text}")
        fig.update_xaxes(fixedrange=True, tickvals=[0,1,2,3],ticktext=["A", "C", "G", "U"])
        fig.update_yaxes(fixedrange=True, tickvals=[0,1,2,3],ticktext=["A", "C", "G", "U"])
        
        return fig

    def create_composition_plot(self) -> go.Figure:
        """
        Create a stacked bar plot visualizing composition rates for different bases.

        Returns:
        go.Figure:
            A Plotly figure representing the composition rates stacked bar plot.

        The method generates a stacked bar plot that visualizes the composition rates
        for different bases in various conditions. The x-axis represents different conditions,
        and the y-axis represents the relative abundance of each base. The plot is customized
        for clear visualization, including colors, tooltips, and layout adjustments.
        """
        median_base_rates_scaled, median_base_rates, x_vals = self.prepare_data_composition()

        bar_a = go.Bar(x=x_vals, y=median_base_rates_scaled["A"], name="A", marker=dict(color="#2ca02c"), customdata=median_base_rates["A"], hovertemplate="Scaled median rate: %{y}<br>Unscaled median rate: %{customdata}")
        bar_c = go.Bar(x=x_vals, y=median_base_rates_scaled["C"], name="C", marker=dict(color="#1f77b4"), customdata=median_base_rates["C"], hovertemplate="Scaled median rate: %{y}<br>Unscaled median rate: %{customdata}")
        bar_g = go.Bar(x=x_vals, y=median_base_rates_scaled["G"], name="G", marker=dict(color="#ff7f0e"), customdata=median_base_rates["G"], hovertemplate="Scaled median rate: %{y}<br>Unscaled median rate: %{customdata}")
        bar_u = go.Bar(x=x_vals, y=median_base_rates_scaled["U"], name="U", marker=dict(color="#d62728"), customdata=median_base_rates["U"], hovertemplate="Scaled median rate: %{y}<br>Unscaled median rate: %{customdata}")
        
        fig = self.update_plot(go.Figure(), ylab="Relative abundance", height=600, width=1000)
        fig.add_traces([bar_a, bar_c, bar_g, bar_u])
        fig.update_layout(barmode="stack")
        fig.update_traces(marker=dict(line=dict(color="black", width=1.5)))

        return fig
    
    def create_error_rate_plot(self) -> Tuple[go.Figure, str]:
        """
        Create a grouped box plot visualizing error rates and quality scores.

        Returns:
        Tuple[go.Figure, str]:
            A tuple containing the Plotly figure representing the box plot and a string
            with a table summarizing statistical test results.

        The method generates a grouped box plot that visualizes error rates and quality scores
        for different conditions. The x-axis represents different conditions, and the y-axis represents
        the error rate or quality score. The plot includes statistical tests and a table summarizing
        p-values for various comparisons.
        """
        data_mismatch, data_deletion, data_insertion, data_ref_skip, data_quality, x_ticks = self.prepare_data_errorrates()

        p_val_table = self.create_test_overview([data_mismatch, data_deletion, data_insertion, data_ref_skip, data_quality])

        box_mismatch = go.Box(x=data_mismatch["x"], y=data_mismatch["y"], name="Mismatch rate", offsetgroup=0, line=dict(color="black"), marker=dict(outliercolor="black", size=2), fillcolor="#8c564b")
        box_deletion = go.Box(x=data_deletion["x"], y=data_deletion["y"], name="Deletion rate", offsetgroup=1, line=dict(color="black"), marker=dict(outliercolor="black", size=2), fillcolor="#e377c2")
        box_insertion = go.Box(x=data_insertion["x"], y=data_insertion["y"], name="Insertion rate", offsetgroup=2, line=dict(color="black"), marker=dict(outliercolor="black", size=2), fillcolor="#7f7f7f")
        box_ref_skip = go.Box(x=data_ref_skip["x"], y=data_ref_skip["y"], name="Reference skip", offsetgroup=3, line=dict(color="black"), marker=dict(outliercolor="black", size=2), fillcolor="#bcbd22")
        
        box_quality = go.Box(x=data_quality["x"], y=data_quality["y"], name="Mean quality", offsetgroup=0, line=dict(color="black"), marker=dict(outliercolor="black", size=2), fillcolor="#17becf")

        fig = self.update_plot(make_subplots(rows=1, cols=2, column_widths=[0.75, 0.25]), height=800, width=1200)
        fig.add_traces([box_mismatch, box_deletion, box_insertion, box_ref_skip], rows=[1,1,1,1], cols=[1,1,1,1])
        fig.add_trace(box_quality, row=1, col=2)
        fig.update_layout(boxmode="group")
        fig.update_yaxes(title_text="Error rate", row = 1, col = 1)
        fig.update_yaxes(title_text="Quality score", row = 1, col = 2)
        fig.update_xaxes(tickvals=[0,1,2,3],ticktext=x_ticks)

        return fig, p_val_table
    
    def create_nb_plot(self) -> go.Figure:
        """
        Create a neighbor position plot for a specified modification type.

        Returns:
            go.Figure: A Plotly figure representing the bar and pie plot.

        This method creates a neighbor position plot for a specified modification type. It prepares data for
        the plot, generates stacked bar charts for neighbor positions with and without modification errors,
        as well as a pie chart showing the distribution of surrounding errors, using Plotly, and returns the
        resulting plot as HTML code.

        """
        x_vals, y_vals, x_vals_mod, y_vals_mod, pie_labs, pie_vals = self.prepare_nb_counts()

        fig = self.update_plot(make_subplots(rows=1, cols=2, column_widths=[0.75, 0.25], specs=[[{"type": "bar"}, {"type": "pie"}]]), width=1200)

        fig.add_trace(go.Bar(x=x_vals, y=y_vals, marker=dict(color="#d62728", line=dict(color='#000000', width=1.5)), name="nb. error"), row=1, col=1)
        fig.add_trace(go.Bar(x=x_vals_mod, y=y_vals_mod, marker=dict(color="#dd8452", line=dict(color='#000000', width=1.5)), name="nb. mod error"), row=1, col=1)
        fig.update_layout(barmode="stack")

        fig.update_xaxes(title = f"Relative position to {self.current_category} positions", row=1, col=1)
        fig.update_yaxes(title = "Count", row=1, col=1)

        fig.add_trace(go.Pie(labels=pie_labs, values=pie_vals, name="", 
                            marker=dict(colors=["#4c72b0", "#d62728", "#dd8452"], line=dict(color='#000000', width=1.5))), row=1, col=2)

        fig.update_layout(showlegend=False)
        
        return fig    
    
    ##############################################################################################################
    #                                               Output methods                                               #
    ##############################################################################################################
    def write_svg(self, fig: go.Figure, name: str, output_path) -> None:
        """
        Write a Plotly figure to an SVG file.

        Parameters:
        - fig (go.Figure): Plotly figure to be saved.
        - name (str): Name to be used for the output SVG file.

        Returns:
        None

        The method takes a Plotly figure and a name as input and writes the figure to an SVG file.
        The output file will be saved with a name based on the provided 'name' parameter and the
        instance's 'output_path'.
        """
        outpath = f"{os.path.splitext(output_path)[0]}_{name}.svg"
        fig.write_image(outpath)

    def figs_to_str(self, plot_figs: List[go.Figure]) -> List[str]:
        """
        Convert a list of Plotly figures to a list of HTML strings.

        Parameters:
        - plot_figs (List[go.Figure]): List of Plotly figures to be converted.

        Returns:
        List[str]: List of HTML strings representing the Plotly figures.

        The method takes a list of Plotly figures as input and converts each figure
        to its corresponding HTML representation. The resulting list contains HTML strings
        for each figure, suitable for embedding in a web page or displaying in an HTML environment.
        """
        return list(map(lambda x: to_html(x, include_plotlyjs=False), plot_figs))

    def write_template(self, plot_figs: List[go.Figure], tables: List[str], output_path: str) -> None:
        """
        Write a summary HTML report with interactive plots and tables.

        Parameters:
        - plot_figs (List[go.Figure]): List of Plotly figures to be included in the report.
        - tables (List[str]): List of HTML tables to be included in the report.

        Returns:
        None

        The method takes a list of Plotly figures and HTML tables, and generates an HTML report
        containing interactive plots and tables. The report includes sections for mapping positions,
        mismatch types, base compositions, error rates, and neighboring errors. If specified, SVG files
        of individual plots are also saved.
        """
        name = f"<i>{self.current_category}</i>"
        time = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')

        if self.export_svg:
            for fig, name in zip(plot_figs, ["poi_map", "poi_mismatch_types", "poi_base_comp", 
                                             "poi_error_rates", "poi_neighbours"]):
                self.write_svg(fig, name, output_path)

        plots = self.figs_to_str(plot_figs)

        css_string, plotly_js_string = hs.load_html_template_str()
        
        template = f"""
            <!DOCTYPE html>
            <html lang="en">
            <head>
                <meta charset="UTF-8">
                <meta http-equiv="X-UA-Compatible" content="IE=edge">
                <meta name="viewport" content="width=device-width, initial-scale=1.0">
                <title>Neet - {self.current_category} summary</title>                        
                <style>{css_string}</style>
            </head>

            <body>
                <script>{plotly_js_string}</script>

                <header>
                    <h1>Positions of interest: {name}</h1>
                    <p>Produced by <a href="https://github.com/dietvin/neet">Neet</a> on <b>{time}</b></p>
                </header>
                
                <section>
                    <p class="intro-text">
                        This summary file was created from the extracted features in file <b>{self.in_path}</b> 
                        with <b>{self.current_category}</b> positions extracted from file <b>{self.bed_path}</b>. 
                        The plots are interactive and allow further information by hovering, zooming and panning.
                    </p>
                </section>

                
                <section>
                    <button class="collapsible">{name} positions on reference sequence(s)</button>
                    <div class="collapsible-content">
                        <h2 class="hiddentitle" id="map"></h2>

                        <h3>Mapping of {name} positions across the reference sequences</h3>
                        <div class="plot-container">
                            {plots[0]}
                        </div>
                        <p>
                            Positions of {name} are indicated as blue lines. Fasta file '{self.ref_path}' was used to extract
                            reference sequence(s). Hover on positions for exact coordinates. 
                        </p>
                    </div>
                </section>

                <section>
                    <button class="collapsible">Mismatch types</button>
                    <div class="collapsible-content">
                        <h2 class="hiddentitle" id="mismatch_types"></h2>
                        
                        <h3>Confusion matrix of {name} positions containing mismatch types</h3>
                        <div class="plot-container">
                            {plots[1]}
                        </div>
                        <p>
                            Overview of all types of (mis-)matches in the data subset corresponding to {name}.
                        </p>
                    </div>
                </section>

                <section>
                    <button class="collapsible">Base compositions</button>
                    <div class="collapsible-content">
                        <h2 class="hiddentitle" id="base_comp"></h2>

                        <h3>Base compositions for different {name} {f"and {self.current_counterpart} " if self.current_counterpart else ""}subsets</h3>
                        <div class="plot-container">
                            {plots[2]}
                        </div>
                        <p>
                            Each A/C/G/U element in the bars corresponds to the median count of a given base in a subset. The four medians were scaled to add up to one. 
                            {name} match: positions labelled {name} from the bed file, where the called base is equal to the reference base. 
                            {f"{self.current_counterpart} match: positions with reference base {self.current_counterpart}, where the called base is equal." if self.current_counterpart else ""}
                            {name} mismatch: positions labelled {name} from the bed file, where the called base differs from the reference base. 
                            {f"{self.current_counterpart} mismatch: positions with reference base {self.current_counterpart} , where the called base differs." if self.current_counterpart else ""}
                        </p>
                    </div>
                </section>

                <section>
                    <button class="collapsible">Error rates</button>
                    <div class="collapsible-content">
                        <h2 class="hiddentitle" id="error_rates"></h2>

                        <h3>Error rates for different {name} {f"and {self.current_counterpart} " if self.current_counterpart else ""}subsets</h3>
                        <div class="plot-container">
                            {plots[3]}
                        </div>
                        
                        <p>
                            Left: Distributions of mismatch, deletion, insertion and reference skip rates for different subsets. Right: Distribution of mean quality scores for different subsets. 
                            {name} match: positions labelled {name} from the bed file, where the called base is equal to the reference base. 
                            {f"{self.current_counterpart} match: positions with reference base {self.current_counterpart}, where the called base is equal." if self.current_counterpart else ""}
                            {name} mismatch: positions labelled {name} from the bed file, where the called base differs from the reference base. 
                            {f"{self.current_counterpart} mismatch: positions with reference base {self.current_counterpart} , where the called base differs." if self.current_counterpart else ""}
                        </p>
                        <p>
                            The table below shows p-values for statistical tests between the four groups for different error rate
                            the mean quality distributions. 
                        </p>
                        <div class="table-box">
                            {tables[0]}
                        </div>

                    </div>
                </section>

                <section>
                    <button class="collapsible">Neighbouring errors</button>
                    <div class="collapsible-content">
                        <h2 class="hiddentitle" id="nb_errors"></h2>

                        <h3>Count of positions with high error rate in the surrounding of {name} positions</h3>
                        <div class="plot-container">
                            {plots[4]}
                        </div>
                        <p>
                            Left: Occurences of high mismatch rates two bases up- and downstream from {name} positions.
                            Right: Pie chart shows the (relative count) of different types of central {name} positions. 
                            Red indicates errors in the surrounding positions where the position in question
                            is also of type {name}. The count of surrounding error positions that do not fall under {name}
                            are colored orange. Blue corresponds to {name} positions where no surrounding errors are found.
                        </p>

                    </div>
                </section>

                <script>
                    var coll = document.getElementsByClassName("collapsible");
                    var i;

                    for (i = 0; i < coll.length; i++) {{
                    coll[i].addEventListener("click", function() {{
                        this.classList.toggle("active");
                        var content = this.nextElementSibling;
                        if (content.style.display === "none") {{
                        content.style.display = "block";
                        }} else {{
                        content.style.display = "none";
                        }}
                    }});
                    }}
                </script>

            </body>
            <footer></footer>
            </html> 
        """
        with open(output_path, "w") as out:
            out.write(template)

    def write_tsv(self) -> None:
        """
        Write the processed data to a tab-separated values (TSV) file.

        The method writes the processed data, stored in the class attribute 'data', to a TSV file.
        The file will be named based on the input file path, and the suffix '_w_bed_info.tsv' will
        be added to indicate that bed information has been appended to the original data.

        Parameters:
        None

        Returns:
        None

        Example:
        If self.in_path = "input_folder/example.txt", the output TSV file will be named
        "input_folder/example_w_bed_info.tsv".
        """
        with open(f"{os.path.splitext(self.in_path)[0]}_w_bed_info.tsv", "w") as out:
            out.write("\t".join(self.data.keys())+"\n")
            for vals in zip(*self.data.values()):
                vals = [str(val) for val in vals]
                out.write("\t".join(vals)+"\n")