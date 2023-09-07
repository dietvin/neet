from typing import List, Dict, Tuple, Any
import os, sys, warnings
import numpy as np
import pandas as pd
import plotly.graph_objects as go
import plotly.express as px
from plotly.subplots import make_subplots
from plotly.io import to_html
from intervaltree import IntervalTree
import collections, datetime, argparse

class POIAnalyzer():

    in_path: str
    bed_path: str
    ref_path: str
    output_path: str

    perc_mismatch_col: str
    output_tsv: bool

    data: pd.DataFrame
    bed_categories: List[str]
    bed_categories_canoncial_bases: List[str]

    dtypes = {'chr': str, 'site': int, 'n_reads': int, 'ref_base': str, 'majority_base': str, 'n_a': int, 'n_c': int,
            'n_g': int, 'n_t': int, 'n_del': int, 'n_ins': int, 'n_ref_skip': int, 'n_a_rel': float, 'n_c_rel': float,
            'n_g_rel': float, 'n_t_rel': float, 'n_del_rel': float, 'n_ins_rel': float, 'n_ref_skip_rel': float,
            'perc_mismatch': float, 'perc_mismatch_alt': float, 'motif': str, 'q_mean': float, 'q_std': float,
            'neighbour_error_pos': str}
    
    def __init__(self, in_path: str, out_path: str, bed_path: str, ref_path: str,
                 categories: str,
                 canonical_counterpart: str,
                 output_tsv: bool = True, 
                 use_perc_mismatch_alt: bool = False) -> None:
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
            use_perc_mismatch_alt (bool, optional): Whether to use 'perc_mismatch_alt' instead of 'perc_mismatch' column. Default is False.

        Returns:
            None

        Note:
            This constructor initializes an instance of the GenomicDataProcessor class. It sets various attributes
            based on the provided arguments and performs necessary data processing and validation steps.

        """
        self.process_path(in_path, out_path, bed_path, ref_path)
        self.load_data()
        self.perc_mismatch_col = "perc_mismatch_alt" if use_perc_mismatch_alt else "perc_mismatch"
        self.get_bed_categories(categories)
        self.get_corrensponding_base(canonical_counterpart)
        self.output_tsv = output_tsv

    ##############################################################################################################
    #                                           Initialization methods                                           #
    ##############################################################################################################

    def process_path(self, in_path: str, out_path: str, bed_path: str, ref_path: str) -> None:
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
        self.check_path(in_path, [".tsv"])
        self.in_path = in_path
        self.check_path(bed_path, [".bed"])
        self.bed_path = bed_path
        self.check_path(ref_path, [".fasta", ".fa", ".fn"])
        self.ref_path = ref_path
        self.output_path = self.process_outpath(out_path)

    def check_path(self, path: str, extensions: List[str]) -> None:
        """
        Check the existence and file extension of a specified path.

        Args:
            path (str): The path to be checked for existence and file extension.
            extensions (List[str]): A list of valid file extensions.

        Returns:
            None

        Raises:
            FileNotFoundError: If the specified 'path' does not exist.
            Warning: If the file extension of 'path' is not in the list of valid 'extensions'. The warning is issued
                    if the extension does not match but execution can continue.

        Note:
            This method checks whether the specified 'path' exists. If it does not exist, a FileNotFoundError
            is raised. It also checks the file extension of 'path' and compares it against the list of valid
            'extensions'. If the file extension does not match any of the valid extensions, a warning is issued.
            The warning is intended to notify the user if the file extension does not match expectations but
            allows for continued execution if this is deliberate.

        """
        if not os.path.exists(path): # does file exist?
            raise FileNotFoundError(f"Input file not found. File '{path}' does not exist.")
        file_type = os.path.splitext(path)[1]
        if not file_type in extensions:
            warnings.warn(f"Found file extension {file_type}. Expected file extension to be one of: {extensions}. If this is deliberate, ignore warning.", Warning)

    def process_outpath(self, out: str) -> str:
        """
        Process the provided output path, ensuring it is a valid directory path.

        Args:
            out (str): The output path to be processed.

        Returns:
            str: The processed output path, guaranteed to be a valid directory path.

        Raises:
            Exception: If the provided path is an existing file or if there is an error creating the directory.

        Note:
            This method checks if the provided path exists as a file or directory. If it exists as a file,
            it raises an exception, as a directory path is required. If it exists as a directory, it ensures
            that the path ends with a forward slash ('/') for consistency. If the path doesn't exist, it
            attempts to create the directory, raising an exception if it fails.

        """
        if os.path.isfile(out):
            raise Exception(f"Provided path '{out}' is a file. Please specify a directory as output path.")
        elif os.path.isdir(out):
            if not out.endswith("/"):
                out += "/"
        else:
            try: 
                os.makedirs(out)
            except Exception as e:
                raise Exception(f"Could not create directory '{out}'. Error: {e}")
        return out

    def load_data(self) -> None:
        """
        Load data from a CSV file with specified data types and add BED information.

        Reads a CSV file from the specified input path using pandas. The data is loaded with tab ('\t') as
        the separator, and custom data types specified in 'dtypes' attribute are applied to columns. After
        loading the data, it can optionally be enriched with information from a BED file specified by 'bed_path'
        attribute.

        Returns:
            None

        Note:
            This method reads data from a CSV file located at 'in_path' attribute using pandas library.
            It also applies custom data types from the 'dtypes' attribute to columns during data loading.
            If 'bed_path' is provided and not None, it invokes the 'add_bed_info' method to add additional
            information from a BED file to the loaded data. This method operates in-place and updates
            the 'data' attribute of the class with the loaded and potentially enriched data.

        """
        data = pd.read_csv(self.in_path, sep="\t", dtype=self.dtypes)
        self.data = self.add_bed_info(data, self.bed_path)
    
    def add_bed_info(self, data: pd.DataFrame, bed_path: str):
        """
        Add BED information to a DataFrame based on a BED file.

        Args:
            data (pd.DataFrame): The DataFrame containing genomic data to enrich with BED information.
            bed_path (str): The path to the BED file containing genomic intervals and associated names.

        Returns:
            pd.DataFrame: The DataFrame with an additional 'bed_name' column containing BED information.

        Note:
            This method takes a DataFrame 'data' containing genomic data and adds an additional column 'bed_name'
            to it. The 'bed_name' column is populated by querying a previously built IntervalTree from the
            specified 'bed_path' for each row in the DataFrame. The method returns the enriched DataFrame.

        """
        tree = self.build_interval_tree(bed_path)
        data["bed_name"] = data.apply(lambda x: self.get_name_for_position((x[0], x[1]), tree), axis=1)
        return data

    def get_name_for_position(self, position, interval_tree) -> str|None:
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
        results = interval_tree[position[1]:position[1]+1]  # Query the interval tree
        for interval in results:
            chromosome, name = interval.data
            if chromosome == position[0]:
                return name
        return None

    def build_interval_tree(self, bed_file):
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

    def get_bed_categories(self, cat_str: str) -> None:
        """
        Extract and validate bed categories from a comma-separated string.

        Args:
            cat_str (str): A comma-separated string containing bed categories.

        Returns:
            None

        Raises:
            Exception: If any category in 'cat_str' is not found in the bed file.

        Note:
            This method parses a comma-separated string 'cat_str' containing bed categories. It ensures
            that each category in the string exists in the 'bed_name' column of the loaded genomic data
            ('data' attribute). If any category is not found, an exception is raised. If all categories
            are valid, they are stored in the 'bed_categories' attribute for later use.

        """
        categories = cat_str.split(",")
        unique_cat = list(self.data.bed_name.unique())

        for category in categories:
            if category not in unique_cat:
                raise Exception(f"Given name '{category}' was not found in the bed file.")
    
        self.bed_categories = categories

    def get_corrensponding_base(self, base_str: str) -> None:
        """
        Extract and validate corresponding bases for bed categories.

        Args:
            base_str (str): A comma-separated string containing bases corresponding to bed categories.

        Returns:
            None

        Raises:
            Exception: If the number of bases does not match the number of bed categories.

        Note:
            This method parses a comma-separated string 'base_str' containing bases corresponding to bed categories.
            It validates that the number of bases matches the number of bed categories previously set in the
            'bed_categories' attribute. If the counts do not match, an exception is raised. If the counts match,
            the corresponding bases are stored in the 'bed_categories_canonical_bases' attribute for later use.

        """        
        counterparts = base_str.split(",")
        if len(counterparts) != len(self.bed_categories):
            raise Exception(f"For the {len(self.bed_categories)} categories {len(counterparts)} corresponding bases were given. Each category must have a base it corresponds to.")
        self.bed_categories_canoncial_bases = counterparts

    ##############################################################################################################
    #                                          Data preparation methods                                          #
    ##############################################################################################################

    def prepare_data_mism_types(self, data: pd.DataFrame):
        """
        Prepare data for mismatch types analysis.

        Args:
            data (pd.DataFrame): The DataFrame containing genomic data.

        Returns:
            Tuple[pd.DataFrame, pd.DataFrame, float]: A tuple containing:
                - matrix_data (pd.DataFrame): A matrix representing mismatch types.
                - matrix_labels (pd.DataFrame): Labels for the mismatch types matrix.
                - matrix_max_mism (float): The maximum value in the mismatch types matrix.

        Note:
            This method prepares data for mismatch types analysis. It filters data to include only rows where
            both 'ref_base' and 'majority_base' are not 'N'. It then computes a cross-tabulation matrix ('matrix_data')
            of 'ref_base' and 'majority_base' columns. Labels for the matrix are stored in 'matrix_labels', with
            diagonal elements set to None. The 'matrix_max_mism' value is the maximum value in the matrix.
            The matrix values are also transformed into percentages.

        """
        matrix_data = data.loc[(data["ref_base"] != "N") & (data["majority_base"] != "N"),["ref_base", "majority_base"]]
        matrix_data = pd.crosstab(matrix_data['ref_base'], matrix_data['majority_base'])
        matrix_labels = matrix_data.copy()
        for i in range(matrix_data.shape[0]):
            matrix_labels.iloc[i,i] = None
        matrix_max_mism = matrix_labels.max().max()
        matrix_labels = np.round(matrix_labels / matrix_labels.sum().sum() * 100, 2).astype(str)
        matrix_labels = (matrix_labels + "%").replace("nan%", "")

        return matrix_data, matrix_labels, matrix_max_mism

    def prepare_data_composition(self, data: pd.DataFrame, mod: str, canonical: str) -> Tuple[Dict[str, List[float]], Dict[str, int]]:
        """
        Prepare data for composition analysis based on modification and canonical sequences.

        Args:
            data (pd.DataFrame): The DataFrame containing genomic data.
            mod (str): The name of the modified sequence category.
            canonical (str): The canonical sequence category.

        Returns:
            Tuple[Dict[str, List[float]], Dict[str, int]]: A tuple containing:
                - composition_data (Dict[str, List[float]]): A dictionary containing composition values for A, C, G, and T.
                - composition_counts (Dict[str, int]): A dictionary containing counts of modification and canonical sequences for A, C, G, and T.

        Note:
            This method prepares data for composition analysis based on modification and canonical sequences.
            It filters data to extract relevant sequences for 'mod' and 'canonical' categories, both matching
            and mismatching. Median values are computed for each category, and the relative composition of A, C, G, and T
            is calculated. The resulting data is returned in dictionaries for composition values and sequence counts.

        """
        mod_match = data.loc[(data.bed_name==mod) & (data.ref_base==data.majority_base), ["n_a_rel", "n_c_rel", "n_g_rel", "n_t_rel"]]
        mod_mismatch = data.loc[(data.bed_name==mod) & (data.ref_base!=data.majority_base), ["n_a_rel", "n_c_rel", "n_g_rel", "n_t_rel"]]
        unm_match = data.loc[(data.ref_base == canonical) & (data.bed_name.isna()) & (data.ref_base==data.majority_base), ["n_a_rel", "n_c_rel", "n_g_rel", "n_t_rel"]]
        unm_mismatch = data.loc[(data.ref_base == canonical) & (data.bed_name.isna()) & (data.ref_base!=data.majority_base), ["n_a_rel", "n_c_rel", "n_g_rel", "n_t_rel"]]

        n_mod_match = mod_match.shape[0]
        n_mod_mismatch = mod_mismatch.shape[0]
        n_unm_match = unm_match.shape[0]
        n_unm_mismatch = unm_mismatch.shape[0]

        mod_match = mod_match.median() 
        mod_mismatch = mod_mismatch.median()
        unm_match = unm_match.median()
        unm_mismatch = unm_mismatch.median()

        mod_match = mod_match / sum(mod_match)
        mod_mismatch = mod_mismatch / sum(mod_mismatch)
        unm_match = unm_match / sum(unm_match)
        unm_mismatch = unm_mismatch / sum(unm_mismatch)

        a_vals = [mod_match["n_a_rel"], unm_match["n_a_rel"], mod_mismatch["n_a_rel"], unm_mismatch["n_a_rel"]]
        c_vals = [mod_match["n_c_rel"], unm_match["n_c_rel"], mod_mismatch["n_c_rel"], unm_mismatch["n_c_rel"]]
        g_vals = [mod_match["n_g_rel"], unm_match["n_g_rel"], mod_mismatch["n_g_rel"], unm_mismatch["n_g_rel"]]
        t_vals = [mod_match["n_t_rel"], unm_match["n_t_rel"], mod_mismatch["n_t_rel"], unm_mismatch["n_t_rel"]]

        return {"A": a_vals, "C": c_vals, "G": g_vals, "T": t_vals}, {"A": n_mod_match, "C": n_mod_mismatch, "G": n_unm_match, "T": n_unm_mismatch}
    
    def prepare_data_errorrates(self, data: pd.DataFrame, mod: str, canonical: str) -> Tuple[Tuple[List[float], List[float]], 
                                                                                    Tuple[List[float], List[float]], 
                                                                                    Tuple[List[float], List[float]], 
                                                                                    Tuple[List[float], List[float]], 
                                                                                    Tuple[List[float], List[float]]]:
        """
        Prepare data for error rates analysis based on modification and canonical sequences.

        Args:
            data (pd.DataFrame): The DataFrame containing genomic data.
            mod (str): The name of the modified sequence category.
            canonical (str): The canonical sequence category.

        Returns:
            Tuple[Tuple[List[float], List[float]],
                Tuple[List[float], List[float]],
                Tuple[List[float], List[float]],
                Tuple[List[float], List[float]],
                Tuple[List[float], List[float]]]:
            A tuple containing five sub-tuples, each with two lists:
            - data_mismatch (Tuple[List[float], List[float]]): Mismatch percentage data for different categories.
            - data_deletion (Tuple[List[float], List[float]]): Deletion percentage data for different categories.
            - data_insertion (Tuple[List[float], List[float]]): Insertion percentage data for different categories.
            - data_ref_skip (Tuple[List[float], List[float]]): Reference skip percentage data for different categories.
            - data_quality (Tuple[List[float], List[float]]): Quality mean data for different categories.

        Note:
            This method prepares data for error rates analysis based on modification and canonical sequences.
            It filters data to extract relevant sequences for 'mod' and 'canonical' categories, both matching
            and mismatching. Various error rate-related columns are selected for analysis, and data is organized
            into sub-tuples containing lists for each category. The resulting data is returned in the specified
            tuple format for plotting.

        """
        cols = [self.perc_mismatch_col, "n_del_rel", "n_ins_rel", "n_ref_skip_rel", "q_mean"]
        mod_match = data.loc[(data.bed_name==mod) & (data.ref_base==data.majority_base), cols]
        unm_match = data.loc[(data.ref_base == canonical) & (data.bed_name.isna()) & (data.ref_base==data.majority_base), cols]
        mod_mismatch = data.loc[(data.bed_name==mod) & (data.ref_base!=data.majority_base), cols]
        unm_mismatch = data.loc[(data.ref_base == canonical) & (data.bed_name.isna()) & (data.ref_base!=data.majority_base), cols]
        
        datasets = [mod_match, unm_match, mod_mismatch, unm_mismatch]
        x_values = [f"<i>{mod}</i> match", f"{canonical} match", f"<i>{mod}</i> mismatch", f"{canonical} mismatch"]

        def get_x_y_vals(col: str) -> Tuple[List[float], List[float]]: 
            y_vals = []
            x_vals = []
            for dataset, x_val in zip(datasets, x_values):
                y = list(dataset[col])
                x = [x_val] * len(y)
                y_vals += y
                x_vals += x
            return x_vals, y_vals

        data_mismatch = get_x_y_vals(self.perc_mismatch_col)
        data_deletion = get_x_y_vals("n_del_rel")
        data_insertion = get_x_y_vals("n_ins_rel")
        data_ref_skip = get_x_y_vals("n_ref_skip_rel")
        data_quality = get_x_y_vals("q_mean")
        return (data_mismatch, data_deletion, data_insertion, data_ref_skip, data_quality)

    def prepare_nb_counts(self, data: pd.DataFrame, bed_category: str) -> Tuple[List[int], List[int], List[int], List[int], List[str], List[int]]:
        """
        Prepare data for analyzing neighboring error positions in a specified bed category.

        Args:
            data (pd.DataFrame): The DataFrame containing genomic data.
            bed_category (str): The bed category for which neighboring error positions should be analyzed.

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
            if type(pos_str) == float: # in case the value is nan
                return []
            else:
                if pos_str.endswith(","): 
                    pos_str = pos_str[:-1] 
                return list(map(int, pos_str.split(",")))

        def extract_neighbour_mod(df: pd.DataFrame) -> Tuple[List[int], List[int]]:
            """
            Receives a pd.DataFrame containing columns chr, site and neighbour_error_pos (with neighbour errors in List[int] format)
            For each row, check each in the neighbour_error_pos list if the corresponding coordinate is also in a modified one. 
            E.g.: At chr1:215 [-1,2] --> check if chr1:214 and chr1:217 is also a modified position (just checking the same mod. type)
            """
            positions = []
            positions_mod = []
            for _, row in df.iterrows():
                chrom = row[0]
                site = row[1]
                distances = row[2]
                for distance in distances:
                    if ((df.chr==chrom) & (df.site == site+distance)).any():
                        positions.append(distance)
                    else:
                        positions_mod.append(distance)
            return positions, positions_mod

        
        subset = data.loc[data["bed_name"]==bed_category, ["chr", "site", "neighbour_error_pos"]]
        n_no_nb = subset["neighbour_error_pos"].isna().sum()

        subset["neighbour_error_pos"] = subset["neighbour_error_pos"].apply(to_numeric)
        dist_vals, dist_vals_mod = extract_neighbour_mod(subset)
        
        value_counts = dict(collections.Counter(dist_vals))
        value_counts_mod = dict(collections.Counter(dist_vals_mod))

        x_vals = list(value_counts.keys())
        y_vals = list(value_counts.values())
        n_has_nb = sum(y_vals)

        x_vals_mod = list(value_counts_mod.keys())
        y_vals_mod = list(value_counts_mod.values())
        n_has_nb_mod = sum(y_vals_mod)

        
        pie_labs = ["No surrounding errors", "Has surrounding errors", "Has surrounding error due to mod."]
        pie_vals = [n_no_nb, n_has_nb, n_has_nb_mod]
        
        return x_vals, y_vals, x_vals_mod, y_vals_mod, pie_labs, pie_vals

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

    def create_map_plot(self, mod_type: str): 
        """
        Create a map plot for a specified modification type.

        Args:
            mod_type (str): The modification type for which the map plot should be created.

        Returns:
            str: HTML code representing the map plot.

        Note:
            This method creates a map plot for a specified modification type. It reads reference sequences from a FASTA file,
            extracts relevant data from the DataFrame, and plots chromosome lengths as bars and modification sites as markers.
            The resulting map plot is converted to HTML code.

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

            def stdout_progress(n: int):
                sys.stdout.write(f"\rSequences found: {n}")
                sys.stdout.flush()

            with open(path, "r") as ref:
                refs = {}
                line = next(ref)
                if not line.startswith(">"):
                    raise Exception(f"Fasta format error. The first line of fasta file '{path}' does not contain a header (starting with '>').")
                
                chr_name = line[1:].strip().split(" ")[0]
                seq = ""

                chr_found = 1
                stdout_progress(chr_found)

                for line in ref:
                    if line.startswith(">"):
                        refs[chr_name] = seq
                        chr_name = line[1:].strip().split(" ")[0]
                        seq = ""

                        chr_found += 1
                        stdout_progress(chr_found)

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

        data = self.data.loc[self.data.bed_name==mod_type] 
        ref_dict = get_references(self.ref_path)

        all_keys = ref_dict.keys()
        present_chr = set(all_keys & data["chr"].unique())
        present_chr = np.array(sorted(present_chr, key=custom_sort_key))
        chr_lens = np.array([len(ref_dict[x]) for x in present_chr])

        width = 1200
        bargap = 0.9
        bar_width = 1-bargap+0.15
        n_bars = len(present_chr)
        scatter_size = width/n_bars*bar_width

        fig = self.update_plot(go.Figure(), height = 1000, width = width, ylab="Coordinate")
        fig.add_trace(go.Bar(x=list(present_chr), y=list(chr_lens), marker=dict(color="lightgrey", line=dict(color="black", width=2)), name="Chromosomes", showlegend=False))
        fig.update_layout(bargap=0.5, yaxis=dict(range=[0,max(chr_lens)+0.1*max(chr_lens)]))

        fig.add_trace(go.Scatter(x=list(data["chr"]), y=list(data["site"]), mode='markers', marker=dict(symbol='line-ew', color="#1f77b4", size=scatter_size, line=dict(width=1.1, color="#1f77b4")), name=f"psU sites", hovertemplate="Chr%{x}:%{y}"))
        fig.update_xaxes(fixedrange=True)

        return to_html(fig, include_plotlyjs=False)

    def create_mism_types_plot(self, mod: str):
        """
        Create a mismatch types plot for a specified modification type.

        Args:
            mod (str): The modification type for which the mismatch types plot should be created.

        Returns:
            str: HTML code representing the mismatch types plot.

        Note:
            This method creates a mismatch types plot for a specified modification type. It prepares data
            for the plot, generates a heatmap using Plotly Express, and returns the resulting plot as HTML code.

        """
        data = self.data.loc[self.data.bed_name == mod]
        data.to_csv("/home/vincent/masterthesis/data/45s_rrna/processed/dRNA_cytoplasm/test.csv")
        matrix_data, matrix_labels, matrix_max_mism = self.prepare_data_mism_types(data)

        fig = px.imshow(matrix_data, labels=dict(x="Called base", y="Reference base", color="Count"), zmin=0, zmax=1.2*matrix_max_mism, color_continuous_scale="portland")
        fig = self.update_plot(fig, None, "Called base", "Reference base", width=800)
        fig.update_traces(text=matrix_labels, texttemplate="%{text}")
        fig.update_xaxes(fixedrange=True)
        fig.update_yaxes(fixedrange=True)
        
        return to_html(fig, include_plotlyjs=False)

    def create_composition_plot(self, mod: str, canonical: str) -> str:
        """
        Create a composition plot for a specified modification type and canonical sequence.

        Args:
            mod (str): The modification type for which the composition plot should be created.
            canonical (str): The canonical sequence for comparison.

        Returns:
            str: HTML code representing the composition plot.

        Note:
            This method creates a composition plot for a specified modification type and canonical sequence.
            It prepares data for the plot, generates a stacked bar chart using Plotly, and returns the resulting
            plot as HTML code.

        """
        y_data, category_count = self.prepare_data_composition(self.data, mod, canonical)

        x_vals = [f"<i>{mod}</i> match<br>(n = {category_count['A']})", 
                f"{canonical} match<br>(n = {category_count['G']})", 
                f"<i>{mod}</i> mismatch<br>(n = {category_count['C']})", 
                f"{canonical} mismatch<br>(n = {category_count['T']})"]

        bar_a = go.Bar(x=x_vals, y=y_data["A"], name="A", marker=dict(color="#2ca02c"), hovertemplate="%{y}")
        bar_c = go.Bar(x=x_vals, y=y_data["C"], name="C", marker=dict(color="#1f77b4"), hovertemplate="%{y}")
        bar_g = go.Bar(x=x_vals, y=y_data["G"], name="G", marker=dict(color="#ff7f0e"), hovertemplate="%{y}")
        bar_t = go.Bar(x=x_vals, y=y_data["T"], name="T", marker=dict(color="#d62728"), hovertemplate="%{y}")
        
        fig = self.update_plot(go.Figure(), ylab="Relative abundance", height=600, width=1000)
        fig.add_traces([bar_a, bar_c, bar_g, bar_t])
        fig.update_layout(barmode="stack")
        fig.update_traces(marker=dict(line=dict(color="black", width=1.5)))

        return to_html(fig, include_plotlyjs=False)

    def create_error_rate_plot(self, mod: str, canonical: str) -> str:
        """
        Create an error rate plot for a specified modification type and canonical sequence.

        Args:
            mod (str): The modification type for which the error rate plot should be created.
            canonical (str): The canonical sequence for comparison.

        Returns:
            str: HTML code representing the error rate plot.

        Note:
            This method creates an error rate plot for a specified modification type and canonical sequence.
            It prepares data for the plot, generates box plots for mismatch rate, deletion rate, insertion rate,
            reference skip rate, and mean quality score using Plotly, and returns the resulting plot as HTML code.

        """
        data_mismatch, data_deletion, data_insertion, data_ref_skip, data_quality = self.prepare_data_errorrates(self.data, mod, canonical)

        box_mismatch = go.Box(x=data_mismatch[0], y=data_mismatch[1], name="Mismatch rate", offsetgroup=0, line=dict(color="black"), marker=dict(outliercolor="black", size=2), fillcolor="#8c564b")
        box_deletion = go.Box(x=data_deletion[0], y=data_deletion[1], name="Deletion rate", offsetgroup=1, line=dict(color="black"), marker=dict(outliercolor="black", size=2), fillcolor="#e377c2")
        box_insertion = go.Box(x=data_insertion[0], y=data_insertion[1], name="Insertion rate", offsetgroup=2, line=dict(color="black"), marker=dict(outliercolor="black", size=2), fillcolor="#7f7f7f")
        box_ref_skip = go.Box(x=data_ref_skip[0], y=data_ref_skip[1], name="Reference skip", offsetgroup=3, line=dict(color="black"), marker=dict(outliercolor="black", size=2), fillcolor="#bcbd22")
        
        box_quality = go.Box(x=data_quality[0], y=data_quality[1], name="Mean quality", offsetgroup=0, line=dict(color="black"), marker=dict(outliercolor="black", size=2), fillcolor="#17becf")

        fig = self.update_plot(make_subplots(rows=1, cols=2, column_widths=[0.75, 0.25]), height=800, width=1200)
        fig.add_traces([box_mismatch, box_deletion, box_insertion, box_ref_skip], rows=[1,1,1,1], cols=[1,1,1,1])
        fig.add_trace(box_quality, row=1, col=2)
        fig.update_layout(boxmode="group")
        fig.update_yaxes(title_text="Error rate", row = 1, col = 1)
        fig.update_yaxes(title_text="Quality score", row = 1, col = 2)

        return to_html(fig, include_plotlyjs=False)

    def create_nb_plot(self, mod: str) -> str:
        """
        Create a neighbor position plot for a specified modification type.

        Args:
            mod (str): The modification type for which the neighbor position plot should be created.

        Returns:
            str: HTML code representing the neighbor position plot.

        Note:
            This method creates a neighbor position plot for a specified modification type. It prepares data for
            the plot, generates stacked bar charts for neighbor positions with and without modification errors,
            as well as a pie chart showing the distribution of surrounding errors, using Plotly, and returns the
            resulting plot as HTML code.

        """
        x_vals, y_vals, x_vals_mod, y_vals_mod, pie_labs, pie_vals = self.prepare_nb_counts(self.data, mod)

        fig = self.update_plot(make_subplots(rows=1, cols=2, column_widths=[0.75, 0.25], specs=[[{"type": "bar"}, {"type": "pie"}]]), width=1200)

        fig.add_trace(go.Bar(x=x_vals_mod, y=y_vals_mod, marker=dict(color="#d62728", line=dict(color='#000000', width=1.5)), name="nb. mod error"), row=1, col=1)
        fig.add_trace(go.Bar(x=x_vals, y=y_vals, marker=dict(color="#dd8452", line=dict(color='#000000', width=1.5)), name="nb. error"), row=1, col=1)
        fig.update_layout(barmode="stack")

        fig.update_xaxes(title = f"Relative position to {mod} positions", row=1, col=1)
        fig.update_yaxes(title = "Count", row=1, col=1)

        fig.add_trace(go.Pie(labels=pie_labs, values=pie_vals, name="", 
                            marker=dict(colors=["#4c72b0", "#dd8452", "#d62728"], line=dict(color='#000000', width=1.5))), row=1, col=2)

        fig.update_layout(showlegend=False)
        
        return to_html(fig, include_plotlyjs=False)
    
    
    ##############################################################################################################
    #                                               Output methods                                               #
    ##############################################################################################################

    def write_template(self, plots: List[str], category: str, corresponding_base: str) -> None:
        """
        Generate an HTML report template with interactive plots.

        Args:
            plots (List[str]): A list of HTML code representing interactive plots.
            category (str): The category or modification type for which the report is generated.
            corresponding_base (str): The corresponding base associated with the category.

        Returns:
            None: The template is written to an HTML file.
        """
        name = f"<i>{category}</i>"
        time = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')

        basename = os.path.splitext(os.path.basename(self.in_path))[0]
        out_path = f"{self.output_path}{category}_{basename}_summary.html"

        template = f"""
            <!DOCTYPE html>
            <html lang="en">
            <head>
                <meta charset="UTF-8">
                <meta http-equiv="X-UA-Compatible" content="IE=edge">
                <meta name="viewport" content="width=device-width, initial-scale=1.0">
                <title>Neet summary</title>
                <style>
                    /* Reset some default browser styles */
                    body, h1, h2, h3, h4, p, ul, li {{
                        margin: 0;
                        padding: 0;
                    }}

                    /* Apply modern font and line height */
                    body {{
                        font-family: 'Arial', sans-serif;
                        line-height: 1.6;
                    }}

                    /* Style the header */
                    header {{
                        background-color: #333;
                        color: white;
                        padding: 1em 0;
                        text-align: center;
                        width: 100%; /* Make the header span the full width */
                        box-shadow: 0 2px 4px rgba(0, 0, 0, 0.1);
                    }}

                    header h1 {{
                        font-size: 2.5em;
                        margin-bottom: 0.3em;
                    }}

                    header p {{
                        font-size: 1.2em;
                        opacity: 0.8;
                    }}


                    footer {{
                        background-color: #333;
                        color: white;
                        padding: 1em 0;
                        text-align: center;
                        width: 100%; /* Make the header span the full width */
                        box-shadow: 0 2px 4px rgba(0, 0, 0, 0.1);
                    }}


                    /* Center the content */
                    body {{
                        display: flex;
                        justify-content: center;
                        align-items: center;
                        min-height: 100vh;
                        flex-direction: column;
                        background-color: #f5f5f5;
                    }}

                    /* Style sections */
                    section {{
                        background-color: white;
                        border-radius: 15px;
                        border-style: solid;
                        border-color: #333;
                        box-shadow: 0 2px 4px rgba(0, 0, 0, 0.1);
                        margin: 1em 0;
                        max-width: 1420px;
                        width: 90%;
                        text-align: left;
                    }}

                    section h2 {{
                        background-color: #333;
                        border-radius: 10px 10px 0 0;
                        box-shadow: 0 2px 4px rgba(0, 0, 0, 0.1);
                        text-align: center;
                        color: white;
                        opacity: 0.9;       
                    }}

                    section h3 {{
                        background-color: #333;
                        border-radius: 5px;
                        box-shadow: 0 2px 4px rgba(0, 0, 0, 0.1);
                        text-align: center;
                        margin-top: 0.5em;
                        margin-left: 0.5em;
                        margin-right: 0.5em;
                        color: white;
                        opacity: 0.8;
                    }}

                    section h4 {{
                        background-color: #333;
                        border-radius: 5px;
                        box-shadow: 0 2px 4px rgba(0, 0, 0, 0.1);
                        text-align: center;
                        margin-top: 0.5em;
                        margin-left: 1.5em;
                        margin-right: 1.5em;
                        color: white;
                        opacity: 0.7;       
                    }}

                    /* Style lists */
                    ul {{
                        list-style-type: disc;
                        padding-left: 1em;
                        margin-left: 2.5em;
                        margin-top: 1em;
                        margin-bottom: 1em;
                    }}

                    li {{
                        margin-bottom: 0.25em;
                    }}

                    /* Style links */
                    a {{
                        color: #007bff;
                        text-decoration: none;
                    }}

                    a:hover {{
                        text-decoration: underline;
                    }}

                    p {{
                        font-size: 1.1rem;
                        margin-top: 0.5em;
                        margin-bottom: 2em;
                        margin-left: 2em;
                        margin-right: 1em;

                    }}

                    .plot-container {{
                        display: flex;
                        justify-content: center;
                        margin-top: 2em;
                        margin-left: 2em;
                        margin-right: 2em;
                        margin-bottom: 0.5em;
                        padding: 1em;
                    }}

                    .intro-text {{
                        font-size: 1.4rem;
                        margin-top: 1.5em;
                        margin-bottom: 1.5em;
                        margin-left: 1.5em;
                        margin-right: 1.5em;
                    }}

                </style>
            </head>

            <body>
                <script src="https://cdn.plot.ly/plotly-latest.min.js"></script>

                <header>
                    <h1>Positions of interest: {name}</h1>
                    <p>Produced by <a href="https://github.com/dietvin/neet">Neet</a> on <b>{time}</b></p>
                </header>
                
                <section>
                    <p class="intro-text">
                        This summary file was created from the extracted features in file <b>{self.in_path}</b> 
                        with <b>{category}</b> positions extracted from file <b>{self.bed_path}</b>. 
                        The plots are interactive and allow further information by hovering, zooming and panning.
                    </p>
                </section>

                <section>
                    <h2>{name} positions on reference sequence(s)</h2>
                    <h3>Mapping of {name} positions across the reference sequences</h3>
                    <div class="plot-container">
                        {plots[0]}
                    </div>
                    <p>
                        Positions of {name} are indicated as blue lines. Fasta file '{self.ref_path}' was used to extract
                        reference sequence(s). Hover on positions for exact coordinates. 
                    </p>
                </section>

                <section>
                    <h2>Mismatch types</h2>
                    <h3>Confusion matrix of {name} positions containing mismatch types</h3>
                    <div class="plot-container">
                        {plots[1]}
                    </div>
                    <p>
                        Overview of all types of (mis-)matches in the data subset corresponding. Different kind of matches
                        indicate possible mislabellings in the bed file, as {category} should only be found at {corresponding_base}
                        positions.
                    </p>
                </section>

                <section>
                    <h2>Base compositions</h2>
                    <h3>Base compositions for different {name} and {corresponding_base} subsets</h3>
                    <div class="plot-container">
                        {plots[2]}
                    </div>
                    <p>
                        Each A/C/G/T element in the bars corresponds to the median count of a given base
                        in a subset. The four medians were scaled to add up to one. {name} match: positions
                        labelled {name} from the bed file, where the called base is equal to the reference
                        base. {corresponding_base} match: positions with reference base {corresponding_base},
                        where the called base is equal. {name} mismatch: positions labelled {name} from the 
                        bed file, where the called base differs from the reference base. {corresponding_base} 
                        mismatch: positions with reference base {corresponding_base}, where the called base 
                        differs.
                    </p>
                </section>

                <section>
                    <h2>Error rates </h2>
                    <h3>Error rates for different {name} and {corresponding_base} subsets</h3>
                    <div class="plot-container">
                        {plots[3]}
                    </div>
                    <p>
                        Left: Distributions of mismatch, deletion, insertion and reference skip rates for different
                        subsets. Right: Distribution of mean quality scores for different subsets. {name} match: positions
                        labelled {name} from the bed file, where the called base is equal to the reference
                        base. {corresponding_base} match: positions with reference base {corresponding_base},
                        where the called base is equal. {name} mismatch: positions labelled {name} from the 
                        bed file, where the called base differs from the reference base. {corresponding_base} 
                        mismatch: positions with reference base {corresponding_base}, where the called base 
                        differs.
                    </p>
                </section>

                <section>
                    <h2>Neighbouring errors</h2>
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
                </section>
            </body>
            <footer></footer>
            </html> 
        """
        with open(out_path, "w") as out:
            out.write(template)

    def write_tsv(self) -> None:
        """
        Write the current DataFrame to a TSV (tab-separated values) file with additional bed information.

        The output filename is derived from the input filename by appending "_w_bed_info.tsv".

        Returns:
            None: The DataFrame is saved as a TSV file.
        """
        self.data.to_csv(f"{os.path.splitext(self.in_path)[0]}_w_bed_info.tsv", sep="\t", index=False, header=True)

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
        for category, corr_base in zip(self.bed_categories, self.bed_categories_canoncial_bases):
            self.process_category(category, corr_base)

    def process_category(self, category: str, corresponding_base: str):
        """
        Process and analyze the data for a specific category and corresponding base.

        This method generates multiple plots and saves them along with a summary HTML template. Optionally, it also
        writes the processed data to a TSV file.

        Parameters:
            category (str): The category of interest.
            corresponding_base (str): The corresponding base for the category.

        Returns:
            None: The plots and summary are saved to files, and optionally, the data is saved in TSV format.
        """
        plot_mod_map = self.create_map_plot(category)
        plot_mism_types = self.create_mism_types_plot(category)
        plot_comp = self.create_composition_plot(category, corresponding_base)
        plot_err_rate = self.create_error_rate_plot(category, corresponding_base)
        plot_nb = self.create_nb_plot(category)

        plots = [plot_mod_map, plot_mism_types, plot_comp, plot_err_rate, plot_nb]
        self.write_template(plots, category, corresponding_base)
        if self.output_tsv:
            self.write_tsv()

def setup_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="Neet - Position-of-Interest Analyzer", description="Analyze features of one or more types of positions of interest.")
    parser.add_argument('-i', '--tsv', type=str, required=True,
                        help="""
                            Path to the input file. Must be of type tsv, as returned by the PileupExtractor.
                            """)
    parser.add_argument('-o', '--output', type=str, required=True,
                        help="""
                            Path to output a output directory, in which all output files will be stored.
                            """)
    parser.add_argument('-b', '--bed', type=str, required=True,
                        help="""
                            Path to the bed file containing information in the fourth column.
                            """)
    parser.add_argument('-r', '--ref', type=str, required=True,
                        help="""
                            Path to the reference file. Must be of type fasta.
                            """)
    parser.add_argument('-c', '--bed_categories', type=str, required=True, 
                        help="""
                            One or more categories from the bed file to aggregate the data by. 
                            Must be in the format: 'cat1' or 'cat1,cat2,cat3'
                            """)
    parser.add_argument('-cc', '--counterparts', type=str, required=True, 
                        help="""
                            Canonical base corresponding to each category specified in --bed_categories. 
                            Same format as --bed_categories flag.
                            """)
    parser.add_argument('--update_tsv', action="store_true", 
                        help="""
                            If specified, the a new tsv file gets written containing the information from the bed 
                            file in the last column. Suffix '_w_bed_info' will be added to newly created file.
                            """)
    parser.add_argument('--use_perc_mismatch_alt', action="store_true", 
                        help="""
                            If specified, uses the mismatch from the perc_mismatch_alt colum.
                            """)
    return parser
    
if __name__ == "__main__":

    parser = setup_parser()
    args = parser.parse_args()

    poi_analyzer = POIAnalyzer(in_path=args.tsv,
                               out_path=args.output,
                               bed_path=args.bed,
                               ref_path=args.ref,
                               categories=args.bed_categories,
                               canonical_counterpart=args.counterparts,
                               output_tsv=args.update_tsv,
                               use_perc_mismatch_alt=args.use_perc_mismatch_alt)
    poi_analyzer.main()
                
    # poi_analyzer = POIAnalyzer("/home/vincent/masterthesis/data/45s_rrna/processed/dRNA_cytoplasm/dRNA_cytoplasm_extracted.tsv",
    #                            "/home/vincent/masterthesis/data/45s_rrna/processed/dRNA_cytoplasm",
    #                            "/home/vincent/masterthesis/data/45s_rrna/rRNA_modifications_conv_cleaned.bed",
    #                            "/home/vincent/masterthesis/data/45s_rrna/RNA45SN1_cleaned.fasta",
    #                            "psu,Gm,Um", "T,G,T", True, False)
    # poi_analyzer.main()
