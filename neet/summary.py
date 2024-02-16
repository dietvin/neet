import plotly.graph_objects as go
import plotly.express as px
from plotly.subplots import make_subplots
from plotly.io import to_html
from . import helper_functions as hs
import os, warnings, sys, datetime, argparse
from statistics import mean
from collections import defaultdict
from typing import List, Tuple, Dict

class SummaryCreator:
    input_path: str
    output_path: str
    n_bins: int | None
    perc_mis_col: str
    data: Dict[str, List[str|int|float]]
    export_svg: bool

    def __init__(self, in_path: str, out_path: str, n_bins: int|None, use_perc_mismatch_alt: bool, export_svg: bool) -> None:
        self.process_paths(in_path, out_path)
        self.n_bins = n_bins if n_bins != -1 else None
        self.perc_mis_col = "mismatch_rate_alt" if use_perc_mismatch_alt else "mismatch_rate"
        self.export_svg = export_svg
        
    #################################################################################################################
    #                                   Functions called during initialization                                      #
    #################################################################################################################

    def process_paths(self, in_path: str, out_path: str) -> None:
        """
        Process the input and output paths for the SummaryCreator instance.

        Parameters:
        - in_path (str): The input file path.
        - out_path (str): The output file path or directory.

        Raises:
        - FileNotFoundError: If the specified input file path does not exist or if the specified output directory does not exist.
        - Warning: If the output file has an extension other than '.html'. A warning is issued, but the function continues execution.

        Returns:
        - None
        """
        # process input path
        hs.check_input_path(in_path, [".tsv"])
        self.input_path = in_path
        # process output path(s)
        self.output_path = hs.process_outpath(out_path, f"{os.path.splitext(os.path.basename(in_path))[0]}_summary.html", [".html"])
      
    def load_data(self):
        """
        Load data from the specified input file into the SummaryCreator instance.

        Reads the data from a tab-separated values (tsv) file as created by the PileupExtractor module
        and stores it in the 'data' attribute of the class.

        Returns:
        - None
        """
        cols = ["chr", "n_reads", "ref_base", "majority_base", "deletion_rate", "insertion_rate", "refskip_rate", "mismatch_rate", "q_mean", "motif"]

        col_idx = {'chr': 0, 'site': 1, 'n_reads': 2, 'ref_base': 3, 'majority_base': 4, 'n_a': 5, 'n_c': 6, 'n_g': 7, 'n_t': 8, 'n_del': 9, 'n_ins': 10, 'n_ref_skip': 11, 'a_rate': 12, 'c_rate': 13, 'g_rate': 14, 'u_rate': 15, 'deletion_rate': 16, 'insertion_rate': 17, 'refskip_rate': 18, 'mismatch_rate': 19, 'mismatch_rate_alt': 20, 'motif': 21, 'q_mean': 22, 'q_std': 23, 'neighbour_error_pos': 24}
        dtypes = {'chr': str, 'site': int, 'n_reads': int, 'ref_base': str, 'majority_base': str, 'n_a': int, 'n_c': int, 'n_g': int, 'n_t': int, 'n_del': int, 'n_ins': int, 'n_ref_skip': int, 'a_rate': float, 'c_rate': float, 'g_rate': float, 'u_rate': float, 'deletion_rate': float, 'insertion_rate': float, 'refskip_rate': float, 'mismatch_rate': float, 'mismatch_rate_alt': float, 'motif': str, 'q_mean': float, 'q_std': float, 'neighbour_error_pos': str}
        
        with open(self.input_path, "r") as file:
            next(file)
            data = dict(zip(cols, [[] for _ in cols]))
            for line in file:
                line = line.strip().split("\t")
                for col in cols: 
                    data[col].append(dtypes[col](line[col_idx[col]]))
            self.data = data

    ######################################################################################################################
    #                                               Main processing method                                               #
    ######################################################################################################################
    def main(self) -> None:
        """
        Orchestrates the creation of a summary report from the input data file.

        This function performs the following steps:
        1. Loads data from the input file.
        2. Creates various summary plots including general statistics, chromosome-wise information,
        general mismatch statistics, specific mismatch type summaries, and motif summaries.
        3. Writes an HTML summary file containing the generated plots.

        Note: The function includes print statements to provide updates on the progress of each step.

        Raises:
            Exception: If any error occurs during the execution, an exception is raised with an
                    accompanying error message.

        Returns:
            None
        """
        hs.print_update(f"Starting creation of summary from file '{self.input_path}'.")
        hs.print_update("  - loading data... ", line_break=False)
        self.load_data()
        n_positions = len(self.data["chr"])
        n_chromosomes = len(set(self.data["chr"]))

        hs.print_update(f"Done. Found {n_positions} sites along {n_chromosomes} sequences.", with_time=False)

        plots = []
        
        hs.print_update("  - creating general summary... ", line_break=False)
        plots.append(self.create_general_plot())
        hs.print_update(f"Done.", with_time=False)

        hs.print_update("  - creating chromosome-wise summary... ", line_break=False)
        plots.append(self.create_chr_plot())
        hs.print_update(f"Done.", with_time=False)

        hs.print_update("  - creating general mismatch summary... ", line_break=False)
        plots.append(self.create_mism_general_plot())
        hs.print_update(f"Done.", with_time=False)

        hs.print_update("  - creating specific mismatch type summary... ", line_break=False)
        plots += self.create_mism_types_plots()
        hs.print_update(f"Done.", with_time=False)

        hs.print_update("  - creating motif summary... ", line_break=False)
        plots.append(self.create_motif_plot())
        hs.print_update(f"Done.", with_time=False)

        hs.print_update(f"  - creating HTML summary file at {self.output_path}")
        self.write_to_html(n_positions, n_chromosomes, plots)
        hs.print_update("Finished.")
    

    ################################################################################################################
    #                                                Helper methods                                                #
    ################################################################################################################

    def update_plot(self, fig, title: str|None = None, xlab: str|None = None, ylab: str|None = None, height: int = 500, width: int = 800):
        """
        Update the layout of the given Plotly figure.

        Parameters:
        - fig (plotly.graph_objs.Figure): The Plotly figure to be updated.
        - title (str|None): Title for the plot.
        - xlab (str|None): Label for the x-axis.
        - ylab (str|None): Label for the y-axis.
        - height (int): Height of the figure.
        - width (int): Width of the figure.

        Returns:
        - plotly.graph_objs.Figure: The updated Plotly figure.
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

    def bin_data(self, data: List[int|str|float]) -> List[int|float]:
        """
        Bin the input data into segments and return the mean value of each segment.

        Parameters:
        - data (List[int|str|float]): The input data to be binned.

        Raises:
        - Exception: If the number of bins is larger than the length of the data.

        Returns:
        - List[int|float]: List of mean values for each bin.
        """
        total_length = len(data)
        num_segments = self.n_bins

        if total_length == num_segments:
            return data
        elif total_length < num_segments:
            return data
        segment_length = segment_length = total_length // num_segments
        remainder = total_length % num_segments
        
        start_index = 0
        end_index = 0
        segments = []
        for i in range(num_segments):
            end_index = start_index + segment_length + (1 if i < remainder else 0)
            segments.append(data[start_index:end_index])
            start_index = end_index

        return [mean(segment) for segment in segments]

    ################################################################################################################
    #                                           Data preparation methods                                           #
    ################################################################################################################

    def prepare_data_general(self) -> Tuple[List[int|float], List[float]]:
        """
        Prepare the data for plotting general information.

        Returns:
        - Tuple[List[int|float], List[float]]: Tuple containing binned or original 'n_reads' and 'q_mean' data.
        """
        if self.n_bins is not None:
            n_reads = self.bin_data(self.data["n_reads"])
            quality = self.bin_data(self.data["q_mean"])
            return n_reads, quality
        else:
            return self.data["n_reads"], self.data["q_mean"]

    def prepare_data_chr(self) :
        data_grouped = defaultdict(lambda: {"n_reads": [], "q_mean": []})

        for chr_val, n_reads, q_mean in zip(self.data['chr'],
                                            self.data['n_reads'],
                                            self.data['q_mean']):
            data_grouped[chr_val]['n_reads'].append(n_reads)
            data_grouped[chr_val]['q_mean'].append(q_mean)
        data_grouped = dict(data_grouped)

        n_sites_x = list(data_grouped.keys())
        n_sites_y = [len(chr_data["n_reads"]) for chr_data in data_grouped.values()]

        chroms = []
        n_reads = []
        q_mean = []

        for chrom, chrom_data in data_grouped.items():
            if (self.n_bins is not None) & (len(chrom_data["n_reads"]) > self.n_bins):
                n_reads += self.bin_data(chrom_data["n_reads"])
                q_mean += self.bin_data(chrom_data["q_mean"])
                chroms += [chrom] * self.n_bins
            else:
                n_reads += chrom_data["n_reads"]
                q_mean += chrom_data["q_mean"]
                chroms += [chrom] * len(chrom_data["n_reads"])

        return n_sites_x, n_sites_y, chroms, n_reads, q_mean

    def prepare_data_mism_general(self) -> Dict[str, Tuple[int|List[float], int|List[float], int|List[float]]]:
        """
        Prepare mismatch data for plotting general mismatch information.

        Returns:
        - Dict[str, Tuple[int|List[float], int|List[float]]]: Dictionary containing binned or original mismatch data.
        """
        data_mis = {"mismatch_rate": [], "deletion_rate": [], "insertion_rate": [], "refskip_rate": []}
        data_del = {"mismatch_rate": [], "deletion_rate": [], "insertion_rate": [], "refskip_rate": []}
        data_mat = {"mismatch_rate": [], "deletion_rate": [], "insertion_rate": [], "refskip_rate": []}

        for d in zip(self.data["ref_base"], 
                     self.data["majority_base"], 
                     self.data["mismatch_rate"], 
                     self.data["deletion_rate"], 
                     self.data["insertion_rate"],
                     self.data["refskip_rate"]):
            if d[0] == d[1]:
                data_mat["mismatch_rate"].append(d[2])
                data_mat["deletion_rate"].append(d[3])
                data_mat["insertion_rate"].append(d[4])
                data_mat["refskip_rate"].append(d[5])
            elif d[1] == "-":
                data_del["mismatch_rate"].append(d[2])
                data_del["deletion_rate"].append(d[3])
                data_del["insertion_rate"].append(d[4])
                data_del["refskip_rate"].append(d[5])
            else:
                data_mis["mismatch_rate"].append(d[2])
                data_mis["deletion_rate"].append(d[3])
                data_mis["insertion_rate"].append(d[4])
                data_mis["refskip_rate"].append(d[5])

        data_dict = {}
        data_dict["overall"] = (len(data_mat["mismatch_rate"]), len(data_mis["mismatch_rate"]), len(data_del["mismatch_rate"]))
        for mismatch_type in ["mismatch_rate", "deletion_rate", "insertion_rate", "refskip_rate"]:
            if self.n_bins is not None:
                if len(data_mat[mismatch_type]) > self.n_bins:
                    data_mat[mismatch_type] = self.bin_data(data_mat[mismatch_type])
                if len(data_mis[mismatch_type]) > self.n_bins:
                    data_mis[mismatch_type] = self.bin_data(data_mis[mismatch_type])
                if len(data_del[mismatch_type]) > self.n_bins:
                    data_del[mismatch_type] = self.bin_data(data_del[mismatch_type])
            data_dict[mismatch_type] = (data_mat[mismatch_type], data_mis[mismatch_type], data_del[mismatch_type])

        return data_dict
    
    def prepare_data_mism_matrix(self, mis_count: Dict[str, int]) -> Tuple[List[List[int]], List[List[str]], int]:
        """
        Prepare data for a confusion matrix based on mismatch counts.

        Parameters:
        - mis_count (Dict[str, int]): Dictionary containing mismatch counts for each base pair combination.

        Returns:
        - Tuple[List[List[int]], List[List[str]], int]: Tuple containing matrix data, matrix labels, and the maximum value in the matrix.
        """
        # prepare data for the confusion matrix
        matrix_data = [[None]*5 for _ in range(5)]
        matrix_labels = [[None]*5 for _ in range(5)]
        bases = ["A", "C", "G", "U", "-"]
        # fill count matrix
        for i in range(5): 
            for j in range(5):
                count = mis_count[f"{bases[i]} - {bases[j]}"] if i<4 else 0 # 0 in case of "- - X"
                matrix_data[i][j] = count
                if i != j:
                    matrix_labels[i][j] = count
        # fill the matrix containing corresponding labels
        vals_flat = [element for sublist in matrix_labels for element in sublist if (element is not None)]
        vals_sum = sum(vals_flat)
        for i in range(5): 
            for j in range(5):
                if matrix_labels[i][j]:
                    matrix_labels[i][j] = f"{round(matrix_labels[i][j] / vals_sum * 100, 2)}%"
                else: 
                    matrix_labels[i][j] = ""
        # get the max value
        val_max = max(vals_flat)

        return matrix_data, matrix_labels, val_max

    def prepare_data_mism_pie(self, mis_count: Dict[str, int]) -> Tuple[List[int], List[str]]:
        """
        Prepare data for a pie chart based on mismatch counts.

        Parameters:
        - mis_count (Dict[str, int]): Dictionary containing mismatch counts for each base pair combination.

        Returns:
        - Tuple[List[int], List[str]]: Tuple containing pie chart data and labels.
        """
        # prepare data for the pie chart
        pie_data = list(mis_count.values())
        pie_labels = list(mis_count.keys())
        pie_data = []
        pie_labels = []
        # remove values for matches 
        for label, data in mis_count.items():
            if label not in ["A - A", "C - C", "G - G", "U - U"]:
                pie_data.append(data)
                pie_labels.append(label)
        return pie_data, pie_labels
    
    def prepare_data_mism_box(self, mis_error_rates: Dict[str, List[str|float]]) -> Dict[str, List[str|float]]:
        """
        Prepare data for box plots based on mismatch error rates by mismatch types.

        Parameters:
        - mis_error_rates (Dict[str, List[str|float]]): Dictionary containing mismatch error rates.

        Returns:
        - Dict[str, List[str|float]]: Dictionary containing the prepared data for the box plot.
        """
        if self.n_bins is not None:
            # group the data by mismatch type and if more samples than n_bins are present in a group, 
            # subset the data to n_bins data points
            error_rates_by_type = defaultdict(lambda: {"mismatch_rate": [], 
                                                    "deletion_rate": [], 
                                                    "insertion_rate": [], 
                                                    "refskip_rate": []})
            mismatch_types = ["A - C","A - G","A - U",
                            "C - A","C - G","C - U",
                            "G - A","G - C","G - U",
                            "U - A","U - C","U - G",
                            "A - -", "C - -", "G - -",
                            "U - -"]
            for mis_type in mismatch_types:
                mask = [mis_type == m for m in mis_error_rates["mismatch_type"]]
                for err_type in ["mismatch_rate", "deletion_rate", "insertion_rate", "refskip_rate"]:
                    # use mask to filter data corresponding to err_type
                    err_type_data = [mis_error_rates[err_type][i] for i in range(len(mask)) if mask[i]]

                    if len(err_type_data) > self.n_bins:
                        error_rates_by_type[mis_type][err_type] = self.bin_data(err_type_data)
                    else:
                        error_rates_by_type[mis_type][err_type] = err_type_data
            error_rates_by_type = dict(error_rates_by_type)

            # transform the data into the initial format
            mis_error_rates = {"mismatch_type": [], "mismatch_rate": [], "deletion_rate": [], "insertion_rate": [], "refskip_rate": []}
            for err_type, err_data in error_rates_by_type.items():
                mis_error_rates["mismatch_type"] += [err_type] * len(err_data["mismatch_rate"])
                mis_error_rates["mismatch_rate"] += err_data["mismatch_rate"]
                mis_error_rates["deletion_rate"] += err_data["deletion_rate"]
                mis_error_rates["insertion_rate"] += err_data["insertion_rate"]
                mis_error_rates["refskip_rate"] += err_data["refskip_rate"]   
            return mis_error_rates
        return mis_error_rates

    def prepare_data_mism_types(self) -> Tuple[List[List[int]], List[List[str]], int, List[int], List[str], Dict[str, List[str|float]]]:
        """
        Prepare general data for mismatch type plots.

        Returns:
        - Tuple[List[List[int]], List[List[str]], int, List[int], List[str], Dict[str, List[str|float]]]: 
        Tuple containing data for confusion matrix, pie chart, and box plot for mismatch types analysis.
        """
        mis_types = [f"{f} - {t}" for f in ["A", "C", "G", "U"] for t in ["A", "C", "G", "U", "-"]]
        mis_count = dict(zip(mis_types, [0]*len(mis_types)))

        mis_error_rates = {"mismatch_type": [],"mismatch_rate": [], "deletion_rate": [], "insertion_rate": [], "refskip_rate": []}

        for ref, maj, *error_rates in zip(self.data["ref_base"], 
                                          self.data["majority_base"], 
                                          self.data["mismatch_rate"], 
                                          self.data["deletion_rate"], 
                                          self.data["insertion_rate"], 
                                          self.data["refskip_rate"]):
            if (ref != "N") & (maj != "N"):
                mis_count[f"{ref} - {maj}"] += 1
                # prepare error rate data for distribution plots by each error type
                if ref != maj:
                    mis_error_rates["mismatch_type"].append(f"{ref} - {maj}")
                    mis_error_rates["mismatch_rate"].append(error_rates[0])
                    mis_error_rates["deletion_rate"].append(error_rates[1])
                    mis_error_rates["insertion_rate"].append(error_rates[2])
                    mis_error_rates["refskip_rate"].append(error_rates[3])

        return *self.prepare_data_mism_matrix(mis_count), *self.prepare_data_mism_pie(mis_count), self.prepare_data_mism_box(mis_error_rates)

    def prepare_data_motifs(self) -> Dict[str, Dict[str, List[str|float]]]:
        """
        Prepare data for motif-wise error rate plotting.

        Returns:
        - Dict[str, Dict[str, List[str|float]]]: Dictionary containing error rates for each motif.
        """
        motif_error_rates = defaultdict(lambda: {"mismatch_rate": [], "deletion_rate": [], "insertion_rate": [], "refskip_rate": []})

        # extract the data for each 3 base-pair motif; store error rates in dict by motifs
        motif_center_idx = len(self.data["motif"][0]) // 2
        for motif, *error_rates in zip(self.data["motif"], self.data["mismatch_rate"], self.data["deletion_rate"], self.data["insertion_rate"], self.data["refskip_rate"]):
            motif_3bp = motif[motif_center_idx-1:motif_center_idx+2]
            for i, err_type in enumerate(["mismatch_rate", "deletion_rate", "insertion_rate", "refskip_rate"]):
                motif_error_rates[motif_3bp][err_type].append(error_rates[i])
        motif_error_rates = dict(motif_error_rates)

        # subsample data for each error rate for each motif
        if self.n_bins is not None:
            for motif in motif_error_rates.keys():
                if len(motif_error_rates[motif]["mismatch_rate"]) > self.n_bins:
                    for err_type in ["mismatch_rate", "deletion_rate", "insertion_rate", "refskip_rate"]:
                        motif_error_rates[motif][err_type] = self.bin_data(motif_error_rates[motif][err_type])

        # sort dict by motif
        motif_error_rates = dict(sorted(motif_error_rates.items()))

        # set up data for each subplot (i.e. center A, C, G, U)
        center_base_error_rates = defaultdict(lambda: {"motif": [], "mismatch_rate": [], "deletion_rate": [], "insertion_rate": [], "refskip_rate": []})
        for motif, error_rates in motif_error_rates.items():
            center_base = motif[1]
            if center_base == "N": continue # only subplots for center A, C, G, U
            num_sites = len(error_rates["mismatch_rate"])
            center_base_error_rates[center_base]["motif"] += [motif] * num_sites
            for i, err_type in enumerate(["mismatch_rate", "deletion_rate", "insertion_rate", "refskip_rate"]):
                center_base_error_rates[center_base][err_type] += error_rates[err_type]

        return dict(center_base_error_rates)
    
    ##############################################################################################################
    #                                             Plotting functions                                             #
    ##############################################################################################################

    def create_error_placeholder(self, e: Exception):
        """
        Create a placeholder plot with an error message.

        Parameters:
        - e (Exception): The exception that occurred.

        Returns:
        - fig: Plotly Figure: Placeholder figure with an error message.
        """
        hs.print_update(f"An error occured: {str(e)}. Replacing plot with empty placeholder. ", with_time=False, line_break=False)
        fig = self.update_plot(make_subplots(rows=1, cols=1))
        fig.add_trace(go.Scatter(x=[0], y=[0], mode='text', text=[f"An error occured: {str(e)}"]))
        fig.update_layout(
            dragmode=False,  # Disable panning
            hovermode='closest',  # Maintain hover behavior
            uirevision='true'  # Disable double-click autoscale
        )
        return fig


    def create_general_plot(self) -> go.Figure:
        """
        Create a general plot displaying total coverage and total quality distributions.

        Returns:
        - go.Figure: Plotly Figure for general plot.
        """
        try:
            n_reads_data, quality_data = self.prepare_data_general()

            fig = make_subplots(rows=1, cols=2, 
                                subplot_titles=["Total coverage distribution", "Total quality distribtion"], 
                                horizontal_spacing=0.1)
            fig = self.update_plot(fig, width=1400)
            fig.update_annotations(font_size=25) # subplot title sizes

            fig.add_traces([go.Box(y=n_reads_data, name="Total"), go.Box(y=quality_data, name="Total")], rows=[1,1], cols=[1,2])

            fig.update_traces(line=dict(color="black"), 
                            marker=dict(outliercolor="black"), 
                            fillcolor="lightgrey")
            fig.update_layout(showlegend=False)
            fig.update_xaxes(showticklabels=False, ticks=None)
            fig.update_yaxes(title_text="Coverage", row=1, col=1)
            fig.update_yaxes(title_text="Mean quality", row=1, col=2)

        except Exception as e:
            fig = self.create_error_placeholder(e)

        return fig
    
    def create_chr_plot(self) -> go.Figure:
        """
        Create a chromosome plot displaying the number of positions, number of reads, and mean quality.

        Returns:
        - go.Figure: Plotly Figure for chromosome plot.
        """
        try:
            def custom_sort_key(item):
                if item.isdigit():  # Check if the item is a digit
                    return (int(item),)  # Convert to integer and sort numerically
                else:
                    return (float('inf'), item)  # Place non-digits at the end

            n_pos_x, n_pos_y, chrom, n_reads, quality = self.prepare_data_chr()

            fig = make_subplots(rows=3, cols=1, 
                                specs=[[{"type": "bar"}], [{"type": "box"}], [{"type": "box"}]], 
                                shared_xaxes=True, vertical_spacing=0.02)
            fig = self.update_plot(fig, height=1000, width=1400)

            fig.add_trace(go.Bar(x=n_pos_x, y=n_pos_y))
            fig.add_trace(go.Box(x=chrom, y=n_reads), row=2, col=1)
            fig.add_trace(go.Box(x=chrom, y=quality), row=3, col=1)

            fig.update_traces(marker=dict(color='lightgrey', line=dict(color='black', width=2)), selector=dict(type='bar'))
            fig.update_traces(marker=dict(color="black", outliercolor="black", size=2), fillcolor="lightgrey", selector=dict(type='box'))
            
            fig.update_xaxes(categoryorder='array', categoryarray=sorted(n_pos_x, key=custom_sort_key))
            for i in range(1, 3):
                fig.update_xaxes(tickvals=[], ticktext=[], row=i, col=1)

            fig.update_xaxes(title_text="Chromosome", row=3, col=1)
            fig.update_yaxes(title_text="Number of positions", row=1, col=1)
            fig.update_yaxes(title_text="Number of reads", row=2, col=1)
            fig.update_yaxes(title_text="Mean quality", row=3, col=1)
            fig.update_layout(showlegend=False)

        except Exception as e:
            fig = self.create_error_placeholder(e)

        return fig

    def create_mism_general_plot(self) -> go.Figure:
        """
        Create a general plot displaying mismatch frequencies, deletion frequencies, 
        insertion frequencies and reference skip frequencies at matched and mismatched sites, 
        and an overall pie chart.

        Returns:
        - go.Figure: Plotly Figure for general mismatch plot.
        """
        try:
            def boxplot(data, group: str = "match", showlegend: bool = False):
                if group == "match":
                    name = "Match"
                    fillcolor = "#4c72b0"
                    legendgroup = "match"
                elif group == "mismatch":
                    name = "Mismatch"
                    fillcolor = "#dd8452"
                    legendgroup = "mism"
                else:
                    name = "Deletion"
                    fillcolor = "#2ca02c"
                    legendgroup = "del"
                return go.Box(y=data, name=name, marker=dict(color="black", outliercolor="black", size=2), 
                            fillcolor=fillcolor, showlegend=showlegend, legendgroup=legendgroup)

            data_processed = self.prepare_data_mism_general()

            fig = make_subplots(rows=1, cols=5, 
                                specs=[[{"type": "box"}, {"type": "box"}, {"type": "box"}, {"type": "box"}, {"type": "pie"}]], 
                                shared_yaxes=True,
                                horizontal_spacing=0.01,
                                column_titles=("Mismatch frequency", "Deletion frequency", "Insertion frequency", "Reference skip frequ.", None))
            fig = self.update_plot(fig, height=600, width=1400)

            for i, (dname, legend) in enumerate(zip(["mismatch_rate", "deletion_rate", "insertion_rate", "refskip_rate"], [True, False, False, False]), start=1):
                box_mat = boxplot(data_processed[dname][0], group="match", showlegend=legend)
                box_mis = boxplot(data_processed[dname][1], group="mismatch", showlegend=legend)
                box_del = boxplot(data_processed[dname][2], group="deletion", showlegend=legend)
                fig.add_traces([box_mat, box_mis, box_del], rows=1, cols=i)

            pie = go.Pie(labels=["Match", "Mismatch", "Deletion"], 
                        values=[data_processed["overall"][0], data_processed["overall"][1], data_processed["overall"][2]], 
                        hoverinfo="label+percent", 
                        textinfo="value", 
                        textfont_size=20, 
                        marker=dict(colors=["#4c72b0", "#dd8452", "#2ca02c"], line=dict(color="#000000", width=2)), 
                        showlegend=False,
                        sort=False)
            fig.add_trace(pie, row=1, col=5)

            fig.update_layout(legend=dict(orientation="h", yanchor="bottom", y=1.05, xanchor="right", x=1, bgcolor='#f5f5f5', bordercolor='#000000', borderwidth=2))
            fig.update_annotations(font_size=21)
            fig.update_yaxes(showticklabels=False, ticks=None, row=1, col=2)
            fig.update_yaxes(showticklabels=False, ticks=None, row=1, col=3)
            fig.update_yaxes(showticklabels=False, ticks=None, row=1, col=4)

        except Exception as e:
            fig = self.create_error_placeholder(e)

        return fig

    def create_mism_types_plots(self) -> List[go.Figure]:
        """
        Create mismatch type plots including a confusion matrix, a pie chart, and a box plot.

        Returns:
        - List[go.Figure]: List of Plotly Figures for mismatch type plots.
        """
        try:
            matrix_data, matrix_labels, matrix_max_mism, pie_data, pie_labels, box_data = self.prepare_data_mism_types()
        except Exception as e:
            fig = self.create_error_placeholder(e)
            fig = fig
            return [fig, fig, fig]

        try:
            fig = px.imshow(matrix_data, labels=dict(x="Called base", y="Reference base", color="Count"), zmin=0, zmax=1.2*matrix_max_mism, color_continuous_scale="portland")
            fig = self.update_plot(fig, None, "Called base", "Reference base", width=800)
            fig.update_traces(text=matrix_labels, texttemplate="%{text}")
            fig.update_xaxes(fixedrange=True, tickvals=[0,1,2,3,4], ticktext=["A", "C", "G", "U", "-"])
            fig.update_yaxes(fixedrange=True, tickvals=[0,1,2,3,4], ticktext=["A", "C", "G", "U", "-"])
            matrix = fig
        except Exception as e:
            fig = self.create_error_placeholder(e)
            matrix = fig

        try:
            fig = go.Figure(go.Pie(labels=pie_labels, values=pie_data, 
                        hoverinfo='label+percent', textinfo='value', textfont_size=20, 
                        marker=dict(line=dict(color='#000000', width=2))))
            fig = self.update_plot(fig, width=800)
            fig.update_layout(legend=dict(bgcolor='#f5f5f5', bordercolor='#000000', borderwidth=2))
            pie = fig
        except Exception as e:
            fig = self.create_error_placeholder(e)
            pie = fig

        try:
            fig = go.Figure()
            fig = self.update_plot(fig, ylab="Error rate", height=600, width=1400)

            x = box_data["mismatch_type"]
            for err_type, c in zip(["mismatch_rate", "deletion_rate", "insertion_rate", "refskip_rate"], ["#55a868", "#c44e52", "#8172b3", "#937860"]):
                y = box_data[err_type]
                fig.add_trace(go.Box(x=x, y=y, name=err_type, 
                                    line=dict(color="black"), 
                                    marker=dict(outliercolor="black", size=2), 
                                    fillcolor=c))

            fig.update_layout(boxmode="group", legend=dict(orientation="h", yanchor="bottom", y=1.05, xanchor="right", x=1.0, bgcolor='#f5f5f5', bordercolor='#000000', borderwidth=2))
            fig.update_xaxes(categoryorder='array', categoryarray=pie_labels)
            box = fig
        except Exception as e:
            fig = self.create_error_placeholder(e)
            box = fig

        return [matrix, pie, box]

    def create_motif_plot(self) -> go.Figure:
        """
        Create a plot showing error rates for different motifs.

        Returns:
        - go.Figure: Plotly Figure for the motif plot.
        """
        try:
            x_order = {"A": ["CAC", "CAG", "CAU", "GAC", "GAG", "GAU", "UAC", "UAG", "UAU"], 
                    "C": ["ACA", "ACG", "ACU", "GCA", "GCG", "GCU", "UCA", "UCG", "UCU"],
                    "G": ["AGA", "AGC", "AGU", "CGA", "CGC", "CGU", "UGA", "UGC", "UGU"],
                    "U": ["AUA", "AUC", "AUG", "CUA", "CUC", "CUG", "GUA", "GUC", "GUG"]}
            
            fig = make_subplots(rows=2, cols=2, 
                                specs=[[{"type": "box"}, {"type": "box"}], [{"type": "box"}, {"type": "box"}]], 
                                shared_yaxes=True, 
                                vertical_spacing=0.1, horizontal_spacing=0.05)
            fig = self.update_plot(fig, height=900, width=1400)

            motif_data = self.prepare_data_motifs()

            for center_base, row, col in zip(["A", "C", "G", "U"], [1,1,2,2], [1,2,1,2]):
                d = motif_data[center_base]
                traces = []
                for i, (err_type, name, color) in enumerate(zip(["mismatch_rate", "deletion_rate", "insertion_rate", "refskip_rate"], ["Mismatch", "Deletion", "Insertion", "Reference skip"], ["#55a868", "#c44e52", "#8172b3", "#937860"])):
                    trace = go.Box(x=d["motif"], y=d[err_type], name=name, 
                                line=dict(color="black", width=1), 
                                marker=dict(outliercolor="black", size=1), 
                                fillcolor=color,
                                legendgroup=i,
                                showlegend=True if center_base == "A" else False,
                                offsetgroup=i)
                    traces.append(trace)

                fig.add_traces(traces, rows=row, cols=col)
                fig.update_xaxes(categoryorder='array', categoryarray=x_order[center_base], row=row, col=col)

            fig.update_layout(boxmode="group", boxgroupgap=0, legend=dict(orientation="h", yanchor="bottom", y=1.02, xanchor="right", x=1, bgcolor='#f5f5f5', bordercolor='#000000', borderwidth=2))
            fig.update_xaxes(title = "3bp Motif", row=2, col=1)
            fig.update_xaxes(title = "3bp Motif", row=2, col=2)
            fig.update_yaxes(title = "Error rate", row=1, col=1)
            fig.update_yaxes(title = "Error rate", row=2, col=1)

        except Exception as e:
            fig = self.create_error_placeholder(e)            
        return fig

    ################################################################################################################
    #                                               Create HTML file                                               #
    ################################################################################################################
    def write_svg(self, fig: go.Figure, name: str) -> None:
        """
        Write the Plotly Figure to an SVG file.

        Args:
        - fig (go.Figure): Plotly Figure.
        - name (str): Name to include in the output SVG file.

        Returns:
        - None
        """
        outpath = f"{os.path.splitext(self.output_path)[0]}_{name}.svg"
        fig.write_image(outpath)

    def figs_to_str(self, plot_figs: List[go.Figure]) -> List[str]:
        """
        Convert a list of Plotly Figures to their HTML string representation.

        Args:
        - plot_figs (List[go.Figure]): List of Plotly Figures.

        Returns:
        - List[str]: List of HTML strings representing the Plotly Figures.
        """
        plot_str = list(map(lambda x: to_html(x, include_plotlyjs=False), plot_figs))
        return plot_str

    def write_to_html(self, n_positions, n_chr, plot_figs: List[go.Figure]) -> None:
        """
        Generate an HTML summary page with collapsible sections for different types of analysis.
        
        Parameters:
        - n_positions (int): Total number of positions extracted.
        - n_chr (int): Total number of chromosomes.
        - plot_figs (List[go.Figure]): List of Plotly figures for different analyses.

        The function exports SVG versions of the Plotly figures if `export_svg` is True.

        The HTML template includes collapsible sections for general statistics, mismatch statistics,
        error rates by motifs, and more. Each section contains relevant Plotly charts and informative text.

        The final HTML template is written to the specified `output_path`.

        Returns:
        None
        """
        if self.export_svg:
            for fig, name in zip(plot_figs, ["summary_general_info", "summary_chr_info", "summary_mismatch_stats", 
                                             "summary_mismatch_matrix", "summary_mismatch_pie", 
                                             "summary_error_rates_by_mismatch", "summary_error_rates_by_motif"]):
                self.write_svg(fig, name)
        plots = self.figs_to_str(plot_figs)

        css_string, plotly_js_string = hs.load_html_template_str()

        time = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
        data_point_descr = f"Each data point corresponds to the average value along all positions in one of {self.n_bins} bins" if self.n_bins else "Each data point corresponds to one extracted position"
        template = f"""
                    <!DOCTYPE html>
                    <html lang="en">
                    <head>
                        <meta charset="UTF-8">
                        <meta http-equiv="X-UA-Compatible" content="IE=edge">
                        <meta name="viewport" content="width=device-width, initial-scale=1.0">
                        <title>Neet - Position extractor summary</title>
                        <style>{css_string}</style>
                    </head>

                    <body>
                        <script>{plotly_js_string}</script>

                        <header>
                            <h1>Pileup extractor summary</h1>
                            <p>Produced by <a href="https://github.com/dietvin/neet">Neet</a> on <b>{time}</b></p>
                        </header>
                    
                        <section>
                            <p class="intro-text">
                                This summary file was created from the extracted features in file <b>{self.input_path}</b>. 
                                {f"Data was averaged into <b>{self.n_bins}</b> bins to allow for better performance." if self.n_bins else ""}
                                In total <b>{n_positions}</b> positions were extracted along <b>{n_chr}</b> {"sequences" if n_chr > 1 else "sequence"}. 
                                The plots are interactive and allow further information by hovering, zooming and panning.
                            </p>
                        </section>

                        <section>
                            <button class="collapsible">General statistics</button>
                            <div class="collapsible-content">
                                <h2 class="hiddentitle" id="general_statistics"></h2>

                                <h3>Over all extracted position</h3>
                                <div class="plot-container">
                                    {plots[0]}
                                </div>
                                <p>
                                    Distribution of the coverage (<b>left</b>) and mean quality (<b>right</b>) of those positions.  
                                    The mean quality at a given position x is calculated from the quality scores from all mapped reads
                                    at this position. {data_point_descr}.
                                </p>

                                <h3>Split by each chromosome</h3>        
                                <div class="plot-container">
                                    {plots[1]}
                                </div>
                                <p>
                                    <b>Top</b>: Number of positions that were extracted on each chromosome. <b>Middle</b>: Distribution of the number of reads on a chromosome. 
                                    Each data point corresponds to one position. <b>Bottom</b>: Distribution of the quality scores averaged over all reads mapped
                                    at a given position. {data_point_descr}.
                                </p>
                            </div>
                        </section>

                        <section>
                            <button class="collapsible">Mismatch statistics</button>
                            <div class="collapsible-content">
                                <h2 class="hiddentitle" id="mismatch_statistics"></h2>

                                <h3>Mismatch, deletetion and insertion rates for matched and mismatched positions</h3>
                                <div class="plot-container">
                                    {plots[2]}
                                </div>
                                <p>
                                    Overview of the number of mismatches and what types of errors contribute to them. Match refers to the positions where the correct base was called. 
                                    Mismatch refers to the positions where the wrong base was called. The pie chart on the right shows the number of matched and mismatched positions
                                    along all chromosomes. The boxplots on the left show the distributions of the mismatch (<b>leftmost</b>), deletion (<b>second from left</b>), insertion (<b>second from right</b>)
                                    and reference skip (<b>right</b>) rates at matched and mismatched positions. {data_point_descr}.
                                </p>

                                <h3>Abundances of different type of mismatches</h3>

                                <h4>As confusion matrix</h4>
                                <div class="plot-container">
                                    {plots[3]}
                                </div>
                                <p>
                                    Abundance of matches by base (diagonal) and all types of mismatches <i>from Reference base to Called base</i>. Warmer colors indicate higher counts.
                                </p>
                            
                                <h4>As pie chart</h4>
                                <div class="plot-container">
                                    {plots[4]}
                                </div>
                                <p>
                                    Relative abundances of mismatch types (<i>[FROM] - [TO]</i>). Section labels show absolute count. Relative count is shown on hovering.
                                </p>

                                <h3>Mismatch, deletion and insertion rates by type of mismatch</h3>
                                <div class="plot-container">
                                    {plots[5]}
                                </div>
                                <p>
                                    Distribtions of Mismatch, Deletion, Insertion and Refernce skip rates for each observed mismatch type (<i>[FROM] - [TO]</i>). 
                                    Each data point corresponds to one position.
                                </p>
                            </div>
                        </section>

                        <section>
                            <button class="collapsible">Error rate by motifs</button>
                            <div class="collapsible-content">
                                <h2 class="hiddentitle" id="error_motif"></h2>

                                <h3>Mismatch, insertion, deletion and reference skip rates for 3bp motifs</h3>
                                <div class="plot-container">
                                    {plots[6]}
                                </div>
                                <p>
                                    Distributions of Mismatch, deletion and insertion rates for different three base motifs with center A (<b>top left</b>), C (<b>top right</b>),
                                    G (<b>bottom left</b>) and U (<b>bottom right</b>). {data_point_descr}.
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
        with open(self.output_path, "w") as o:
            o.write(template)