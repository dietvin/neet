import pandas as pd
import numpy as np
import plotly.graph_objects as go
import plotly.express as px
from plotly.subplots import make_subplots
from plotly.io import to_html
import helper_functions as hs

import os, warnings, sys, datetime, argparse
from typing import List, Tuple, Dict

class SummaryCreator:
    """
    A class for creating summary outputs based on input data.

    Attributes:
        input_path (str): Path to the input data file.
        output_path (str): Path to the output summary file.
        n_bins (int|None): Number of bins for summarization (optional).
        data (pd.DataFrame): Loaded input data as a pandas DataFrame.

    Methods:
        __init__(self, in_path: str, out_path: str, n_bins: int|None) -> None:
            Initializes the SummaryCreator object.
        process_path(self, in_path: str, out_path: str) -> None:
            Processes input and output paths.
        check_path(self, path: str, extensions: List[str]) -> None:
            Checks if the given path exists and has expected extensions.
        process_outpath(self, out: str) -> str:
            Processes the output path or filename.
        load_data(self) -> None:
            Loads data from the input path into the DataFrame.
    """
    input_path: str
    output_path: str
    n_bins: int | None
    perc_mis_col: str
    data: pd.DataFrame
    export_svg: bool

    dtypes = {'chr': str, 'site': int, 'n_reads': int, 'ref_base': str, 'majority_base': str, 'n_a': int, 'n_c': int,
            'n_g': int, 'n_t': int, 'n_del': int, 'n_ins': int, 'n_ref_skip': int, 'n_a_rel': float, 'n_c_rel': float,
            'n_g_rel': float, 'n_t_rel': float, 'n_del_rel': float, 'n_ins_rel': float, 'n_ref_skip_rel': float,
            'perc_mismatch': float, 'perc_mismatch_alt': float, 'motif': str, 'q_mean': float, 'q_std': float,
            'neighbour_error_pos': str}

    def __init__(self, in_path: str, out_path: str, n_bins: int|None = 5000, use_perc_mismatch_alt: bool = False, export_svg: bool = False) -> None:
        """
        Initializes a SummaryCreator object.

        Args:
            in_path (str): Path to the input data file.
            out_path (str): Path to the output summary file.
            n_bins (int|None): Number of bins for summarization (optional).
        """
        self.process_path(in_path, out_path)
        self.n_bins = n_bins
        self.perc_mis_col = "perc_mismatch_alt" if use_perc_mismatch_alt else "perc_mismatch"
        self.export_svg = export_svg
        
    #################################################################################################################
    #                                   Functions called during initialization                                      #
    #################################################################################################################

    def process_path(self, in_path: str, out_path: str) -> None:
        """
        Processes input and output path.

        Args:
            in_path (str): Path to the input data file.
            out_path (str): Path to the output summary file.
        """
        # process input path
        self.check_path(in_path, [".tsv"])
        self.input_path = in_path
        # process output path(s)
        self.output_path = self.process_outpath(out_path)

    def check_path(self, path: str, extensions: List[str]) -> None:
        """
        Checks if the given path exists and has expected extensions.

        Args:
            path (str): Path to check.
            extensions (List[str]): List of valid file extensions.

        Raises:
            FileNotFoundError: If the file does not exist.
            Warning: If the file extension is unexpected.
        """
        if not os.path.exists(path): # does file exist?
            raise FileNotFoundError(f"Input file not found. File '{path}' does not exist.")
        file_type = os.path.splitext(path)[1]
        if not file_type in extensions:
            warnings.warn(f"Found file extension {file_type}. Expected file extension to be one of: {extensions}. If this is deliberate, ignore warning.", Warning)

    def process_outpath(self, out: str) -> str:
        """
        Processes the output path or filename.

        Args:
            out (str): Output path or filename.

        Returns:
            str: Processed output path.

        Raises:
            FileNotFoundError: If the directory does not exist.
            Warning: If the output file extension is unexpected.
        """
        if os.path.isdir(out): # check if outpath is directory, if the directory exists and create output file path(s) according to input file name(s)
            if not os.path.exists(out):
                raise FileNotFoundError(f"Directory not found. Output directory '{out}' does not exist.")
            if not out.endswith("/"):
                out += "/"
            basename = os.path.splitext(os.path.basename(self.input_path))[0]
            out_path = f"{out}{basename}_summary.html"
            return out_path

        else: # check if outpath is a list of filename(s), if the basedirectory exists and if the given file extension(s) is .tsv
            dirname = os.path.dirname(out)
            if not os.path.exists(dirname):
                raise FileNotFoundError(f"Path to output file not found. '{dirname}' does not exist.")
            file_extension = os.path.splitext(out)[1]
            if file_extension != ".html":
                warnings.warn(f"Given output file has extension '{file_extension}'. Note that the output file will be of type '.tsv'.")
            return out
      
    def load_data(self) -> None:
        self.data = pd.read_csv(self.input_path, sep="\t", dtype=self.dtypes)

    ######################################################################################################################
    #                                               Main processing method                                               #
    ######################################################################################################################
    def create_summary(self) -> None:
        """
        Generate a summary report containing various plots and statistics based on the provided data.

        This method orchestrates the creation of multiple types of plots and summary visualizations, such
        as coverage distributions, chromosome-specific statistics, mismatch statistics, error type distributions,
        motif-related statistics, and more. The generated plots and statistics are combined into an HTML report.

        The final HTML report is saved to the specified output path.

        Returns:
            None: The method generates and writes an HTML report containing summary plots and statistics
                to the specified output path.
        """
        hs.print_update(f"Starting creation of summary from file '{self.input_path}'.")
        hs.print_update("1 - loading data")
        self.load_data()
        hs.print_update("2 - creating general summary")

        n_positions = self.data.shape[0]
        n_chromosomes = len(self.data["chr"].unique())
        plots = []
        
        plots.append(self.create_general_plot())
        hs.print_update("3 - creating chromosome-wise summary")
        plots.append(self.create_chr_plot())
        hs.print_update("4 - creating general mismatch summary")
        plots.append(self.create_mism_general_plot())
        hs.print_update("5 - creating specific mismatch type summary")
        plots += self.create_mism_types_plots()
        hs.print_update("6 - creating motif summary")
        plots.append(self.create_motif_plot())
        hs.print_update(f"7 - creating HTML summary file at {self.output_path}")
        self.write_to_html(n_positions, n_chromosomes, plots)
        hs.print_update("Finished.")
    

    ################################################################################################################
    #                                                Helper methods                                                #
    ################################################################################################################

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

    def bin_data(self, data: np.ndarray|pd.Series, num_segments: int = 5000) -> np.ndarray:
        """
        Bins input data into segments and calculates the mean value of each segment.

        Args:
            data (np.ndarray|pd.Series): Input data to be binned and summarized.
            num_segments (int): Number of segments to divide the data into (default: 5000).

        Returns:
            np.ndarray: Array of mean values for each segment.
        """
        # Calculate the length of each segment
        data_len = len(data)
        segment_length = data_len // num_segments
        remainder = data_len % num_segments
        
        # Split the array into segments with adjusted lengths
        segments = np.split(data, [segment_length * i + min(i, remainder) for i in range(1, num_segments)])
        # Calculate the mean value of each segment using numpy.mean()
        segment_means = np.array([segment.mean() for segment in segments])
        
        return segment_means

    ################################################################################################################
    #                                           Data preparation methods                                           #
    ################################################################################################################

    def prepare_data_general(self) -> Tuple[np.ndarray, np.ndarray]:
        """
        Prepares data for the general summary plot.

        Returns:
            Tuple[np.ndarray, np.ndarray]: Processed data arrays for total coverage and mean quality.
        """
        data = self.data
        n_reads = data["n_reads"].to_numpy()
        quality = data["q_mean"].to_numpy()
        if self.n_bins:
            n_reads = self.bin_data(n_reads, self.n_bins)
            quality = self.bin_data(quality, self.n_bins)
        return n_reads, quality

    def prepare_data_chr(self):
        """
        Prepares data for the chromosome-specific summary plot.

        Returns:
            Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]: Processed data arrays for the plot.
        """        
        data = self.data
        def bin_two_cols(d_group: Tuple[str,pd.DataFrame]):
            chr_val = d_group[0]
            d = d_group[1]
            del(d_group)
            return np.full(shape=self.n_bins, fill_value=chr_val), self.bin_data(d["n_reads"], self.n_bins), self.bin_data(d["q_mean"], self.n_bins)

        data = data[["chr", "n_reads", "q_mean"]]
        data_grouped = data.groupby("chr")

        n_pos = data_grouped.size()
        n_pos_x, n_pos_y = n_pos.index, n_pos.values
        del(n_pos)

        if self.n_bins:
            chrom_data = []
            n_reads_data = []
            quality_data = []
            
            for subset in data_grouped:
                # get the binned information for one chromosome
                tmp = bin_two_cols(subset)  
                # append it to the information from other chromosomes
                chrom_data.append(tmp[0])
                n_reads_data.append(tmp[1])
                quality_data.append(tmp[2])

            # concat information from all chromosomes to one array
            chrom_data = np.concatenate(chrom_data)
            n_reads_data = np.concatenate(n_reads_data)
            quality_data = np.concatenate(quality_data)
        else:
            chrom_data = data["chr"].to_numpy()
            n_reads_data = data["n_reads"].to_numpy()
            quality_data = data["q_mean"].to_numpy()

        return n_pos_x, n_pos_y, chrom_data, n_reads_data, quality_data

    def prepare_data_mism_general(self) -> Dict[str, Tuple[np.ndarray, np.ndarray]]:
        """
        Prepares data for the general mismatch summary plot.

        Returns:
            Dict[str, Tuple[np.ndarray, np.ndarray]]: Processed data for mismatch summary.
        """
        data = self.data[["ref_base", "majority_base", self.perc_mis_col, "n_del_rel", "n_ins_rel", "n_ref_skip_rel"]]
        data_mis = data.loc[data["ref_base"]!=data["majority_base"]]
        data_mat = data.loc[data["ref_base"]==data["majority_base"]]
        del(data)

        data_dict = {}
        data_dict["overall"] = (data_mat.shape[0], data_mis.shape[0])
        for i, j in zip(["mis_rate", "del_rate", "ins_rate", "ref_skip_rate"], [self.perc_mis_col, "n_del_rel", "n_ins_rel", "n_ref_skip_rel"]):
            data_dict[i] = (self.bin_data(data_mat[j], self.n_bins), self.bin_data(data_mis[j], self.n_bins)) if self.n_bins else (data_mat[j].to_numpy(), data_mis[j].to_numpy())

        return data_dict
    
    def prepare_data_mism_types(self) -> Tuple[pd.DataFrame, pd.DataFrame, int, np.ndarray, np.ndarray, pd.DataFrame]:
        """
        Prepares data for the mismatch types summary plots.

        Returns:
            Tuple[pd.DataFrame, pd.DataFrame, int, np.ndarray, np.ndarray, pd.DataFrame]: Processed data for different mismatch types.
        """   
        # prepare data for the confusion matrix
        matrix_data = self.data.loc[(self.data["ref_base"] != "N") & (self.data["majority_base"] != "N"),["ref_base", "majority_base"]]
        matrix_data = pd.crosstab(matrix_data['ref_base'], matrix_data['majority_base'])
        matrix_labels = matrix_data.copy()
        for i in range(matrix_data.shape[0]):
            matrix_labels.iloc[i,i] = None
        matrix_max_mism = matrix_labels.max().max()
        matrix_labels = np.round(matrix_labels / matrix_labels.sum().sum() * 100, 2).astype(str)
        matrix_labels = (matrix_labels + "%").replace("nan%", "")

        # prepare data for the pie chart
        pie_data = matrix_data.reset_index().melt(id_vars="ref_base", value_name="count").sort_values("count", ascending=False)
        pie_data = pie_data.loc[pie_data["ref_base"] != pie_data["majority_base"]].reset_index(drop=True)
        pie_labels = np.array([f"{a} - {b}" for a,b in zip(pie_data["ref_base"], pie_data["majority_base"])])
        pie_data = pie_data["count"].to_numpy()

        # prepare data for error type distributions boxplots
        box_data = self.data.loc[(self.data["ref_base"] != self.data["majority_base"]) & (self.data["ref_base"] != "N")].rename(columns={self.perc_mis_col: "Mismatch", "n_del_rel": "Deletion", "n_ins_rel": "Insertion", "n_ref_skip_rel": "Reference skip"})
        box_data = box_data.melt(value_vars=["Mismatch", "Deletion", "Insertion", "Reference skip"], var_name="Error type", id_vars=["ref_base", "majority_base"])
        box_data["mismatch_type"] = [f"{i} - {j}" for i,j in zip(box_data["ref_base"], box_data["majority_base"])]

        return matrix_data, matrix_labels, matrix_max_mism, pie_data, pie_labels, box_data

    def prepare_data_motifs(self):
        """
        Prepares data for the motif summary plot.

        Returns:
            pd.DataFrame: Processed data for motif summary.
        """
        data = self.data[["ref_base", "motif", self.perc_mis_col, "n_del_rel", "n_ins_rel", "n_ref_skip_rel"]]
        motif_center_idx = len(data.motif[0]) //2
        data.loc[:, "motif_3b"] = data["motif"].map(lambda x: x[motif_center_idx-1:motif_center_idx+2])
        return data

    ##############################################################################################################
    #                                             Plotting functions                                             #
    ##############################################################################################################

    def create_error_placeholder(self, e: Exception):
        """
        Creates a placeholder plot in case there is an error during creation.
        Displays the error message.
        """
        hs.print_update(f"An error occured: {str(e)}. Replacing plot with empty placeholder.")
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
        Creates and returns an HTML representation of the general summary plot.

        Returns:
            str: HTML representation of the plot.
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
        Creates and returns an HTML representation of the chromosome-specific summary plot.

        Returns:
            str: HTML representation of the plot.
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
                                shared_xaxes=True)
            fig = self.update_plot(fig, height=1000, width=1400)

            fig.add_trace(go.Bar(x=n_pos_x, y=n_pos_y))
            fig.add_trace(go.Box(x=chrom, y=n_reads), row=2, col=1)
            fig.add_trace(go.Box(x=chrom, y=quality), row=3, col=1)

            fig.update_traces(marker=dict(color='lightgrey', line=dict(color='black', width=2)), selector=dict(type='bar'))
            fig.update_traces(marker=dict(color="black", outliercolor="black", size=2), fillcolor="lightgrey", selector=dict(type='box'))
            
            fig.update_xaxes(categoryorder='array', categoryarray=sorted(n_pos_x, key=custom_sort_key))
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
        Creates and returns an HTML representation of the general mismatch summary plot.

        Returns:
            str: HTML representation of the plot.
        """
        try:
            def boxplot(data: np.ndarray, group: str = "match", showlegend: bool = False):
                if group == "match":
                    name = "Match"
                    fillcolor = "#4c72b0"
                    legendgroup = "match"
                else:
                    name = "Mismatch"
                    fillcolor = "#dd8452"
                    legendgroup = "mism"
                return go.Box(y=data, name=name, marker=dict(color="black", outliercolor="black", size=2), 
                            fillcolor=fillcolor, showlegend=showlegend, legendgroup=legendgroup)

            data_processed = self.prepare_data_mism_general()

            fig = make_subplots(rows=1, cols=5, 
                                specs=[[{"type": "box"}, {"type": "box"}, {"type": "box"}, {"type": "box"}, {"type": "pie"}]], 
                                shared_yaxes=True,
                                horizontal_spacing=0.01,
                                column_titles=("Mismatch frequency", "Deletion frequency", "Insertion frequency", "Reference skip frequ.", None))
            fig = self.update_plot(fig, height=600, width=1400)

            for i, (dname, legend) in enumerate(zip(["mis_rate", "del_rate", "ins_rate", "ref_skip_rate"], [True, False, False, False]), start=1):
                box_mat = boxplot(data_processed[dname][0], group="match", showlegend=legend)
                box_mis = boxplot(data_processed[dname][1], group="mismatch", showlegend=legend)
                fig.add_traces([box_mat, box_mis], rows=1, cols=i)

            pie = go.Pie(labels=["Match", "Mismatch"], 
                        values=[data_processed["overall"][0], data_processed["overall"][1]], 
                        hoverinfo='label+percent', 
                        textinfo='value', 
                        textfont_size=20, 
                        marker=dict(colors=['#4c72b0', '#dd8452'], line=dict(color='#000000', width=2)), 
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
        Creates and returns HTML representations of mismatch types summary plots.

        Returns:
            List[str]: List of HTML representations of the plots.
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
            fig.update_xaxes(fixedrange=True)
            fig.update_yaxes(fixedrange=True)
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

            for err_type, c in zip(["Mismatch", "Deletion", "Insertion", "Reference skip"], ["#55a868", "#c44e52", "#8172b3", "#937860"]):
                d = box_data.loc[box_data["Error type"] == err_type]
                fig.add_trace(go.Box(x=d["mismatch_type"], y=d["value"], name=err_type, 
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
        Creates and returns an HTML representation of the motif summary plot.

        Returns:
            str: HTML representation of the plot.
        """
        try:
            def create_motif_trace(data: pd.DataFrame, center_base: str, showlegend: bool = False):
                data = data.loc[(data.ref_base == center_base) & (~data.motif_3b.isin(["AAA", "AAC", "AAG", "AAT", "CAA", "GAA", "TAA", "CCC", "CCA", "CCG", "CCT", "ACC", "GCC", "TCC", "GGG", "GGA", "GGC", "GGT", "AGG", "CGG", "TGG", "TTT", "TTA", "TTC", "TTG", "ATT", "CTT", "GTT"]))]
                data = data.melt(value_vars=[self.perc_mis_col, "n_del_rel", "n_ins_rel"], var_name="Error type", id_vars=["ref_base","motif_3b"])

                d = []
                for motif in data["motif_3b"].unique():
                    if self.n_bins:
                        d_m_tmp = pd.DataFrame({"value": self.bin_data(data.loc[(data["motif_3b"]==motif) & (data["Error type"]==self.perc_mis_col), "value"], num_segments=self.n_bins),
                                                "motif": motif,
                                                "Error type": "Mismatch"})
                        d_d_tmp = pd.DataFrame({"value": self.bin_data(data.loc[(data["motif_3b"]==motif) & (data["Error type"]=="n_del_rel"), "value"], num_segments=self.n_bins),
                                                "motif": motif,
                                                "Error type": "Deletion"}) 
                        d_i_tmp = pd.DataFrame({"value": self.bin_data(data.loc[(data["motif_3b"]==motif) & (data["Error type"]=="n_ins_rel"), "value"], num_segments=self.n_bins),
                                                "motif": motif,
                                                "Error type": "Insertion"}) 
                    else:
                        d_m_tmp = pd.DataFrame({"value": data.loc[(data["motif_3b"]==motif) & (data["Error type"]==self.perc_mis_col), "value"],
                                                "motif": motif,
                                                "Error type": "Mismatch"})
                        d_d_tmp = pd.DataFrame({"value": data.loc[(data["motif_3b"]==motif) & (data["Error type"]=="n_del_rel"), "value"],
                                                "motif": motif,
                                                "Error type": "Deletion"}) 
                        d_i_tmp = pd.DataFrame({"value": data.loc[(data["motif_3b"]==motif) & (data["Error type"]=="n_ins_rel"), "value"],
                                                "motif": motif,
                                                "Error type": "Insertion"}) 
                    d += [d_m_tmp, d_d_tmp, d_i_tmp]
                
                d = pd.concat(d)
                
                traces = []
                fillcolors = ["#55a868", "#c44e52", "#8172b3", "#937860"]
                for i, er_type in enumerate(d["Error type"].unique()):
                    d_err = d.loc[d["Error type"] == er_type]
                    trace = go.Box(x=d_err["motif"], y=d_err["value"], name=er_type, 
                                line=dict(color="black", width=1), 
                                marker=dict(outliercolor="black", size=1), 
                                fillcolor=fillcolors[i],
                                legendgroup=i,
                                showlegend=showlegend,
                                offsetgroup=i)
                    traces.append(trace)

                return traces

            data = self.prepare_data_motifs()
            
            x_order = {"A": ["CAC", "CAG", "CAT", "GAC", "GAG", "GAT", "TAC", "TAG", "TAT"], 
                    "C": ["ACA", "ACG", "ACT", "GCA", "GCG", "GCT", "TCA", "TCG", "TCT"],
                    "G": ["AGA", "AGC", "AGT", "CGA", "CGC", "CGT", "TGA", "TGC", "TGT"],
                    "T": ["ATA", "ATC", "ATG", "CTA", "CTC", "CTG", "GTA", "GTC", "GTG"]}

            fig = make_subplots(rows=2, cols=2, 
                                specs=[[{"type": "box"}, {"type": "box"}], [{"type": "box"}, {"type": "box"}]], 
                                shared_yaxes=True, 
                                vertical_spacing=0.05, horizontal_spacing=0.05)
            fig = self.update_plot(fig, height=900, width=1400)

            fig.add_traces(create_motif_trace(data, "A", showlegend=True), rows=1, cols=1)
            fig.add_traces(create_motif_trace(data, "C"), rows=1, cols=2)
            fig.add_traces(create_motif_trace(data, "G"), rows=2, cols=1)
            fig.add_traces(create_motif_trace(data, "T"), rows=2, cols=2)

            fig.update_layout(boxmode="group", boxgroupgap=0, legend=dict(orientation="h", yanchor="bottom", y=1.02, xanchor="right", x=1, bgcolor='#f5f5f5', bordercolor='#000000', borderwidth=2))

            fig.update_xaxes(categoryorder='array', categoryarray=x_order["A"], row=1, col=1)
            fig.update_xaxes(categoryorder='array', categoryarray=x_order["C"], row=1, col=2)
            fig.update_xaxes(title = "3bp Motif", categoryorder='array', categoryarray=x_order["G"], row=2, col=1)
            fig.update_xaxes(title = "3bp Motif", categoryorder='array', categoryarray=x_order["T"], row=2, col=2)

            fig.update_yaxes(title = "Error rate", row=1, col=1)
            fig.update_yaxes(title = "Error rate", row=2, col=1)

        except Exception as e:
            fig = self.create_error_placeholder(e)            
        return fig


    ################################################################################################################
    #                                               Create HTML file                                               #
    ################################################################################################################
    def write_svg(self, fig: go.Figure, name: str) -> None:
        outpath = f"{os.path.splitext(self.output_path)[0]}_{name}.svg"
        fig.write_image(outpath)

    def figs_to_str(self, plot_figs: List[go.Figure]) -> List[str]:
        plot_str = list(map(lambda x: to_html(x, include_plotlyjs=False), plot_figs))
        return plot_str

    def write_to_html(self, n_positions, n_chr, plot_figs: List[go.Figure]) -> None:
        """
        Generate an HTML document containing summary information and plots created from extracted data.

        Args:
            n_positions (int): Total number of extracted positions.
            n_chr (int): Total number of sequences (chromosomes) analyzed.
            plots (List[str]): List of HTML strings containing plot visualizations.

        Returns:
            None: The method generates an HTML file based on the provided template and plots,
                and writes it to the specified output path.
        """
        if self.export_svg:
            for fig, name in zip(plot_figs, ["summary_general_info", "summary_chr_info", "summary_mismatch_stats", 
                                             "summary_mismatch_matrix", "summary_mismatch_pie", 
                                             "summary_error_rates_by_mismatch", "summary_error_rates_by_motif"]):
                self.write_svg(fig, name)
        plots = self.figs_to_str(plot_figs)

        time = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
        data_point_descr = f"Each data point corresponds to the average value along all positions in one of {self.n_bins} bins" if self.n_bins else "Each data point corresponds to one extracted position"
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
                            <h2 class="collapsible-header">General statistics</h2>

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
                        </section>

                        <section>
                            <h2>Mismatch statistics</h2>

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
                        </section>
                        
                        <section>
                            <h2>Error rate by motifs</h2>
                            
                            <h3>Mismatch, insertion and deletion rates for 3bp motifs</h3>
                            <div class="plot-container">
                                {plots[6]}
                            </div>
                            <p>
                                Distributions of Mismatch, deletion and insertion rates for different three base motifs with center A (<b>top left</b>), C (<b>top right</b>),
                                G (<b>bottom left</b>) and T (<b>bottom right</b>). {data_point_descr}.
                            </p>

                        </section>
                    </body>
                    <footer></footer>
                    </html> 
                    """
        with open(self.output_path, "w") as o:
            o.write(template)

def setup_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="Neet - summary creator", description="Create overview in HTML format containing interactive plots.")
    parser.add_argument('-i', '--input', type=str, required=True,
                        help="""
                            Path to the input TSV file. 
                            """)
    parser.add_argument('-o', '--output', type=str, required=True,
                        help="""
                            Path to output a output directory, in which all output files will be stored.
                            """)
    parser.add_argument('-b', '--n_bins', type=int, required=False, default=5000,
                        help="""Number of bins to split the data into when creating the summary plots. This does not affect the extracted data.
                            Used only to improve performance and clarity of the created plots. Note that setting the value to a low number 
                            can lead to misleading results. Set to '-1' to disable binning. Default: 5000
                            """)
    parser.add_argument('--plot_alt', action="store_true", 
                        help="""
                            Specify whether to use the perc_mismatch or perc_mismatch_alt values for plot creation.
                            """)
    parser.add_argument('--export_svg', action="store_true", 
                        help="""
                            Specify whether to export the created plots as svg files.
                            """)

    return parser


if __name__=="__main__":

    parser = setup_parser()
    args = parser.parse_args()
    
    sc = SummaryCreator(args.input, args.output, args.n_bins, args.plot_alt, export_svg=args.export_svg)
    sc.create_summary()

    # sc = SummaryCreator("/home/vincent/masterthesis/data/45s_rrna/processed/45s_cytoplasm_extracted.tsv", "/home/vincent/masterthesis/data/45s_rrna/processed/", None)