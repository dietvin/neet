from typing import List, Dict, Tuple, Any
import os, sys, warnings
import numpy as np
import pandas as pd
import plotly.graph_objects as go
import plotly.express as px
from plotly.subplots import make_subplots
from plotly.io import to_html
from intervaltree import IntervalTree
import collections, datetime

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

        self.process_path(in_path, out_path, bed_path, ref_path)
        self.load_data()
        self.perc_mismatch_col = "perc_mismatch_alt" if use_perc_mismatch_alt else "perc_mismatch"
        self.get_bed_categories(categories)
        self.get_corrensponding_base(canonical_counterpart)
        self.output_tsv = output_tsv

    ##############################################################################################################
    #                                             Initialization methods                                             #
    ##############################################################################################################

    def process_path(self, in_path: str, out_path: str, bed_path: str, ref_path: str) -> None:
        self.check_path(in_path, [".tsv"])
        self.in_path = in_path
        self.check_path(bed_path, [".bed"])
        self.bed_path = bed_path
        self.check_path(ref_path, [".fasta", ".fa", ".fn"])
        self.ref_path = ref_path
        self.output_path = self.process_outpath(out_path)

    def check_path(self, path: str, extensions: List[str]) -> None:
        if not os.path.exists(path): # does file exist?
            raise FileNotFoundError(f"Input file not found. File '{path}' does not exist.")
        file_type = os.path.splitext(path)[1]
        if not file_type in extensions:
            warnings.warn(f"Found file extension {file_type}. Expected file extension to be one of: {extensions}. If this is deliberate, ignore warning.", Warning)

    def process_outpath(self, out: str) -> str:
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
            data = pd.read_csv(self.in_path, sep="\t", low_memory=False)
            self.data = self.add_bed_info(data, self.bed_path)
    
    def add_bed_info(self, data: pd.DataFrame, bed_path: str):
        tree = self.build_interval_tree(bed_path)
        data["bed_name"] = data.apply(lambda x: self.get_name_for_position((x[0], x[1]), tree), axis=1)
        return data

    def get_name_for_position(self, position, interval_tree) -> str|None:
        results = interval_tree[position[1]:position[1]+1]  # Query the interval tree
        for interval in results:
            chromosome, name = interval.data
            if chromosome == position[0]:
                return name
        return None

    def build_interval_tree(self, bed_file):
        interval_tree = IntervalTree()
        with open(bed_file, 'r') as file:
            for line in file:
                parts = line.strip().split('\t')
                if len(parts) >= 3:
                    chromosome = parts[0]
                    start = int(parts[1])
                    end = int(parts[2])
                    name = parts[3] if len(parts) >= 4 else None
                    interval_tree[start:end] = (chromosome, name)
        return interval_tree

    def get_bed_categories(self, cat_str: str) -> None:
        categories = cat_str.split(",")
        unique_cat = list(self.data.bed_name.unique())

        for category in categories:
            if category not in unique_cat:
                raise Exception(f"Given name '{category}' was not found in the bed file.")
    
        self.bed_categories = categories

    def get_corrensponding_base(self, base_str: str) -> None:
        """
        Extracts the bases that correspond to each category given for the bed-categories.
        """
        counterparts = base_str.split(",")
        if len(counterparts) != len(self.bed_categories):
            raise Exception(f"For the {len(self.bed_categories)} categories {len(counterparts)} corresponding bases were given. Each category must have a base it corresponds to.")
        self.bed_categories_canoncial_bases = counterparts

    ##############################################################################################################
    #                                             Data preparation methods                                             #
    ##############################################################################################################

    def prepare_data_mism_types(self, data: pd.DataFrame):
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
    #                                             Plotting methods                                             #
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
        data = self.data.loc[self.data.bed_name == mod]
        matrix_data, matrix_labels, matrix_max_mism = self.prepare_data_mism_types(data)

        fig = px.imshow(matrix_data, labels=dict(x="Called base", y="Reference base", color="Count"), zmin=0, zmax=1.2*matrix_max_mism, color_continuous_scale="portland")
        fig = self.update_plot(fig, None, "Called base", "Reference base", width=800)
        fig.update_traces(text=matrix_labels, texttemplate="%{text}")
        fig.update_xaxes(fixedrange=True)
        fig.update_yaxes(fixedrange=True)
        
        return to_html(fig, include_plotlyjs=False)

    def create_composition_plot(self, mod: str, canonical: str) -> str:
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
    #                                             Output methods                                             #
    ##############################################################################################################

    def write_template(self, plots: List[str], category: str, corresponding_base: str) -> None:
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
        self.data.to_csv(f"{os.path.splitext(self.in_path)[0]}_w_bed_info.tsv", sep="\t", index=False, header=True)

    ##############################################################################################################
    #                                             Main processing methods                                             #
    ##############################################################################################################

    def main(self):
        for category, corr_base in zip(self.bed_categories, self.bed_categories_canoncial_bases):
            self.process_category(category, corr_base)

    def process_category(self, category: str, corresponding_base: str):
        plot_mod_map = self.create_map_plot(category)
        plot_mism_types = self.create_mism_types_plot(category)
        plot_comp = self.create_composition_plot(category, corresponding_base)
        plot_err_rate = self.create_error_rate_plot(category, corresponding_base)
        plot_nb = self.create_nb_plot(category)

        plots = [plot_mod_map, plot_mism_types, plot_comp, plot_err_rate, plot_nb]
        self.write_template(plots, category, corresponding_base)
        if self.output_tsv:
            self.write_tsv()

if __name__ == "__main__":
    poi_analyzer = POIAnalyzer("/home/vincent/masterthesis/data/45s_rrna/processed/dRNA_cytoplasm/dRNA_cytoplasm_extracted.tsv",
                               "/home/vincent/masterthesis/data/45s_rrna/processed/dRNA_cytoplasm",
                               "/home/vincent/masterthesis/data/45s_rrna/rRNA_modifications_conv.bed",
                               "/home/vincent/masterthesis/data/45s_rrna/RNA45SN1_cleaned.fasta",
                               "Gm", "G", True, False)
    poi_analyzer.main()
