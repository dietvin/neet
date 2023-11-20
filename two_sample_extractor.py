from typing import List, Dict, Any, Tuple, Set
import os, warnings, sys, datetime, argparse
from collections import Counter
import helper_functions as hs
import plotly.graph_objects as go
from plotly.io import to_html

class PositionExtractor:
    out_dir: str
    export_svg: bool

    in_paths_a: List[str]
    in_paths_b: List[str]

    label_a: str
    label_b: str

    error_feature_idx: int
    error_threshold: float
    coverage_threshold: int

    ref_dict: Dict[str, str]

    def __init__(self, in_paths_a: str, in_paths_b: str, out_dir: str, ref_path: str, error_feature: str, error_threshold: float, coverage_threshold: int, 
                 label_a: str, label_b: str, export_svg: bool = False) -> None:
        """
        Initialize the PositionExtractor object with input parameters.

        Parameters:
        - in_paths_a (str): Comma-separated string of file paths for dataset A.
        - in_paths_b (str): Comma-separated string of file paths for dataset B.
        - out_dir (str): Output directory for the extracted data.
        - ref_path (str): File path to the reference data in FASTA format.
        - error_feature (str): The type of error feature to be analyzed.
        - error_threshold (float): Threshold for the specified error feature.
        - coverage_threshold (int): Threshold for the coverage of positions to be considered.
        - label_a (str): Label for dataset A.
        - label_b (str): Label for dataset B.
        - export_svg (bool, optional): Flag indicating whether to export SVG plots (default is False).

        Raises:
        - Exception: If an invalid error feature is provided.
        """
        self.in_paths_a = self.process_in(in_paths_a)
        self.in_paths_b = self.process_in(in_paths_b)

        self.check_out(out_dir)

        if not out_dir.endswith("/"): out_dir+="/"
        self.out_dir = out_dir
        self.export_svg = export_svg

        self.label_a = label_a
        self.label_b = label_b
        
        try:
            feature_idx = {"n_del_rel": 16,
                           "n_ins_rel": 17,
                           "perc_mismatch": 19,
                           "perc_mismatch_alt": 20}
            self.error_feature_idx = feature_idx[error_feature]
        except KeyError:
            raise Exception(f"Invalid error feature '{error_feature}'. Use one of: n_del_rel, n_ins_rel, perc_mismatch, perc_mismatch_alt")

        self.error_threshold = error_threshold
        self.coverage_threshold = coverage_threshold

        self.ref_dict = self.get_references(ref_path)

    ##############################################################################################################
    #                                           Initialization methods                                           #
    ##############################################################################################################

    def process_in(self, in_paths: str) -> List[str]:
        """
        Process a comma-separated string of input file paths.

        Parameters:
        - in_paths (str): A comma-separated string of input file paths.

        Returns:
        - in_list (List[str]): A list of validated input file paths.

        Raises:
        - Exception: If any input file does not exist or has an unexpected file extension.
        - Exception: If input files are of different kinds (either .bam or .pileup). All files must be of the same kind.
        """
        in_list = in_paths.split(",")
        extensions = []
        for path in in_list:
            ext = os.path.splitext(path)[1]
            extensions.append(ext)
            self.check_path(path, extension=".tsv")

        if len(set(extensions)) > 1:
            raise Exception("Input files of different kind. All files must be .bam or .pileup, not mixed.")
        return in_list

    def check_out(self, out: str): 
        """
        Check if the specified output directory exists, and create it if necessary.

        Parameters:
        - out (str): The path to the output directory.

        Raises:
        - Warning: If the specified directory does not exist, a warning is issued, and the directory is created.
        - Exception: If the directory could not be created.
        """
        if not os.path.isdir(out):
            warnings.warn(f"Directory '{out}' does not exist. Creating directory.")
            try:
                os.makedirs(out)
            except:
                raise Exception(f"Directory '{out}' was not found and could not be created.")

    def check_path(self, path: str, extension: str = ".tsv") -> None:
        """
        Check if the specified file path exists and has the expected file extension.

        Parameters:
        - path (str): The file path to be checked.
        - extension (str, optional): The expected file extension (default is ".tsv").

        Raises:
        - FileNotFoundError: If the specified file path does not exist.
        - Exception: If the file extension of the specified path does not match the expected extension.
        """        
        if not os.path.exists(path): # does file exist?
            raise FileNotFoundError(f"Input file not found. File '{path}' does not exist.")
        file_type = os.path.splitext(path)[1]
        if file_type != extension:
            raise Exception(f"Found file extension {file_type}. File extension to must be {extension}.")

    def get_references(self, path: str) -> Dict[str, str]:
        """
        Extract reference sequences from a FASTA file.

        Parameters:
        - path (str): The path to the FASTA file.

        Returns:
        - refs (Dict[str, str]): A dictionary where keys are chromosome names and values are corresponding sequences.

        Raises:
        - Exception: If there is a format error in the FASTA file, such as a missing header in the first line.
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

    ##############################################################################################################
    #                                             Processing methods                                             #
    ##############################################################################################################

    def extract_positions(self) -> None:
        """
        Extract positions based on coverage and error thresholds, identify systematic positions,
        and categorize them into overlapping and exclusive sets for datasets A and B.
        """
        subsets_a = [self.extract_suffient_cvg_pos(path) for path in self.in_paths_a]
        subsets_b = [self.extract_suffient_cvg_pos(path) for path in self.in_paths_b]

        a_high_cvg = [pos[0] for pos in subsets_a]
        a_high_err = [pos[1] for pos in subsets_a]
        del(subsets_a)

        b_high_cvg = [pos[0] for pos in subsets_b]
        b_high_err = [pos[1] for pos in subsets_b]
        del(subsets_b)

        a_high_cvg_syst = self.find_systematic_pos(a_high_cvg)
        del(a_high_cvg)
        a_high_err_syst = self.find_systematic_pos(a_high_err)
        del(a_high_err)

        b_high_cvg_syst = self.find_systematic_pos(b_high_cvg)
        del(b_high_cvg)
        b_high_err_syst = self.find_systematic_pos(b_high_err)
        del(b_high_err)

        a_and_b = set.intersection(a_high_err_syst, b_high_err_syst)
        a_excl = set.intersection(a_high_err_syst, b_high_cvg_syst-a_and_b)
        b_excl = set.intersection(b_high_err_syst, a_high_cvg_syst-a_and_b)

        self.a_and_b = sorted(a_and_b)
        self.a_excl = sorted(a_excl)
        self.b_excl = sorted(b_excl)

    def extract_suffient_cvg_pos(self, path: str) -> Tuple[List[Tuple[str, int]], List[Tuple[str, int]]]:
        """
        Extract positions with sufficient coverage and error rate from the specified input file.

        Parameters:
        - path (str): Path to the input file.

        Returns:
        - Tuple[List[Tuple[str, int]], List[Tuple[str, int]]]: Two lists of positions - one for high coverage and one for high error rate.

        Note:
        - Uses the error feature index specified during object initialization.
        """
        hs.print_update(f"Processing {path} ... ", line_break=False)
        sys.stdout.flush()
        with open(path, "r") as file:
            high_cvg_sites = []
            high_cvg_high_err_sites = []

            next(file)
            for line in file:
                line = line.split("\t")
                chrom, site, n_reads, n_ref_skip, error_rate = line[0], int(line[1]), int(line[2]), int(line[11]), float(line[self.error_feature_idx])            
                if n_reads-n_ref_skip >= self.coverage_threshold:
                    high_cvg_sites.append((chrom, site)) 
                    if error_rate >= self.error_threshold:
                        high_cvg_high_err_sites.append((chrom, site))

            hs.print_update("done", with_time=False)

            return high_cvg_sites, high_cvg_high_err_sites

    def find_systematic_pos(self, replicate_pos: List[List[Tuple[str, int]]]) -> List[Tuple[str, int]]:
        """
        Identify systematic positions that are common across multiple replicates.

        Parameters:
        - replicate_pos (List[List[Tuple[str, int]]]): List of positions from multiple replicates.

        Returns:
        - List[Tuple[str, int]]: List of systematic positions.
        """
        replicate_pos = [set(pos) for pos in replicate_pos]
        return set.intersection(*replicate_pos)

    def write_bed_files(self) -> None:
        """Write BED files for overlapping and exclusive positions."""
        self.positions_to_bed(self.a_and_b, self.out_dir, f"{self.label_a}_{self.label_b}")
        self.positions_to_bed(self.a_excl, self.out_dir, f"{self.label_a}_excl")
        self.positions_to_bed(self.b_excl, self.out_dir, f"{self.label_b}_excl")

    def positions_to_bed(self, positions: Set[Tuple[str, int]], out_dir: str, name: str) -> None:
        """
        Write positions to a BED file.

        Parameters:
        - positions (Set[Tuple[str, int]]): Set of positions to be written.
        - out_dir (str): Output directory.
        - name (str): Name for the output BED file.
        """
        path = f"{out_dir}{name}.bed"
        hs.print_update(f"Writing {len(positions)} sites to {path}")
        with open(path, "w") as out:
            for position in positions:
                chrom, site = position[0], position[1]
                out.write(f"{chrom}\t{site-1}\t{site}\t{name}\n")

    

    ##############################################################################################################
    #                                          Summary creation methods                                          #
    ##############################################################################################################

    def custom_sort_key(self, item):
        """
        Define a custom sorting key for sorting mixed alphanumeric items.

        This method defines a custom sorting key that can be used with the `sorted()` function or other sorting functions.
        It allows for sorting a list of mixed alphanumeric items in a way that numeric values are sorted numerically,
        and non-numeric values are placed at the end of the sorted list.

        Parameters:
            item (str): The item to be sorted.

        Returns:
            tuple: A tuple used as a sorting key. The tuple consists of two elements:
            - If the item is numeric, its integer value is placed first to sort numerically.
            - If the item is not numeric, it is placed last with a float('inf') value to ensure it comes after numeric values.

        Example:
            To sort a list of mixed alphanumeric items, you can use this custom sorting key as follows:
            ```
            sorted_list = sorted(my_list, key=obj.custom_sort_key)
            ```

        """
        if item.isdigit():  # Check if the item is a digit
            return (int(item),)  # Convert to integer and sort numerically
        else:
            return (float('inf'), item)  # Place non-digits at the end


    def create_plot(self) -> go.Figure:
        """
        Create a Plotly figure representing the chromosome-wise distribution of positions.

        Returns:
        - go.Figure: The Plotly figure.
        """
        def update_plot(fig, title: str|None = None, xlab: str|None = None, ylab: str|None = None, height: int = 500, width: int = 800):
            """
            Updates the layout of a Plotly figure.

            Parameters:
            - fig (go.Figure): The Plotly figure to be updated.
            - title (str | None): Title of the plot (optional).
            - xlab (str | None): Label for the x-axis (optional).
            - ylab (str | None): Label for the y-axis (optional).
            - height (int): Height of the plot (default: 500).
            - width (int): Width of the plot (default: 800).

            Returns:
            - go.Figure: The updated Plotly figure.
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
        
        chr_a_and_b = [pos[0] for pos in self.a_and_b]
        sites_a_and_b = [pos[1] for pos in self.a_and_b]

        chr_a_excl = [pos[0] for pos in self.a_excl]
        sites_a_excl = [pos[1] for pos in self.a_excl]

        chr_b_excl =[pos[0] for pos in self.b_excl]
        sites_b_excl =[pos[1] for pos in self.b_excl]

        all_keys = self.ref_dict.keys()
        present_chr = set((all_keys & chr_a_and_b) | (all_keys & chr_a_excl) | (all_keys & chr_b_excl))
        present_chr = list(sorted(present_chr, key=self.custom_sort_key))
        chr_lens = [len(self.ref_dict[x]) for x in present_chr]

        width = 1200
        bargap = 0.9
        bar_width = 1-bargap+0.15
        n_bars = len(present_chr)
        scatter_size = width/n_bars*bar_width

        fig = update_plot(go.Figure(), height = 1000, width = width)
        fig.add_trace(go.Bar(x=present_chr, y=chr_lens, marker=dict(color="lightgrey", line=dict(color="black", width=2)), name="Chromosomes", showlegend=False))
        fig.update_layout(bargap=0.5, yaxis=dict(range=[0,max(chr_lens)+0.1*max(chr_lens)]))

        fig.add_trace(go.Scatter(x=chr_a_and_b, y=sites_a_and_b, mode='markers', marker=dict(symbol='line-ew', color="#1f77b4", size=scatter_size, line=dict(width=1.1, color="#1f77b4")), name=f"{self.label_a}+{self.label_b}", hovertemplate="Chr%{x}:%{y}"))
        fig.add_trace(go.Scatter(x=chr_a_excl, y=sites_a_excl, mode='markers', marker=dict(symbol='line-ew', color="#ff7f0e", size=scatter_size, line=dict(width=1.1, color="#ff7f0e")), name=f"{self.label_a}", hovertemplate="Chr%{x}:%{y}"))
        fig.add_trace(go.Scatter(x=chr_b_excl, y=sites_b_excl, mode='markers', marker=dict(symbol='line-ew', color="#2ca02c", size=scatter_size, line=dict(width=1.1, color="#2ca02c")), name=f"{self.label_b}", hovertemplate="Chr%{x}:%{y}"))
        fig.update_xaxes(fixedrange=True)

        return fig

    def get_file_paths(self) -> Tuple[str, str]:
        """
        Get formatted HTML lists of input file paths for datasets A and B.

        Returns:
        - Tuple[str, str]: Formatted HTML lists of file paths for datasets A and B.
        """
        def create_list(paths: List[str]):
            path_list = "<ul>"
            for path in paths:
                path_list += f"<li>{path}</li>"
            path_list += "</ul>"
            return path_list
        
        return create_list(self.in_paths_a), create_list(self.in_paths_b)

    def create_count_table(self) -> str:
        """
        Create an HTML table showing the counts of occurrences of chromosomes in each set.

        Returns:
        - str: The HTML representation of the count table.
        """
        # Count occurrences of chromosomes in each set
        a_and_b_counts = Counter(chromosome for chromosome, _ in self.a_and_b)
        a_excl_counts = Counter(chromosome for chromosome, _ in self.a_excl)
        b_excl_counts = Counter(chromosome for chromosome, _ in self.b_excl)

        # Get all unique chromosomes from the sets
        all_chromosomes = set(a_and_b_counts.keys()) | set(a_excl_counts.keys()) | set(b_excl_counts.keys())

        # Create HTML table
        html_table = """
        <table border="1">
        <thead><tr>
            <th></th>
        """

        # Add column headers
        for chromosome in all_chromosomes:
            html_table += f"<th>{chromosome}</th>"
        html_table += "<th>Total</th></tr></thead><tbody>"

        # Add rows and counts
        for set_name, counts in [(f'{self.label_a} + {self.label_b}', a_and_b_counts), (f'{self.label_a} excl.', a_excl_counts), (f'{self.label_b} excl.', b_excl_counts)]:
            html_table += f"<tr><td><b>{set_name}</b></td>"
            total_count = sum(counts.values())
            for chromosome in all_chromosomes:
                html_table += f"<td>{counts[chromosome]}</td>"
            html_table += f"<td>{total_count}</td></tr>"

        html_table += "</tbody></table>"
        return html_table
    
    def write_svg(self, fig: go.Figure, name: str) -> None:
        """
        Write a Plotly figure to an SVG file.

        Parameters:
        - fig (go.Figure): The Plotly figure to be written.
        - name (str): The name of the SVG file.
        """
        outpath = f"{self.out_dir}/{name}.svg"
        fig.write_image(outpath)

    def create_summary(self) -> None:
        """
        Create a summary HTML file containing an overview of shared and exclusive positions in samples A and B.

        This method generates a summary HTML file that includes information about the input files, a count table, and a genome map.

        If export_svg is True, the genome map is saved as an SVG file in the specified output directory.

        The HTML file is saved with the format: "<out_dir>/<label_a>_<label_b>_summary.html".

        Returns:
        - None
        """
        plot = self.create_plot()
        if self.export_svg:
            self.write_svg(plot, "twosample_map")
        
        plot = to_html(plot, include_plotlyjs=False)

        paths_1, paths_2 = self.get_file_paths()
        count_table = self.create_count_table()
        time = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')


        template = f"""
                    <!DOCTYPE html>
                    <html lang="en">
                    <head>
                        <meta charset="UTF-8">
                        <meta http-equiv="X-UA-Compatible" content="IE=edge">
                        <meta name="viewport" content="width=device-width, initial-scale=1.0">
                        <title>Neet - Position extractor summary</title>
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
                            ol {{
                                padding-left: 1em;
                                margin-left: 3.5em;
                                margin-top: 0.5em;
                                margin-bottom: 0.5em;
                            }}
                            
                            ul {{
                                padding-left: 1em;
                                margin-left: 3.5em;
                                margin-top: 0.5em;
                                margin-bottom: 0.5em;
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
                                margin-bottom: 1em;
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
                                margin-top: 0.5em;
                                margin-bottom: 0.5em;
                                margin-left: 1.5em;
                                margin-right: 1.5em;
                            }}

                            .intro-text ul {{
                                font-size: 1.3rem;
                                margin-left: 3.5em;
                                margin-top: 0.5em;
                                margin-bottom: 0.5em;
                            }}

                            .table-box {{
                                margin: 2em;
                                border: 1px solid #ccc;
                                overflow-x: auto; /* Enable horizontal scrolling */
                            }}

                            table {{
                                width: 100%;
                                border-collapse: collapse;
                                border: 1px solid black;
                                font-family: Arial, sans-serif;
                            }}

                            th, td {{
                                border: 1px solid #ccc;
                                padding: 8px;
                                text-align: center;

                            }}

                            thead {{
                                background-color: #333;
                                color: white;
                                font-size: 0.875rem;
                                text-transform: uppercase;
                                letter-spacing: 2%;
                            }}

                            tbody tr:nth-child(odd) {{
                                background-color: #fff;
                            }}

                            tbody tr:nth-child(even) {{
                                background-color: #eee;
                            }}
                        </style>
                    </head>

                    <body>
                        <script src="https://cdn.plot.ly/plotly-latest.min.js"></script>

                        <header>
                            <h1>Position extractor summary</h1>
                            <p>Produced by <a href="https://github.com/dietvin/neet">Neet</a> on <b>{time}</b></p>
                        </header>
                    
                        <section>
                            <p class="intro-text">
                                This summary file contains an overview of the shared and exclusive positions in samples <i>{self.label_a}</i> and <i>{self.label_b}</i>.
                                Files provided for sample <i>{self.label_a}</i>:
                            </p>
                            {paths_1}
                            <p class="intro-text">Files provided for sample <i>{self.label_b}</i>:</p>
                            {paths_2}
                        </section>

                        <section>
                            <h2 class="collapsible-header">General information</h2>
                            <h3>Count table</h3>
                            <p>
                                Positions are only regarded if they are present in all replicates of a sample. The count table displays the number of individual genomic positions for different groups in each row. These groups are: 
                                <ol>
                                    <li><i>{self.label_a}</i> + <i>{self.label_b}</i>:       Positions that are shared between samples {self.label_a} and {self.label_b}</li>
                                    <li><i>{self.label_a}</i> excl.:              Positions that are exclusive to sample {self.label_a}</li>
                                    <li><i>{self.label_b}</i> excl.:              Positions that are exclusive to sample {self.label_b}</li>
                                </ol>
                            </p>
                            <div class="table-box">
                                {count_table}
                            </div>
                        </section>

                        <section>
                            <h2 class="collapsible-header">Genome map</h2>
                            <h3>Genome map displaying shared and exclusive positions between samples <i>{self.label_a}</i> and <i>{self.label_b}</i></h3>
                            <div class="plot-container">
                                {plot}
                            </div>
                            <p>
                                Genome map displaying the positions that are present in both samples, only in sample <i>{self.label_a}</i> and only in sample <i>{self.label_b}</i>.
                                Only chromosomes that contain extracted position(s) are shown.  
                            </p>
                        </section>
                    </body>
                    <footer></footer>
                    </html> 
            """
        with open(f"{self.out_dir}{self.label_a}_{self.label_b}_summary.html", "w") as out:
            out.write(template)

    ##############################################################################################################
    #                                                Main method                                                 #
    ##############################################################################################################

    def main(self) -> None:
        """
        Execute the main workflow of the PositionExtractor class.

        This method performs the following tasks:
        1. Extract positions from input files.
        2. Write BED files for shared and exclusive positions.
        3. Create a summary HTML file containing an overview of shared and exclusive positions.

        Returns:
        - None
        """
        self.extract_positions()
        self.write_bed_files()
        self.create_summary()




def setup_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="Neet - position extractor", description="Identify and map positions of interest between two samples.")
    parser.add_argument('-i', '--sample1', type=str, required=True,
                        help="""
                            Path to the input file(s). If replicates are available, specify paths comma-separated (<repl1.tsv>,<repl2.tsv>,...).
                            Must be of type tsv, as returned by the PileupExtractor.
                            """)
    parser.add_argument('-bn', '--basename', type=str, default="sample1",
                        help="""
                            Basename of the given sample. Used to create the pileup and extracted features files. Default: 'sample1'
                            """)
    parser.add_argument('-i2', '--sample2', type=str, required=True,
                        help="""
                            Path to the input file(s) from a second sample. If replicates are available, specify paths comma-separated (<repl1.tsv>,<repl2.tsv>,...).
                            Must be of type tsv, as returned by the PileupExtractor. 
                            """)
    parser.add_argument('-bn2', '--basename2', type=str, default="sample2",
                        help="""
                            Basename of the second sample. Used to create the pileup and extracted features files. Default: 'sample2'
                            """)
    parser.add_argument('-o', '--output', type=str, required=True,
                        help="""
                            Path to output a output directory, in which all output files will be stored.
                            """)
    parser.add_argument('-r', '--reference', type=str, required=True, help="Path to the reference file")
    parser.add_argument("-f", "--error_feature", type=str, default="perc_mismatch_alt", 
                        help="""
                            Error feature to use during extraction. Can be one of the following: 
                            n_del_rel, n_ins_rel, perc_mismatch, perc_mismatch_alt.
                            """)
    parser.add_argument("-e", "--error_threshold", type=hs.float_between_zero_and_one, default=0.5, 
                        help="""
                            Threshold to identify positions of iterest. Uses the perc_mismatch_alt feature.
                            """)
    parser.add_argument("-c", "--coverage_threshold", type=hs.positive_int, default=40,
                        help="""
                            Minimum coverage of a position to be regarded in the extraction.
                            """)
    parser.add_argument('--export_svg', action="store_true", 
                        help="""
                            If specified, exports created plots as svg files.
                            """)
    return parser

if __name__ == "__main__":
    parser = setup_parser()
    args = parser.parse_args()

    posextr = PositionExtractor(in_paths_a=args.sample1, in_paths_b=args.sample2, 
                                out_dir=args.output, ref_path=args.reference, 
                                label_a=args.basename, label_b=args.basename2, 
                                error_feature=args.error_feature,
                                error_threshold=args.error_threshold, coverage_threshold=args.coverage_threshold,
                                export_svg=args.export_svg)
    posextr.main()