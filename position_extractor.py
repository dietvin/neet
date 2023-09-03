from typing import List, Dict, Any, Tuple
import os, warnings, sys, datetime, argparse
import pandas as pd
import numpy as np
import plotly.graph_objects as go
from plotly.io import to_html

class PositionExtractor:
    out_dir: str

    in_paths_1: List[str]
    in_paths_2: List[str]

    datasets_1: List[pd.DataFrame]
    datasets_2: List[pd.DataFrame]

    label_1: str
    label_2: str

    ref_dict: Dict[str, str]

    common_pos_1: pd.DataFrame
    common_pos_2: pd.DataFrame

    excl_in_1: pd.DataFrame
    excl_in_2: pd.DataFrame
    in_1_and_2: pd.DataFrame

    dtypes = {'chr': str, 'site': int, 'n_reads': int, 'ref_base': str, 'majority_base': str, 'n_a': int, 'n_c': int,
          'n_g': int, 'n_t': int, 'n_del': int, 'n_ins': int, 'n_ref_skip': int, 'n_a_rel': float, 'n_c_rel': float,
          'n_g_rel': float, 'n_t_rel': float, 'n_del_rel': float, 'n_ins_rel': float, 'n_ref_skip_rel': float,
          'perc_mismatch': float, 'perc_mismatch_alt': float, 'motif': str, 'q_mean': float, 'q_std': float,
          'neighbour_error_pos': str}

    def __init__(self, in_paths_1: str, in_paths_2: str, out_dir: str, ref_path: str, label_1: str|None = None, label_2: str|None = None) -> None:
        self.in_paths_1 = self.process_in(in_paths_1)
        self.in_paths_2 = self.process_in(in_paths_2)
        self.datasets_1 = []
        self.datasets_2 = []
        self.read_datasets(self.in_paths_1, self.in_paths_2)

        self.check_out(out_dir)
        if not out_dir.endswith("/"): out_dir+="/"
        self.out_dir = out_dir

        self.label_1 = label_1 if label_1 else "sample1"
        self.label_2 = label_2 if label_2 else "sample2"

        self.ref_dict = self.get_references(ref_path)

    def process_in(self, in_paths: str) -> List[str]:
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
        if not os.path.isdir(out):
            warnings.warn(f"Directory '{out}' does not exist. Creating directory.")
            try:
                os.makedirs(out)
            except:
                raise Exception(f"Directory '{out}' was not found and could not be created.")

    def check_path(self, path: str, extension: str = ".tsv") -> None:
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
        if file_type != extension:
            raise Exception(f"Found file extension {file_type}. File extension to must be {extension}.")

    def read_datasets(self, in_paths_1: List[str], in_paths_2: List[str]) -> None:
        for path in in_paths_1:
            self.datasets_1.append(pd.read_csv(path, sep="\t", dtype=self.dtypes, usecols=["chr", "site"]))
        for path in in_paths_2:
            self.datasets_2.append(pd.read_csv(path, sep="\t", dtype=self.dtypes, usecols=["chr", "site"]))        

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

    def extract_subsets_positions(self) -> None:
        # find positions that are present in all replicates
        self.replicates_find_common_pos()

        # find positions that are exclusive in dataset 1
        self.excl_in_1 = self.find_a_wo_b(self.common_pos_1, self.common_pos_2)
        # find positions that are exclusive in dataset 2
        self.excl_in_2 = self.find_a_wo_b(self.common_pos_2, self.common_pos_1)
        # find positions that are present in 1 AND 2 
        self.in_1_and_2 = self.find_intersect_pos([self.common_pos_1, self.common_pos_2])

    def replicates_find_common_pos(self) -> None:
        if len(self.datasets_1) > 1:
            self.common_pos_1 = self.find_intersect_pos(self.datasets_1)
        else:
            self.common_pos_1 = self.datasets_1[0]
        if len(self.datasets_2) > 1:
            self.common_pos_2 = self.find_intersect_pos(self.datasets_2)
        else:
            self.common_pos_2 = self.datasets_2[0]
    
    def find_a_wo_b(self, df_a: pd.DataFrame, df_b: pd.DataFrame) -> pd.DataFrame:
        merged = pd.merge(df_a, df_b, on=["chr", "site"], how="left", indicator=True)
        return merged.loc[merged["_merge"]=="left_only"][["chr", "site"]]
    
    def find_intersect_pos(self, datasets: List[pd.DataFrame]) -> pd.DataFrame:
        merged = datasets[0]
        for dataset in datasets[1:]:
            merged = pd.merge(merged, dataset, on=["chr", "site"])
        return merged

    def write_bed_files(self) -> None:
        outpath = f"{self.out_dir}{self.label_1}_exclusive.bed"
        self.excl_in_1.to_csv(outpath, sep="\t", header=False, index=False)

        outpath = f"{self.out_dir}{self.label_2}_exclusive.bed"
        self.excl_in_2.to_csv(outpath, sep="\t", header=False, index=False)

        outpath = f"{self.out_dir}{self.label_1}_and_{self.label_2}.bed"
        self.in_1_and_2.to_csv(outpath, sep="\t", header=False, index=False)

    def create_summary(self) -> None:
        plot = self.create_plot()
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
                                This summary file contains an overview of the shared and exclusive positions in samples <i>{self.label_1}</i> and <i>{self.label_2}</i>.
                                Files provided for sample <i>{self.label_1}</i>:
                            </p>
                            {paths_1}
                            <p class="intro-text">Files provided for sample <i>{self.label_2}</i>:</p>
                            {paths_2}
                        </section>

                        <section>
                            <h2 class="collapsible-header">General information</h2>
                            <h3>Count table</h3>
                            <p>
                                The count table displays the number of individual genomic positions for different groups in each row. These groups are: 
                                <ol>
                                    <li>sytematic <i>{self.label_1}</i>:        Includes positions that are present in all replicates of sample <i>{self.label_1}</i></li>
                                    <li>sytematic <i>{self.label_2}</i>:        Includes positions that are present in all replicates of sample <i>{self.label_2}</i></li>
                                    <li><i>{self.label_1}</i> + <i>{self.label_2}</i>:       Positions from group (1) and (2) that are shared between the two</li>
                                    <li><i>{self.label_1}</i> excl.:              Positions from group (1) that are not present in group (2)</li>
                                    <li><i>{self.label_2}</i> excl.:              Positions from group (2) that are not present in group (1)</li>
                                </ol>
                            </p>
                            <div class="table-box">
                                {count_table}
                            </div>
                        </section>

                        <section>
                            <h2 class="collapsible-header">Genome map</h2>
                            <h3>Genome map displaying shared and exclusive positions between samples <i>{self.label_1}</i> and <i>{self.label_2}</i></h3>
                            <div class="plot-container">
                                {plot}
                            </div>
                            <p>
                                Genome map displaying the positions that are present in both samples, only in sample <i>{self.label_1}</i> and only in sample <i>{self.label_2}</i>.
                                Only chromosomes that contain extracted position(s) are shown.  
                            </p>
                        </section>
                    </body>
                    <footer></footer>
                    </html> 
            """
        with open(f"{self.out_dir}{self.label_1}_{self.label_2}_summary.html", "w") as out:
            out.write(template)

    def custom_sort_key(self, item):
        if item.isdigit():  # Check if the item is a digit
            return (int(item),)  # Convert to integer and sort numerically
        else:
            return (float('inf'), item)  # Place non-digits at the end


    def create_plot(self):
        def update_plot(fig, title: str|None = None, xlab: str|None = None, ylab: str|None = None, height: int = 500, width: int = 800):
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
        
        all_keys = self.ref_dict.keys()
        present_chr = set((all_keys & self.excl_in_1.chr.unique()) | (all_keys & self.excl_in_2.chr.unique()) | (all_keys & self.in_1_and_2.chr.unique()))
        present_chr = np.array(sorted(present_chr, key=self.custom_sort_key))
        chr_lens = np.array([len(self.ref_dict[x]) for x in present_chr])

        fig = update_plot(go.Figure(), height = 1000, width = 1200)
        fig.add_trace(go.Bar(x=present_chr, y=chr_lens, marker=dict(color="lightgrey"), name="Chromosomes", showlegend=False))
        fig.add_trace(go.Scatter(x=self.in_1_and_2.chr, y=self.in_1_and_2.site, mode='markers', marker=dict(symbol='line-ew', color="#1f77b4", size=18, line=dict(width=1.1, color="#1f77b4")), name=f"{self.label_1}+{self.label_2}", hovertemplate="Chr%{x}:%{y}"))
        fig.add_trace(go.Scatter(x=self.excl_in_1.chr, y=self.excl_in_1.site, mode='markers', marker=dict(symbol='line-ew', color="#ff7f0e", size=18, line=dict(width=1.1, color="#ff7f0e")), name=f"{self.label_1}", hovertemplate="Chr%{x}:%{y}"))
        fig.add_trace(go.Scatter(x=self.excl_in_2.chr, y=self.excl_in_2.site, mode='markers', marker=dict(symbol='line-ew', color="#2ca02c", size=18, line=dict(width=1.1, color="#2ca02c")), name=f"{self.label_2}", hovertemplate="Chr%{x}:%{y}"))
        fig.update_xaxes(fixedrange=True)

        return to_html(fig, include_plotlyjs=False)

    def get_file_paths(self) -> Tuple[str, str]:
        def create_list(paths: List[str]):
            path_list = "<ul>"
            for path in paths:
                path_list += f"<li>{path}</li>"
            path_list += "</ul>"
            return path_list
        
        return create_list(self.in_paths_1), create_list(self.in_paths_2)

    def create_count_table(self) -> str:
        count_table = pd.concat([self.common_pos_1.groupby("chr").count(), 
                         self.common_pos_2.groupby("chr").count(), 
                         self.in_1_and_2.groupby("chr").count(), 
                         self.excl_in_1.groupby("chr").count(), 
                         self.excl_in_2.groupby("chr").count()], axis=1).sort_values("chr", key=lambda x: x.map(self.custom_sort_key))
        count_table.columns = [f"systematic {self.label_1}", 
                               f"systematic {self.label_2}", 
                               f"{self.label_1} + {self.label_2}", 
                               f"{self.label_1} excl.", 
                               f"{self.label_2} excl."]
        count_table.index = count_table.index.rename("")
        count_table = count_table.T.fillna(0).astype(int)
        count_table["Genome-wide"] = count_table.sum(axis=1)
        return count_table.to_html()

    def main(self) -> None:
        self.extract_subsets_positions()
        self.write_bed_files()
        self.create_summary()

def setup_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="Neet - position extractor", description="Identify and map positions of interest between two samples.")
    parser.add_argument('-i', '--sample1', type=str, required=True,
                        help="""
                            Path to the input file(s). If replicates are available, specify paths comma-separated (<repl1.tsv>,<repl2.tsv>,...).
                            Must be of type tsv, as returned by the PileupExtractor.
                            """)
    parser.add_argument('-bn', '--basename', type=str, required=False, default="sample1",
                        help="""
                            Basename of the given sample. Used to create the pileup and extracted features files. Default: 'sample1'
                            """)
    parser.add_argument('-i2', '--sample2', type=str, required=False,
                        help="""
                            Path to the input file(s) from a second sample. If replicates are available, specify paths comma-separated (<repl1.tsv>,<repl2.tsv>,...).
                            Must be of type tsv, as returned by the PileupExtractor. 
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
    return parser



if __name__ == "__main__":
    parser = setup_parser()
    args = parser.parse_args()

    posextr = PositionExtractor(in_paths_1=args.sample1, in_paths_2=args.sample2, out_dir=args.ouput, ref_path=args.reference, label_1=args.basename, label_2=args.basename2)
    posextr.main()
    
    # posextr = PositionExtractor(in_paths_1="/home/vincent/masterthesis/data/nuc_cyt_3rep/processed/cytoplasm1_extracted.tsv,/home/vincent/masterthesis/data/nuc_cyt_3rep/processed/cytoplasm2_extracted.tsv,/home/vincent/masterthesis/data/nuc_cyt_3rep/processed/cytoplasm3_extracted.tsv",
    #                             in_paths_2="/home/vincent/masterthesis/data/nuc_cyt_3rep/processed/nucleus1_extracted.tsv,/home/vincent/masterthesis/data/nuc_cyt_3rep/processed/nucleus2_extracted.tsv,/home/vincent/masterthesis/data/nuc_cyt_3rep/processed/nucleus3_extracted.tsv",
    #                             out_dir="/home/vincent/masterthesis/data/nuc_cyt_3rep/processed/extracted_pos",
    #                             ref_path="/home/vincent/masterthesis/data/nuc_cyt_3rep/GRCh38_clean.genome.fa",
    #                             label_1="cytoplasm", label_2="nucleus")
    # posextr.main()
