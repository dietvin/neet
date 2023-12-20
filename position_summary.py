import helper_functions as hs
from typing import Tuple, List, Dict
import argparse, os, warnings, datetime
import plotly.graph_objects as go
from plotly.subplots import make_subplots

class PositionSummary:
    positions: List[Tuple[str, int]]
    nb_size: int

    basename_a: str
    basename_b: str
    
    paths_a: List[str]
    paths_b: List[str]

    bed_path: str

    out_dir: str
    """
    Data structure:
    data_sample_a
        position1{[
            replicate1[
                [chr, site, ...] (pos-nb_size)
                ...
                [chr, site, ...] (pos)
                ...
                [chr, site, ...] (pos+nb_size)
                ],
            replicate2[
                ...
                ],
            ...
            replicaten[
                ...
                ]
            ]}
        position2{
        ...
        }
        ...
        positionn{
        ...
        }
    """
    data_sample_a: Dict[Tuple[chr, int], List[List[List[str]]]]
    data_sample_b: Dict[Tuple[chr, int], List[List[List[str]]]]

    colors = {"A": "#2ca02c", "C": "#1f77b4", "G": "#ff7f0e", "U": "#d62728", 
              "match": "#9467bd", "mis": "#8c564b", "del": "#e377c2", 
              "ins": "#7f7f7f", "ref_skip": "#7f7f7f"}

    def __init__(self, paths_a: str, paths_b: str, basename_a: str, basename_b: str, bed_path: str, out_path: str, nb_size: int) -> None:
        
        self.basename_a = basename_a
        self.basename_b = basename_b

        self.check_out(out_path)
        if not out_path.endswith("/"):
            out_path += "/"
        self.out_dir = out_path

        self.check_path(bed_path, extension=".bed")
        self.get_positions(bed_path)
        self.bed_path = bed_path
        self.nb_size = nb_size

        paths_a = self.process_in(paths_a)
        paths_b = self.process_in(paths_b)

        self.paths_a = paths_a
        self.paths_b = paths_b

        self.get_data(paths_a, paths_b)

    ##########################################################################################################
    #                                  Methods called during initialization                                  #
    ##########################################################################################################
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

    def check_out(self, out: str) -> None: 
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

    def get_positions(self, bed_path: str) -> None:
        hs.print_update(f"Extracting positions from {bed_path}... ", line_break=False)
        with open(bed_path, "r") as bed_file:
            positions = []
            for line in bed_file:
                line = line.strip().split("\t")
                chrom, site = line[0], int(line[2])
                positions.append((chrom, site))
            hs.print_update(f"Found {len(positions)} sites.", with_time=False)
            self.positions = positions

    def get_data(self, paths_a: List[str], paths_b: List[str]) -> None:

        def process_sample(paths: List[str]) -> List[List[List[str]]]:
            target_sites = set([(position[0], position[1]+i) for i in range(-self.nb_size, self.nb_size+1) for position in self.positions])

            rep_data = []
            for path in paths:
                with open(path, "r") as file:
                    data_collection = {}
                    hs.print_update(f"Processing file {path}")
                    next(file)
                    for line in file:
                        line = line.strip().split("\t")
                        current_pos = line[0], int(line[1])
                        if (current_pos[0], current_pos[1]) in target_sites:
                            data_collection[current_pos] = line
                    rep_data.append(data_collection)
                                
            position_data = {}
            for position in self.positions:
                data_ordered = []
                for rep in rep_data:
                    data = [None] * (2 * self.nb_size + 1)

                    sites = [position[1]+i for i in range(-self.nb_size, self.nb_size+1)]
                    neighbour_positions = [(position[0], site) for site in sites]
                    for neighbour_position in neighbour_positions:
                        if neighbour_position in rep.keys():
                            data_index = neighbour_position[1] - position[1] + self.nb_size
                            data[data_index] = rep[neighbour_position]
                    data_ordered.append(data)

                position_data[position] = data_ordered
            
            return position_data
        
        self.data_sample_a = process_sample(paths_a)
        self.data_sample_b = process_sample(paths_b)

    ##########################################################################################################
    #                                    Methods called for plot creation                                    #
    ##########################################################################################################
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


    def create_bar_trace_base_composition(self, pos_data: List[str], x: int = 0, center: bool = True, show_legend: bool = False) -> List[go.Bar]:
        a, c, g, u = int(pos_data[5]), int(pos_data[6]), int(pos_data[7]), int(pos_data[8])
        cvg = a+c+g+u
        a_rel, c_rel, g_rel, u_rel = a/cvg, c/cvg, g/cvg, u/cvg
        a_bar = go.Bar(x = [x], y = [a_rel], name = "A", 
                    marker = dict(color = self.colors["A"], line_color="black", line_width=1.5) if center else dict(color = self.colors["A"]), 
                    legendgroup = 1, 
                    opacity = 1 if center else 0.75, 
                    showlegend = show_legend, 
                    customdata = [[[a], [round(a_rel*100, 2)], [cvg]]], 
                    hovertemplate = "n=%{customdata[0]} (%{customdata[1]}%)<br>Cvg=%{customdata[2]}",
                    width=0.75)
        c_bar = go.Bar(x = [x], y = [c_rel], name = "C", 
                    marker = dict(color = self.colors["C"], line_color="black", line_width=1.5) if center else dict(color = self.colors["C"]), 
                    legendgroup = 2, 
                    opacity = 1 if center else 0.75, 
                    showlegend = show_legend, 
                    customdata = [[[c], [round(c_rel*100, 2)], [cvg]]], 
                    hovertemplate = "n=%{customdata[0]} (%{customdata[1]}%)<br>Cvg=%{customdata[2]}",
                    width=0.75)
        g_bar = go.Bar(x = [x], y = [g_rel], name = "G", 
                    marker = dict(color = self.colors["G"], line_color="black", line_width=1.5) if center else dict(color = self.colors["G"]),
                    legendgroup = 3, 
                    opacity = 1 if center else 0.75, 
                    showlegend = show_legend, 
                    customdata = [[[g], [round(g_rel*100, 2)], [cvg]]], 
                    hovertemplate = "n=%{customdata[0]} (%{customdata[1]}%)<br>Cvg=%{customdata[2]}",
                    width=0.75)
        u_bar = go.Bar(x = [x], y = [u_rel], name = "U", 
                    marker = dict(color = self.colors["U"], line_color="black", line_width=1.5) if center else dict(color = self.colors["U"]), 
                    legendgroup = 4, 
                    opacity = 1 if center else 0.75, 
                    showlegend = show_legend, 
                    customdata = [[[u], [round(u_rel*100, 2)], [cvg]]], 
                    hovertemplate = "n=%{customdata[0]} (%{customdata[1]}%)<br>Cvg=%{customdata[2]}",
                    width=0.75)
        return [a_bar, c_bar, g_bar, u_bar]

    def create_bar_trace_error_rates(self, pos_data: List[str], x: int = 0, center: bool = True, show_legend: bool = False) -> List[go.Bar]:
        base_idx_rel = {"A": 12, "C": 13, "G": 14, "T": 15}
        base_idx_abs = {"A": 5, "C": 6, "G": 7, "T": 8}

        n_match = int(pos_data[base_idx_abs[pos_data[3]]])
        n_mismatch = int(pos_data[5]) + int(pos_data[6]) + int(pos_data[7]) + int(pos_data[8])
        n_deletion = int(pos_data[9])
        n_refskip = int(pos_data[11])

        match_rate = float(pos_data[base_idx_rel[pos_data[3]]])
        mismatch_rate = float(pos_data[19])
        deletion_rate = float(pos_data[16])
        refskip_rate = float(pos_data[18])

        cvg = int(pos_data[2])

        match_bar = go.Bar(x = [x], y = [match_rate], name = "Match", 
                        marker = dict(color = self.colors["match"], line_color="black", line_width=1.5) if center else dict(color = self.colors["match"]), 
                        legendgroup = 1, 
                        opacity = 1 if center else 0.75, 
                        showlegend = show_legend, 
                        customdata = [[[n_match], [round(match_rate*100, 2)], [cvg]]], 
                        hovertemplate = "n=%{customdata[0]} (%{customdata[1]}%)<br>Cvg=%{customdata[2]}",
                        width=0.75)
        mismatch_bar = go.Bar(x = [x], y = [mismatch_rate], name = "Mismatch", 
                            marker = dict(color = self.colors["mis"], line_color="black", line_width=1.5) if center else dict(color = self.colors["mis"]), 
                            legendgroup = 2, 
                            opacity = 1 if center else 0.75, 
                            showlegend = show_legend, 
                            customdata = [[[n_mismatch], [round(mismatch_rate*100, 2)], [cvg]]], 
                            hovertemplate = "n=%{customdata[0]} (%{customdata[1]}%)<br>Cvg=%{customdata[2]}",
                            width=0.75)
        deletion_bar = go.Bar(x = [x], y = [deletion_rate], name = "Deletion", 
                            marker = dict(color = self.colors["del"], line_color="black", line_width=1.5) if center else dict(color = self.colors["del"]),
                            legendgroup = 3, 
                            opacity = 1 if center else 0.75, 
                            showlegend = show_legend, 
                            customdata = [[[n_deletion], [round(deletion_rate*100, 2)], [cvg]]], 
                            hovertemplate = "n=%{customdata[0]} (%{customdata[1]}%)<br>Cvg=%{customdata[2]}",
                            width=0.75)
        refskip_bar = go.Bar(x = [x], y = [refskip_rate], name = "Reference skip", 
                            marker = dict(color = self.colors["ref_skip"], line_color="black", line_width=1.5) if center else dict(color = self.colors["ref_skip"]), 
                            legendgroup = 4, 
                            opacity = 1 if center else 0.75, 
                            showlegend = show_legend, 
                            customdata = [[[n_refskip], [round(refskip_rate*100, 2)], [cvg]]], 
                            hovertemplate = "n=%{customdata[0]} (%{customdata[1]}%)<br>Cvg=%{customdata[2]}",
                            width=0.75)
        return [match_bar, mismatch_bar, deletion_bar, refskip_bar]

    def create_bar_trace_insertion_rate(self, pos_data: List[str], x: int = 0, center: bool = True, show_legend: bool = False) -> List[go.Bar]:
        n_ins = int(pos_data[10])
        ins_rate = float(pos_data[17])

        cvg = int(pos_data[2])

        ins_bar = go.Bar(x = [x], y = [ins_rate], name = "Insertion", 
                        marker = dict(color = self.colors["ins"], line_color="black", line_width=1.5) if center else dict(color = self.colors["ins"]), 
                        legendgroup = 1, 
                        opacity = 1 if center else 0.75, 
                        showlegend = show_legend, 
                        customdata = [[[n_ins], [round(ins_rate*100, 2)], [cvg]]], 
                        hovertemplate = "n=%{customdata[0]} (%{customdata[1]}%)<br>Cvg=%{customdata[2]}",
                        width=0.75)
        return [ins_bar]

    def get_rep_bar_traces(self, data: List[List[str]], plot_type: str, showlegend: bool = False) -> List[go.Bar]:
        x_vals = range(-self.nb_size,self.nb_size+1)
        is_center_pos = [False]*self.nb_size + [True] + [False]*self.nb_size
        show_legend = is_center_pos if showlegend else [showlegend]*(2*self.nb_size+1)

        if plot_type == "base_composition":
            create_bar_func = self.create_bar_trace_base_composition
        elif plot_type == "error_rates":
            create_bar_func = self.create_bar_trace_error_rates
        elif plot_type == "insertion_rate":
            create_bar_func = self.create_bar_trace_insertion_rate
        else:
            raise Exception(f"Unknown plot type given: {plot_type}. Must be one of the following: base_composition, error_rates, insertion_rate")
        
        bar_traces = [create_bar_func(d, x, c, l) for d, x, c, l in zip(data, x_vals, is_center_pos, show_legend)]
        return [item for t in bar_traces for item in t]

    def create_position_plots(self, data_a: List[List[List[str]]], data_b: List[List[List[str]]]) -> Tuple[str]:
        n_rep_a = len(data_a)
        n_rep_b = len(data_b)
        n_samples = n_rep_a+n_rep_b

        figs = []

        for plot_type, width in zip(["base_composition", "error_rates", "insertion_rate"], [1200, 700, 500]):
            fig = make_subplots(cols=1, rows=n_samples, shared_xaxes=True, vertical_spacing=0.02)
            fig = self.update_plot(fig, height=200*(n_samples), width=width)

            current_row = 1
            for rep in data_a:
                bars = self.get_rep_bar_traces(rep, plot_type=plot_type, showlegend=True) if current_row==1 else self.get_rep_bar_traces(rep, plot_type=plot_type)
                fig.add_traces(bars, rows=[current_row for _ in range(len(bars))], cols=[1 for _ in range(len(bars))])
                current_row += 1
            for rep in data_b:
                bars = self.get_rep_bar_traces(rep, plot_type=plot_type)
                fig.add_traces(bars, rows=[current_row for _ in range(len(bars))], cols=[1 for _ in range(len(bars))])
                current_row += 1
            fig.update_layout(barmode="stack")

            y_labels = [f"{self.basename_a} {i}" for i in range(1,n_rep_a+1)] + [f"{self.basename_b} {i}" for i in range(1,n_rep_b+1)]
            for i in range(n_samples):
                fig.update_yaxes(title=y_labels[i], row=i+1, col=1)
                
            for i in range(1, n_samples):
                fig.update_xaxes(tickvals=[], ticktext=[], row=i, col=1)
            fig.update_xaxes(tickvals=list(range(-self.nb_size,self.nb_size+1)), ticktext=[d[3] for d in data_a[0]], row=n_samples, col=1)

            figs.append(fig.to_html(include_plotlyjs=False))
        return tuple(figs)
    
    def create_html_section(self, position: Tuple[str, int], plots: Tuple[str, str, str]) -> str:
        chrom, site = position[0], position[1]
    
        plot_base_comp = plots[0]
        plot_err_rates = plots[1]
        plot_ins_rate = plots[2]

        collapsible_section = f"""
            <section>
                <button class="collapsible">{chrom}:{site}</button>

                <div class="collapsible-content">
                <h2 class="hiddentitle" id="{chrom}_{site}"></h2>

                    <h3>Base compositions</h3>
                    <div class="plot-container">
                        {plot_base_comp}
                    </div>

                    <h3>Error rates</h3>
                    <div class="plot-container">
                        <div id="left">
                            {plot_err_rates}
                        </div>
                        <div id="right">
                            {plot_ins_rate}
                        </div>
                    </div>
                </div>
            </section>
        """
        return collapsible_section

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
        
        return create_list(self.paths_a), create_list(self.paths_b)

    def main(self) -> str:
        collapsible_sections = ""
        for (position_a, data_a), (position_b, data_b) in zip(self.data_sample_a.items(), self.data_sample_b.items()):
            plots = self.create_position_plots(data_a, data_b)
            collapsible_sections += self.create_html_section(position_a, plots) 
        
        time = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
        files_a, files_b = self.get_file_paths()
        # collapsible sections adapted from: https://github.com/wdecoster/NanoPlot/blob/master/nanoplot/report.py
        template = f"""
            <!DOCTYPE html>
            <html lang="en">
            <head>
                <meta charset="UTF-8">
                <meta http-equiv="X-UA-Compatible" content="IE=edge">
                <meta name="viewport" content="width=device-width, initial-scale=1.0">
                <title>Neet - Position extractor summary</title>
                <link rel="stylesheet" type="text/css" href="/home/vincent/projects/neet_project/neet/style.css">
            </head>

            <body>
                <script type="text/javascript" src="/home/vincent/projects/neet_project/neet/plotly_js.js"></script>

                <header>
                    <h1>Position extractor summary</h1>
                    <p>Produced by <a href="https://github.com/dietvin/neet">Neet</a> on <b>{time}</b></p>
                </header>
            
                <section>
                    <p class="intro-text">
                        This summary file contains an overview of {len(self.positions)} positions found in <i>{self.bed_path}</i>. 
                    </p>
                    <p class="intro-text">Files provided for sample <i>{self.basename_a}</i>:</p>
                    {files_a}
                    <p class="intro-text">Files provided for sample <i>{self.basename_b}</i>:</p>
                    {files_b}
                </section>

                {collapsible_sections}

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
        outfile = f"{self.out_dir}{self.basename_a}_{self.basename_b}_summary.html"
        with open(outfile, "w") as out:
            hs.print_update(f"Writing summary file to: {outfile}")
            out.write(template)







def setup_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="Neet - Position Summary", description="Create summaries of positions in BED file.")
    parser.add_argument('-i', '--sample1', type=str, required=True,
                        help="""
                            Path to the input file(s). If replicates are available, specify paths comma-separated (<repl1.tsv>,<repl2.tsv>,...).
                            Must be of type tsv, as returned by the PileupExtractor.
                            """)
    parser.add_argument('-bn', '--basename1', type=str, default="sample1",
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
    parser.add_argument('-b', '--bed_file', type=str, required=True,
                        help="""
                            Path to the input BED file.
                            """)
    parser.add_argument('-o', '--output_dir', type=str, required=True,
                        help="""
                            Path to output directory. The output file is created using the basenames from the input samples.
                            """)
    parser.add_argument('-n', '--n_surrounding', type=hs.positive_int, required=False, default=2,
                        help="Number of neighbouring positions to consider when plotting a position X. In the created plots sites X-n, ..., X, ..., X+n are shown. Default: 2")
    return parser

if __name__=="__main__":
    parser = setup_parser()
    args = parser.parse_args()

    position_summary = PositionSummary(paths_a = args.sample1,
                                       paths_b = args.sample2,
                                       basename_a = args.basename1,
                                       basename_b = args.basename2,
                                       bed_path = args.bed_file,
                                       out_path = args.output_dir,
                                       nb_size = args.n_surrounding)
    position_summary.main()   
