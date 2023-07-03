import argparse
import plotly.graph_objects as go
from typing import List, Tuple, Dict
import sys
from helper_functions import float_between_zero_and_one

class CoverageTrackPlotter:
    """
    A class for creating bar chart visualizations based on provided data.

    The CoverageTrackPlotter class allows you to create bar chart visualizations based on data and settings
    provided as input. The resulting plot can be displayed or saved to a file.

    Attributes:
        in_path (str): The path to the input file.
        out_path (str): The path to save the output file.
        threshold (float): The threshold value for categorizing the data.
        chromosome (str): The chromosome to retrieve the sequence from.
        ref (str): The reference sequence for the chromosome.
        data_above_threshold (Dict[str, int]): Dictionary containing data points above the threshold.
        data_below_threshold (Dict[str, int]): Dictionary containing data points below the threshold.

    Methods:
        __init__(self, in_path: str, out_path: str, ref_path: str, threshold: float, chromosome: str = None) -> None:
            Initializes the CoverageTrackPlotter object with the provided parameters.

        get_ref_sequence(self, ref: str) -> str:
            Retrieves the sequence from a reference file based on the given chromosome.

        add_line(self, line: List[str], data_dict: Dict) -> Dict:
            Adds data from a line to a dictionary.

        set_up_data(self) -> Tuple[Dict[str, int]]:
            Sets up data from a file into separate dictionaries based on a threshold.

        get_custom_data(self, data_dict: dict) -> List[List[str|int]]:
            Retrieves custom data from a data dictionary to fill the hover templates.

        create_traces(self) -> List[go.Bar]:
            Creates traces for a bar chart visualization.

        create_plot(self) -> None:
            Creates a plot based on the provided data.

    """

    in_path: str
    out_path: str
    threshold: float
    chromosome: str
    ref: str
    data_above_threshold: Dict[str, int]
    data_below_threshold: Dict[str, int]

    def __init__(self, in_path: str, 
                 out_path: str, 
                 ref_path: str, 
                 threshold: float, 
                 chromosome: str = None) -> None:
        """
        Initializes the CoverageTrackPlotter object with the provided parameters.

        Args:
            in_path (str): The path to the input file.
            out_path (str): The path to save the output file.
            ref_path (str): The path to the reference file.
            threshold (float): The threshold value for categorizing the data.
            chromosome (str, optional): The chromosome to retrieve the sequence from. Defaults to None.

        Returns:
            None

        """
        self.in_path = in_path
        self.out_path = out_path
        self.chromosome = chromosome
        self.ref = self.get_ref_sequence(ref_path)
        self.threshold = threshold
        self.data_above_threshold, self.data_below_threshold = self.set_up_data()

    def get_ref_sequence(self, ref: str) -> str:
        """
        Retrieves the sequence from a reference file based on the given chromosome.

        Args:
            ref (str): The path to the reference file.

        Returns:
            str: The sequence corresponding to the specified chromosome.

        Raises:
            Exception: If the chromosome is not found in the reference file.
        """
        with open(ref, "r") as r:
            for line in r:
                if (line.startswith(">")) & (line[1:].strip() == self.chromosome):
                    return next(r).strip()
            raise Exception(f"Sequence not found. Chromosome '{chr}' not found in '{ref}'.")

    def add_line(self, line: List[str], data_dict: Dict):
        """
        Adds data from a line to a dictionary.

        The function takes a line, which is a list of strings representing various data points,
        and a data dictionary, which is a dictionary to store the data points. The function appends
        the data points from the line to the corresponding lists in the data dictionary.

        Args:
            line (list[str]): A list of strings representing the data points to be added.
            data_dict (dict): A dictionary to store the data points.

        Returns:
            dict: The updated data dictionary.
        """
        data_dict["site"].append(int(line[1]))
        data_dict["reads"].append(int(line[2]))
        data_dict["a"].append(int(line[5]))
        data_dict["c"].append(int(line[6]))
        data_dict["g"].append(int(line[7]))
        data_dict["t"].append(int(line[8]))
        data_dict["del"].append(int(line[9]))
        data_dict["ins"].append(int(line[10]))
        return data_dict


    def set_up_data(self) -> Tuple[Dict[str, int]]:
        """
        Sets up data from a file into separate dictionaries based on a threshold.

        The function reads data from the specified input file and separates it into two dictionaries:
        'above_thr' and 'below_thr'. The data is categorized based on a threshold value. If the 'perc_mis'
        value in a line is greater than or equal to the threshold, the data points are added to the 'above_thr'
        dictionary. Otherwise, the data points are added to the 'below_thr' dictionary.

        Args:
            infile (str): The path to the input file.

        Returns:
            Tuple[Dict[str, int]]: A tuple containing the 'above_thr' and 'below_thr' dictionaries.

        """
        above_thr = {"site": [], "reads": [], "a": [], "c": [], "g": [], "t": [], "del": [], "ins": []}
        below_thr = {"site": [], "reads": [], "a": [], "c": [], "g": [], "t": [], "del": [], "ins": []}
        
        with open(self.in_path, "r") as file:
            next(file)
            for line in file:
                line = line.strip("\n").split("\t")
                if line[0] == self.chromosome:
                    perc_mis = float(line[17])

                    if perc_mis >= self.threshold:
                        above_thr = self.add_line(line, above_thr)
                    else:
                        below_thr = self.add_line(line, below_thr)
        
        return above_thr, below_thr

    def get_custom_data(self, data_dict: dict) -> List[List[str|int]]:
        """
        Retrieves custom data from a data dictionary to fill the hover templates.

        The function takes a data dictionary as input and retrieves custom data from it.
        The custom data consists of a list of lists, where each inner list represents a data point
        and contains the following elements: Chromosome, Reads, A, C, G, T, Deletion (Del), Insertion (Ins), Reference (Ref)

        Args:
            data_dict (dict): A dictionary containing the data points.

        Returns:
            List[List[str|int]]: A list of lists representing the custom data.
        """
        custom_data = [[self.chromosome, 
                        data_dict["reads"][i],
                        data_dict["a"][i],
                        data_dict["c"][i],
                        data_dict["g"][i],
                        data_dict["t"][i], 
                        data_dict["del"][i],
                        data_dict["ins"][i],
                        self.ref[i]]
                        for i in range(len(data_dict["site"]))]    
        return custom_data

    def create_traces(self) -> List[go.Bar]:
        """
        Creates traces for a bar chart visualization.

        The function creates traces for a bar chart visualization based on the data provided.
        The resulting traces represent different data points and are returned as a list.

        Returns:
            List[go.Bar]: A list of traces for a bar chart visualization.

        """
        hover_template_top = "%{customdata[0]}:%{x}<br>" + \
            "Total count: %{customdata[1]}<br>" + \
            "----------"

        hover_template_bottom = "----------<br>" + \
            "DEL: %{customdata[6]}<br>" + \
            "INS: %{customdata[7]}<br>" + \
            "----------<br>" + \
            "Ref: %{customdata[8]}<extra></extra>"

        hover_template = hover_template_top + \
            "<br>A: %{customdata[2]}<br>" + \
            "C: %{customdata[3]}<br>" + \
            "G: %{customdata[4]}<br>" + \
            "T: %{customdata[5]}<br>" + \
            hover_template_bottom

        all_traces = []

        hidden_trace_top = go.Bar(x = self.data_above_threshold["site"], 
                                  y = [0 for _ in self.data_above_threshold["site"]],
                                  name = "", 
                                  customdata = self.get_custom_data(self.data_above_threshold),
                                  hovertemplate = hover_template_top, 
                                  marker_color = "white")
        all_traces.append(hidden_trace_top)

        colors = {"a": "green", "c": "blue", "g": "orange", "t": "red"}
        bar_bases = ["a", "c", "g", "t"]
        for base in bar_bases:
            trace = go.Bar(x = self.data_above_threshold["site"], 
                           y = self.data_above_threshold[base],
                           name = base.upper(),
                           marker_color = colors[base])
            all_traces.append(trace)

        remaining_trace = go.Bar(x = self.data_above_threshold["site"], 
                                 y = self.data_above_threshold["del"],
                                 name = "", 
                                 marker_color = "grey", 
                                 hoverinfo = "skip")

        hidden_trace_bottom = go.Bar(x = self.data_above_threshold["site"], 
                                     y = [0 for _ in self.data_above_threshold["site"]],
                                     customdata = self.get_custom_data(self.data_above_threshold),
                                     hovertemplate = hover_template_bottom, 
                                     marker_color = "white")

        nonerror_hidden_trace = go.Bar(x = self.data_below_threshold["site"], 
                                       y = [0 for _ in self.data_below_threshold["site"]],
                                       customdata = self.get_custom_data(self.data_below_threshold),
                                       hovertemplate = hover_template, 
                                       marker_color = "white")

        nonerror_trace = go.Bar(x = self.data_below_threshold["site"], 
                                y = self.data_below_threshold["reads"], 
                                name = "Total count", 
                                marker_color = "grey", 
                                hoverinfo = "skip")
        
        all_traces += [remaining_trace, hidden_trace_bottom, nonerror_hidden_trace, nonerror_trace]
        return all_traces

    def create_plot(self) -> None:
        """
        Creates a plot based on the provided data.

        The function creates a plot based on the data and settings provided. The resulting plot
        is displayed or saved to a file, depending on the configuration.

        Returns:
            None
        """
        all_traces = self.create_traces()

        fig = go.Figure()
        fig.update_xaxes(range=[1, len(self.ref)])
        fig.add_traces(all_traces)

        fig.update_layout(
            template = "simple_white",
            title="Coverage track",
            xaxis_title='Site',
            yaxis_title='Count',
            barmode='stack',
            hovermode="x unified",
            showlegend=False,
            yaxis = dict(fixedrange=True)
        )

        if self.out_path:
            fig.write_html(self.out_path)
        else:
            fig.show()


if __name__=="__main__":
    parser = argparse.ArgumentParser(prog="CoverageTrackPlotter",
                                     description="Provides different plotting functionalities.")
    subparsers = parser.add_subparsers(dest="subcommand")
    coverage_track_parser = subparsers.add_parser("coverage_track")
    coverage_track_parser.add_argument("-i", "--input", type=str, required=True,
                                       help="Path to the input file")
    coverage_track_parser.add_argument("-o", "--output", type=str, required=False,
                                       help="Path to the output file")
    coverage_track_parser.add_argument("-r", "--reference", type=str, required=True,
                                       help="Path to the reference fasta file")
    coverage_track_parser.add_argument('-t', '--error_threshold', type=float_between_zero_and_one, required=True,
                        help="Percentage of error threshold for coloring positions [0,1]")
    coverage_track_parser.add_argument('-c', '--chromosome', type=str, required=True,
                        help="Chromosome to be plotted")
    
    args = parser.parse_args()

    if args.subcommand == "coverage_track":
        plotter = CoverageTrackPlotter(in_path=args.input,
                          ref_path=args.reference,
                          threshold=args.error_threshold,
                          chromosome=args.chromosome,
                          out_path=args.output)
        plotter.create_plot()
    else:
        parser.print_help()