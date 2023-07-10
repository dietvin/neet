import os
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from typing import Set, Tuple, List, Dict
from helper_functions import check_get_in_path

class CompositionPlotter:
    """A class for creating composition plots based on TSV files and a BED file."""

    positions: Set[Tuple[str, int]]
    filenames: List[str]
    data: List[Dict[str, str|float]]

    def __init__(self, tsv_paths: str, bed_path: str) -> None:
        """
        Initializes a CompositionPlotter object.

        Args:
            tsv_paths (str): Comma-separated paths to TSV files or a directory path containing TSV files.
            bed_path (str): Path to the BED file.

        Returns:
            None
        """
        self.filenames = self.get_filenames(tsv_paths)
        self.positions = self.get_bed_positions(bed_path)
        self.data = self.get_info()


    def get_filenames(self, tsv_paths: str) -> List[str]:
        """
        Retrieves the list of TSV file names.

        Args:
            tsv_paths (str): Comma-separated paths to TSV files or a directory path containing TSV files.

        Returns:
            List[str]: List of TSV file names.
        """
        def is_list_of_paths(tsv_paths: str) -> bool:
            paths = tsv_paths.split(",")
            for path in paths:
                if not os.path.isfile(path.strip()):
                    return False
            return True
        
        if is_list_of_paths(tsv_paths):
            return tsv_paths.split(",")
        else:
            file_names = []
            for file_name in os.listdir(tsv_paths):
                if file_name.endswith(".tsv") and os.path.isfile(os.path.join(tsv_paths, file_name)):
                    file_names.append(os.path.join(tsv_paths, file_name))
            return file_names


    def get_bed_positions(self, bed_file: str) -> Set[Tuple[str, int]]:
        """
        Retrieves positions from the BED file.

        Args:
            bed_file (str): Path to the BED file.

        Returns:
            Set[Tuple[str, int]]: Set of positions (chromosome, position).
        """
        positions = set()
        with open(bed_file, 'r') as bed:
            for line in bed:
                fields = line.strip().split('\t')
                chromosome = fields[0]
                position = int(fields[1])
                positions.add((chromosome, position))
        return positions

    def get_info(self) -> List[Dict[str, str|float]]:
        """
        Retrieves information from TSV files.

        Returns:
            List[Dict[str, str|float]]: List of dictionaries containing information for each position.
        """
        counts = []
        for position in self.positions:
            pos_count = {"pos": f"{position[0]}:{position[1]}",
                        "name": [],
                        "a": [],
                        "c": [],
                        "g": [],
                        "t": []}

            for sample in self.filenames:
                basename = os.path.splitext(os.path.basename(sample))[0]
                pos_count["name"].append(basename)

                vals = self.extract_rows(sample, position)
                pos_count["a"].append(vals[0])
                pos_count["c"].append(vals[1])
                pos_count["g"].append(vals[2])
                pos_count["t"].append(vals[3])
            counts.append(pos_count)
        return counts

    def extract_rows(self, tsv_file: str, ref_position: str) -> Set[float]:
        """
        Extracts rows from the TSV file based on the reference position.

        Args:
            tsv_file (str): Path to the TSV file.
            ref_position (str): Reference position (chromosome, position).

        Returns:
            Set[float]: Set of values extracted from the TSV file.
        """
        with open(tsv_file, 'r') as tsv:
            next(tsv)
            for line in tsv:
                fields = line.strip().split('\t')
                chromosome = fields[0]
                site = int(fields[1])
                position = (chromosome, site)
                if position == ref_position:
                    return float(fields[11]), float(fields[12]), float(fields[13]), float(fields[14])
        return 0, 0, 0, 0

    def create_plots(self) -> None:
        """
        Creates and displays the composition plots using Plotly.

        Returns:
            None
        """
        fig = make_subplots(rows=len(self.data), cols=1, shared_xaxes=True, subplot_titles=[d['pos'] for i, d in enumerate(self.data)])

        legend = True
        for i, data in enumerate(self.data):
            sub_fig = go.Figure()
            sub_fig.add_trace(go.Bar(x=data["name"], y=data["a"], name="A", marker_color="green", legendgroup="group_1", showlegend = legend))
            sub_fig.add_trace(go.Bar(x=data["name"], y=data["c"], name="C", marker_color="blue", legendgroup="group_1", showlegend = legend))
            sub_fig.add_trace(go.Bar(x=data["name"], y=data["g"], name="G", marker_color="orange", legendgroup="group_1", showlegend = legend))
            sub_fig.add_trace(go.Bar(x=data["name"], y=data["t"], name="T", marker_color="red", legendgroup="group_1", showlegend = legend))
            sub_fig.update_layout(
                title=f"{data['pos']}")

            for trace in sub_fig.data:
                fig.add_trace(trace, row=i + 1, col=1)

            fig.update_layout(
                template="simple_white",
                barmode='stack',
                title=f"Relative proportions of bases between sampes",
                yaxis=dict(range=[0, 1], fixedrange=True),
                legend=dict(y=0.5))
            
            if i == 0:
                legend = False

        fig.show()

if __name__ == "__main__":
    file_mod = "/home/vincent/masterthesis/data/nanocompore_data/processed/Oligo_1_extracted.tsv"
    file_unm = "/home/vincent/masterthesis/data/nanocompore_data/processed/Oligo_control_extracted.tsv"
    file_x = "/home/vincent/masterthesis/data/nanocompore_data/processed/Oligo_x_extracted.tsv"
    bed = "/home/vincent/masterthesis/data/nanocompore_data/processed/composition_plot_test.bed"

    plotter = CompositionPlotter(file_mod+","+file_unm+","+file_x, bed)
    plotter.create_plots()
