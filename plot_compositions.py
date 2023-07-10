import os
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from typing import Set, Tuple, List, Dict
from helper_functions import check_get_in_path

class CompositionPlotter:
    positions: Set[Tuple[str, int]]
    filenames: List[str]
    data: List[Dict[str, str|float]]

    def __init__(self, tsv_paths: str, bed_path: str) -> None:
        
        self.filenames = tsv_paths.split(",")
        self.positions = self.get_bed_positions(bed_path)
        self.data = self.get_info()

    def get_bed_positions(self, bed_file: str) -> Set[Tuple[str, int]]:
        positions = set()
        with open(bed_file, 'r') as bed:
            for line in bed:
                fields = line.strip().split('\t')
                chromosome = fields[0]
                position = int(fields[1])
                positions.add((chromosome, position))
        return positions

    def get_info(self) -> List[Dict[str, str|float]]:
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
