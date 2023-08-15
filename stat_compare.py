from typing import Dict, List, Tuple
import warnings, sys, os, datetime

import pandas as pd
import matplotlib as plt
import matplotlib.colors as mcolors

from scipy.stats import wilcoxon

class StatComparer:
    path_1: str
    path_2: str
    outpath: str
    alpha: float
    write_tsv: bool
    data1_tmp: Dict[str, List[float|str]]
    data2_tmp: Dict[str, List[float|str]]
    data1: pd.DataFrame
    data2: pd.DataFrame
    results: pd.DataFrame
    stats: Dict[str, int]

    IDX = {"n_a_rel": 7, 
           "n_c_rel": 8, 
           "n_g_rel": 9, 
           "n_t_rel": 10, 
           "n_del_rel": 11,
           "n_ins_rel": 12, 
           "n_ref_skip_rel": 13, 
           "perc_mismatch": 14, 
           "q_mean": 16}

    def __init__(self, path_1: str, path_2: str, outpath: str, alpha: float = 0.01, write_tsv: bool = True) -> None:
        """
        Initialize the StatComparer instance.

        Args:
            path_1 (str): Path to the first input file.
            path_2 (str): Path to the second input file.
            outpath (str): Path to the output file.
            alpha (float, optional): Significance level for statistical tests. Defaults to 0.01.
        """
        self.process_paths(path_1, path_2, outpath)
        self.alpha = alpha
        self.write_tsv = write_tsv
        self.data1_tmp = {"n_a_rel": [],
                      "n_c_rel": [],
                      "n_g_rel": [],
                      "n_t_rel": [],
                      "n_del_rel": [],
                      "n_ins_rel": [],
                      "n_ref_skip_rel": [],
                      "perc_mismatch": [],
                      "q_mean": [],
                      "chrom": [],
                      "is_mismatch": []}
        self.data2_tmp = {"n_a_rel": [],
                      "n_c_rel": [],
                      "n_g_rel": [],
                      "n_t_rel": [],
                      "n_del_rel": [],
                      "n_ins_rel": [],
                      "n_ref_skip_rel": [],
                      "perc_mismatch": [],
                      "q_mean": [],
                      "chrom": [],
                      "is_mismatch": []}
        self.stats = {}

    def process_paths(self, in1: str, in2: str, out: str) -> None:
        """
        Process input and output paths, ensuring their validity.

        Args:
            in1 (str): Path to the first input file.
            in2 (str): Path to the second input file.
            out (str): Path to the output file.
        """
        self.check_path(in1, [".tsv"])
        self.path_1 = in1
        self.check_path(in2, [".tsv"])
        self.path_2 = in2
        self.process_outpath(out, in1, in2)

    def check_path(self, path: str, extensions: List[str]) -> None:
        """
        Check if a given file path exists and has a valid extension.

        Args:
            path (str): File path to be checked.
            extensions (List[str]): List of valid file extensions.

        Returns:
            None
        """
        if not os.path.exists(path):
            raise FileNotFoundError(f"Input file not found. File '{path}' does not exist.")
        file_type = os.path.splitext(path)[1]
        if not file_type in extensions:
            warnings.warn(f"Found file extension {file_type}. Expected file extension to be one of: {extensions}. If this is deliberate, ignore warning.", Warning)

    def process_outpath(self, out: str, in1: str, in2: str) -> None:
        """
        Process the output path for the result file.

        Args:
            out (str): Output file path.
            in1 (str): Path to the first input file.
            in2 (str): Path to the second input file.
        """
        if os.path.isdir(out): 
            if not os.path.exists(out):
                raise FileNotFoundError(f"Directory not found. Output directory '{out}' does not exist.")
            if not out.endswith("/"):
                out += "/"

            basename1 = os.path.splitext(os.path.basename(in1))[0]
            basename2 = os.path.splitext(os.path.basename(in2))[0]

            self.outpath = f"{out}{basename1}_{basename2}_comp.tsv"

        else: 
            dirname = os.path.dirname(out)
            if not os.path.exists(dirname):
                raise FileNotFoundError(f"Path to output file not found. '{dirname}' does not exist.")
            file_extension = os.path.splitext(out)[1]
            if file_extension != ".tsv":
                warnings.warn(f"Given output file has extension '{file_extension}'. Note that the output file will be of type '.tsv'.")

            self.outpath = out

    def find_common_pos(self) -> None:
        """
        Find common positions between the two input datasets and collect data for comparison.
        """
        with open(self.path_1, 'r') as file1, open(self.path_2, 'r') as file2:
            # skip header
            next(file1)
            next(file2)

            line1 = file1.readline()
            line2 = file2.readline()
            
            n_1 = 0
            n_2 = 0
            n_common = 0

            while line1 and line2:
                chr1, site1, _, ref_base1, maj_base1, *data1 = line1.split('\t')
                chr2, site2, _, ref_base2, maj_base2, *data2 = line2.split('\t')

                if chr1 == chr2 and site1 == site2:
                    self.add_data(data1, data2, chr1, ref_base1!=maj_base1, ref_base2!=maj_base2)
                    line1 = file1.readline()
                    line2 = file2.readline()
                    n_common += 1
                    n_1 += 1
                    n_2 += 1
                elif (chr1, site1) < (chr2, site2):
                    line1 = file1.readline()
                    n_1 += 1
                else:
                    line2 = file2.readline()
                    n_2 += 1
            self.stats = {"n_1": n_1, "n_2": n_2, "n_common": n_common}

        # convert data dicts to DataFrames for easier filtering
        self.data1 = pd.DataFrame(self.data1_tmp)
        self.data2 = pd.DataFrame(self.data2_tmp)
        del(self.data1_tmp)
        del(self.data2_tmp)

    def add_data(self, data1, data2, chrom: str, is_mismatch1: bool, is_mismatch2: bool) -> None:
        """
        Add data from each input dataset to temporary storage for comparison.

        Args:
            data1: Data from the first dataset.
            data2: Data from the second dataset.
            chrom (str): Chromosome identifier.
            is_mismatch1 (bool): Whether there is a mismatch in the first dataset.
            is_mismatch2 (bool): Whether there is a mismatch in the second dataset.

        Returns:
            None
        """
        for col in list(self.data1_tmp.keys())[:-2]:
            self.data1_tmp[col].append(float(data1[self.IDX[col]]))
            self.data2_tmp[col].append(float(data2[self.IDX[col]]))
        self.data1_tmp["chrom"].append(chrom)
        self.data2_tmp["chrom"].append(chrom)
        self.data1_tmp["is_mismatch"].append(is_mismatch1)
        self.data2_tmp["is_mismatch"].append(is_mismatch2)

    def perform_test(self, set1, set2):
        """
        Perform Wilcoxon rank-sum test for multiple columns between two datasets.

        Args:
            set1 (pd.DataFrame): First dataset for comparison.
            set2 (pd.DataFrame): Second dataset for comparison.

        Returns:
            dict: Dictionary containing p-values for each column.
        """
        results = {}
        for col in ["n_a_rel", "n_c_rel", "n_g_rel", "n_t_rel", 
                    "n_del_rel", "n_ins_rel", "n_ref_skip_rel", 
                    "perc_mismatch", "q_mean"]:
            try:
                results[col] = wilcoxon(set1[col], set2[col])[1]
            except ValueError as e:
                warnings.warn(f"Error in column '{col}': {e}. Setting results to None")
                results[col] = pd.NA #(None, None)
        return results

    def stat_comp(self) -> None:
        """
        Perform statistical tests on the data and store the results.
        """
        res_overall = self.perform_test(self.data1.drop(["chrom", "is_mismatch"], axis=1), 
                                        self.data2.drop(["chrom", "is_mismatch"], axis=1))    

        res_by_chrom = {}
        for chrom in self.data1["chrom"].unique():
            res_by_chrom[chrom] = self.perform_test(self.data1.loc[self.data1["chrom"] == chrom], 
                                                    self.data2.loc[self.data2["chrom"] == chrom])  

        res_overall = pd.DataFrame(res_overall, index=["Total"])
        res_by_chrom = pd.DataFrame(res_by_chrom).T

        self.results = pd.concat([res_overall, res_by_chrom])

    def write_results(self) -> None:
        """
        Write the results of statistical tests to the output file.
        """
        self.results.to_csv(self.outpath, sep="\t")

    def write_html(self) -> None:
        table = self.get_results_html()
        time = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')

        template = f"""
        <!DOCTYPE html>
        <html lang="en">
        <head>
            <meta charset="UTF-8">
            <meta http-equiv="X-UA-Compatible" content="IE=edge">
            <meta name="viewport" content="width=device-width, initial-scale=1.0">
            <title>Document</title>
            <style>
            body, h1, h2, h3, h4, p, ul, li {{
                margin: 0;
                padding: 0;
            }}

            body {{
                display: flex;
                justify-content: center;
                align-items: center;
                font-family: 'Arial', sans-serif;
                line-height: 1.6;
                min-height: 100vh;
                flex-direction: column;
                background-color: #f5f5f5;
            }}

            header {{
                background-color: #333;
                color: white;
                padding: 1em 0;
                text-align: center;
                width: 100%;
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

            table {{
                width: 100%;
                border-collapse: collapse;
                font-family: Arial, sans-serif;
            }}

            .table-box {{
                margin: 2em;
                border: 1px solid #ccc;
            }}

            th, td {{
                padding: 8px;
                text-align: center;
                border: 1px solid #ccc;
            }}
            th {{
                background-color: #f2f2f2;
            }}
            td:hover {{
                filter: brightness(85%);
                border: 1px solid #ccc;
            }}
            </style>
        </head>
        <body>
            <header>
                <h1>Statistical comparison summary</h1>
                <p>Produced by <a href="https://github.com/dietvin/neet">Neet</a> on <b>{time}</b></p>
            </header>

            <section>
                <p class="intro-text">
                    Below are the results of the statistical comparisons of different features between <b>{self.path_1}</b> 
                    and <b>{self.path_2}</b>. <b>{self.stats['n_common']}</b> positions that were present in both samples were 
                    used for the comparison. A Wilcoxon rank sum test for connected samples was used with alpha = <b>{self.alpha}</b>.
                </p>
            </section>

            <section>
                <h2>Overview of p-values</h1>
                <p>
                    Table shows the p-values from the statistical comparisons for different extracted features. 
                    Statistical tests were performed for both all samples together (<b>first row</b>), as well as 
                    for each chromosome (<b>remaining rows</b>). Cells are colored if the test suggests significant
                    differences with alpha = <b>{self.alpha}</b>. The darker the color, the smaller the p-value.
                </p>
                <div class="table-box">
                    {table}
                </div>
            </section>

            <footer></footer>
        </body>                
        </html>
        """
        out_name = os.path.splitext(self.outpath)[0]+"_summary.html"
        with open(out_name, "w") as out:
            out.write(template)

    def get_results_html(self) -> str:
        data = self.results
        colormap = plt.cm.Greens_r
        def color_background_font(value):
            try:
                if value < self.alpha:
                    norm_value = (value - data.min().min()) / (self.alpha - data.min().min())
                    color = mcolors.to_hex(colormap(norm_value))
                    font_color = 'white' if norm_value < 0.5 else 'black'
                    return f'background-color: {color}; color: {font_color};'
            except:
                return 'color: {font_color};'
            return ''
        data.rename(columns={"n_a_rel": "Fraction A", "n_c_rel": "Fraction C", "n_g_rel": "Fraction G", "n_t_rel": "Fraction T",
                             "n_del_rel": "Deletion Rate", "n_ins_rel": "Insertion Rate", "n_ref_skip_rel": "Reference skip Rate", "perc_mismatch": "Mismatch Rate",
                             "q_mean": "Mean Q-score"}, inplace=True)
        return data.style.applymap(color_background_font).to_html(float_format=lambda x: '{:e}'.format(x))

    def main(self) -> None:
        """
        Main method to orchestrate the comparison and statistical testing process.
        """
        sys.stdout.write("Searching for shared positions... ")
        self.find_common_pos()
        sys.stdout.write(f"Found {self.stats['n_common']} shared positions.\nPerforming statistical tests... ")
        self.stat_comp()
        if self.write_tsv:
            sys.stdout.write(f"Done.\nWriting results to TSV... ")
            self.write_results()
        sys.stdout.write(f"Done.\nWriting results to HTML... ")
        self.write_html()
        sys.stdout.write(f"Done.\nFinished.")


if __name__=="__main__":
    comparer = StatComparer("/home/vincent/masterthesis/data/epinano_data/processed/curlcake_m6a_extracted.tsv",
                            "/home/vincent/masterthesis/data/epinano_data/processed/curlcake_unm_extracted.tsv",
                            "/home/vincent/masterthesis/data/epinano_data/processed/stat_comp.tsv")
    comparer.main()
    