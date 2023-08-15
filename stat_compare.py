from typing import Dict, List, Tuple
from scipy.stats import wilcoxon
import pandas as pd
import warnings, sys

class StatComparer:
    path_1: str
    path_2: str
    outpath: str
    alpha: float
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

    def __init__(self, path_1: str, path_2: str, outpath: str, alpha: float = 0.01) -> None:
        self.path_1 = path_1
        self.path_2 = path_2
        self.outpath = outpath
        self.alpha = alpha
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

    def process_paths(self):
        pass

    def check_in_paths(self):
        pass

    def get_out_path(self):
        pass

    def find_common_pos(self) -> None:
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
        for col in list(self.data1_tmp.keys())[:-2]:
            self.data1_tmp[col].append(float(data1[self.IDX[col]]))
            self.data2_tmp[col].append(float(data2[self.IDX[col]]))
        self.data1_tmp["chrom"].append(chrom)
        self.data2_tmp["chrom"].append(chrom)
        self.data1_tmp["is_mismatch"].append(is_mismatch1)
        self.data2_tmp["is_mismatch"].append(is_mismatch2)

    def perform_test(self, set1, set2):
        results = {}
        for col in ["n_a_rel", "n_c_rel", "n_g_rel", "n_t_rel", 
                    "n_del_rel", "n_ins_rel", "n_ref_skip_rel", 
                    "perc_mismatch", "q_mean"]:
            try:
                results[col] = wilcoxon(set1[col], set2[col])[1]
            except ValueError as e:
                warnings.warn(f"Error in column '{col}': {e}. Setting results to None")
                results[col] = None #(None, None)
        return results

    def stat_comp(self) -> None:
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
        self.results.to_csv(self.outpath, sep="\t")

    def main(self) -> None:
        sys.stdout.write("Searching for shared positions... ")
        self.find_common_pos()
        sys.stdout.write(f"Found {self.stats['n_common']} shared positions.\nPerforming statistical tests...")
        self.stat_comp()
        self.write_results()


if __name__=="__main__":
    comparer = StatComparer("/home/vincent/masterthesis/data/epinano_data/processed/curlcake_m6a_extracted.tsv",
                            "/home/vincent/masterthesis/data/epinano_data/processed/curlcake_unm_extracted.tsv",
                            "/home/vincent/masterthesis/data/epinano_data/processed/stat_comp.tsv")
    comparer.main()
    