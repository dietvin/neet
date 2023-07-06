import argparse, io, sys, re
from helper_functions import check_get_in_path, check_get_out_path
from typing import List, Tuple, Callable

class Filter:
    input_path: str
    output_path: str

    chr: str
    site: List[int]
    n_reads: Tuple[int, Callable[[int, int], bool]]
    base: List[str]
    mismatched: bool
    perc_mismatched: Tuple[float, Callable[[float, float], bool]]
    q_score: Tuple[float, Callable[[float, float], bool]]
    from_stdin: bool
    to_stdout: bool

    OPERATORS = {
        "<": lambda x, y: x < y,
        "<=": lambda x, y: x <= y,
        ">": lambda x, y: x > y,
        ">=": lambda x, y: x >= y,
        "==": lambda x, y: x == y
        }

    def __init__(self, 
                 input_path: str, 
                 output_path: str, 
                 chr: str, 
                 site: str,
                 n_reads: str, 
                 base: str, 
                 mismatched: bool, 
                 perc_mismatched: str, 
                 motif: str, 
                 q_score: str) -> None:
        self.input_path = check_get_in_path(input_path, 
                                            exp_extensions=[".tsv"], 
                                            warn_expected_text="Expected .tsv file")
        self.output_path = check_get_out_path(output_path, self.input_path)
        self.chromosome = chr
        self.site = self.get_sites(site)
        self.n_reads = self.get_n_reads(n_reads)
        self.base = self.get_bases(base)
        self.mismatched = mismatched
        self.perc_mismatched = self.get_val_fun_float(perc_mismatched)
        self.motif = motif
        self.q_score = self.get_val_fun_float(q_score)
        

    def get_sites(self, site_str: str) -> List[int]:
        """Extracts site IDs from a string representation of sites.

        Args:
            site_str (str): A string representing site IDs. It can be a single ID,
                multiple IDs separated by commas, or a range of IDs separated by a dash.

        Returns:
            list[int]: A list of site IDs extracted from the input string.

        Raises:
            Exception: If the site IDs cannot be extracted from the given string.

        Example:
            >>> obj = MyClass()
            >>> obj.get_sites('1')
            [1]
            >>> obj.get_sites('2,3,4')
            [2, 3, 4]
            >>> obj.get_sites('10-15')
            [10, 11, 12, 13, 14, 15]
        """
        try:
            site = [int(site_str)]
        except:
            if "," in site_str:
                site = list(map(int, site_str.split(",")))
            elif "-" in site_str:
                site = site_str.split("-")
                start = int(site[0])
                end = int(site[1])
                site = list(range(start, end+1))
            else:
                raise Exception(f"Could not extract sites from given string '{site_str}'")

        return site                

    def get_n_reads(self, n_reads_str: str) -> Tuple[int, Callable[[int, int], bool]]:
        """
        Extracts the number of reads and the corresponding comparison function from the given string.

        Args:
            n_reads_str (str): The string representing the number of reads.

        Returns:
            Tuple[int, callable]: A tuple containing the extracted number of reads as an integer
                and the corresponding comparison function.

        Raises:
            Exception: If the string does not match any of the supported formats.
        """
        try:
            val = int(n_reads_str)
            fun = self.OPERATORS.get("==")
        except:
            match = re.match(r'([<>]=?|==)(\d+)', n_reads_str)
            if match:
                value = int(match.group(2))
                op = match.group(1)
                fun = self.OPERATORS.get(op)
            else:
                raise Exception(f"Could not extract information from given string '{n_reads_str}'")
        return val, fun
    
    def get_bases(self, base_str: str) -> List[str]:
        """
        Extracts the bases from the given string.

        Args:
            base_str (str): The string representing the bases.

        Returns:
            List[str]: A list of unique bases extracted from the string.

        Raises:
            Exception: If the string contains unexpected bases.
        """
        bases = base_str.split(",")
        bases = list(set(bases)) # remove duplicates
        for base in bases:
            if base not in ["A", "C", "G", "T", "N"]:
                raise Exception(f"Given string for --base flag '{base_str}' contains unexpected base(s). \
                                (Allowed bases: A, C, G, T, N)")
        return bases

    def get_val_fun_float(self, string: str) -> Tuple[float, Callable[[float, float], bool]]:
        """
        Extracts the percentage of mismatched values / mean q_score and the corresponding comparison function
        from the given string.

        Args:
            string (str): The string representing the percentage of mismatched values.

        Returns:
            Tuple[float, Callable[[float, float], bool]]: A tuple containing the extracted percentage
                of mismatched values as a float and the corresponding comparison function.

        Raises:
            Exception: If the string does not match any of the supported formats.
        """
        try:
            val = float(string)
            fun = self.OPERATORS.get("==")
        except:
            match = re.match(r'([<>]=?|==)(\d+\.?\d*)', string)
            if match:
                val = float(match.group(2))
                op = match.group(1)
                fun = self.OPERATORS.get(op)
            else:
                raise Exception(f"Could not extract information from given string '{string}'")
        return val, fun

    def filter_tsv(self) -> None:
        out = None if self.to_stdout else open(self.output_path, "w")

        with open(self.input_path, "r") as file:
            for row in file:
                if self.passes_filter():
                    self.output_line(row, out)
        
        if not self.to_stdout:
            out.close()

    def output_line(self, line: str, output: io.TextIOWrapper = None) -> None:
        if output:
            output.write(line)
        else:
            sys.stdout.write(line)

    def passes_filter(self) -> bool:
        pass

if __name__ == "__main__":
    parser = argparse.ArgumentParser(prog="TSV filter",
                                     description="Filter TSV output from PileupExtractor by given values.")
    parser.add_argument("-i", "--input", type=str, required=False,
                        help="Path to input TSV file. If none is given, read from stdin.")
    parser.add_argument("-o", "--output", type=str, required=False,
                        help="Path to outptu TSV file. If none is given, write from stdout.")
    parser.add_argument("-c", "--chromosome", type=str, required=False,
                        help="Filter by given chromosome.")
    parser.add_argument("-s", "--site", type=str, required=False,
                        help="Filter by given site(s) or range. For single site: 'x'; for multiple sites: 'x,y,z,...'; for range: 'x-y'")
    parser.add_argument("-n", "--n_reads", type=str, required=False,
                        help="Filter by coverage. To filter coverage >= x: 'x'; coverage <= x: '<=x'; coverage == x: '==x'")
    parser.add_argument("-b", "--base", type=str, required=False,
                        help="Filter by reference base(s). To filter single base (e.g. A): 'A'; multiple bases (e.g. A, C & T): 'A,C,T'")
    parser.add_argument("-m", "--mismatched", action="store_true", type=str, required=False,
                        help="Filter mismatched positions.")
    parser.add_argument("-p", "--percent_mismatched", type=str, required=False,
                        help="Filter by percent of mismatched reads. To filter perc_mismatched >= x: 'x'; perc_mismatched <= x: '<=x'")
    parser.add_argument("-f", "--motif", type=str, required=False,
                        help="Filter by motif around position.")
    parser.add_argument("-q", "--q_score", type=str, required=False,
                        help="Filter by mean quality. To filter q_mean >= x: 'x'; q_mean <= x: '<=x'")

    