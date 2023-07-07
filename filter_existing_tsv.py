import argparse, io, sys, re
from helper_functions import check_get_in_path, check_get_out_path
from typing import List, Tuple, Callable

class Filter:
    """
    Class for filtering TSV files based on specified criteria.

        Attributes:
        input_path (str): Path to the input TSV file.
        output_path (str): Path to the output TSV file.
        chr (str): Filter by chromosome.
        site (List[int]): Filter by site(s) or range.
        n_reads (Tuple[int, Callable[[int, int], bool]]): Filter by coverage.
        base (List[str]): Filter by reference base(s).
        mismatched (bool): Filter mismatched positions.
        perc_mismatched (Tuple[float, Callable[[float, float], bool]]): Filter by percent of mismatched reads.
        motif (str): Filter by motif around position.
        q_score (Tuple[float, Callable[[float, float], bool]]): Filter by mean quality.
        OPERATORS (dict): Dictionary mapping comparison operators to their corresponding lambda functions.

    Methods:
        get_sites(site_str: str) -> List[int]:
            Extracts site IDs from a string representation of sites.
        get_n_reads(n_reads_str: str) -> Tuple[int, Callable[[int, int], bool]]:
            Extracts the number of reads and the corresponding comparison function from the given string.
        get_bases(base_str: str) -> List[str]:
            Extracts the bases from the given string.
        get_val_fun_float(string: str) -> Tuple[float, Callable[[float, float], bool]]:
            Extracts the percentage of mismatched values and the corresponding comparison function from the given string.
        filter_tsv() -> None:
            Filters the TSV file based on the specified criteria.
        output_line(line: str, output: io.TextIOWrapper = None) -> None:
            Writes a line to the output file or stdout.
        passes_filter(row: str) -> bool:
            Checks if a row passes the specified filter criteria.
    """
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
        """
        Initializes the Filter object.

        Args:
            input_path (str): Path to the input TSV file.
            output_path (str): Path to the output TSV file.
            chr (str): Filter by chromosome.
            site (str): Filter by site(s) or range.
            n_reads (str): Filter by coverage.
            base (str): Filter by reference base(s).
            mismatched (bool): Filter mismatched positions.
            perc_mismatched (str): Filter by percent of mismatched reads.
            motif (str): Filter by motif around position.
            q_score (str): Filter by mean quality.
        """
        self.input_path = check_get_in_path(input_path, 
                                            exp_extensions=[".tsv"], 
                                            warn_expected_text="Expected .tsv file") if input_path else None
        self.output_path = check_get_out_path(output_path, self.input_path) if output_path else None

        self.from_stdin = False if input_path else True
        self.to_stdout = False if output_path else True

        self.chr = chr
        self.site = self.get_sites(site)
        self.n_reads = self.get_n_reads(n_reads)
        self.base = self.get_bases(base)
        self.mismatched = mismatched
        self.perc_mismatched = self.get_val_fun_float(perc_mismatched)
        self.motif = motif.upper() if motif else None
        self.q_score = self.get_val_fun_float(q_score)
        

    def get_sites(self, site_str: str) -> List[int]:
        """
        Extracts site IDs from a string representation of sites.

        Args:
            site_str (str): A string representing site IDs. It can be a single ID,
                multiple IDs separated by commas, or a range of IDs separated by a dash.

        Returns:
            list[int]: A list of site IDs extracted from the input string.
            None if given string is None

        Raises:
            Exception: If the site IDs cannot be extracted from the given string.
        """
        if site_str is None:
            return None
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
            None if given string is None

        Raises:
            Exception: If the string does not match any of the supported formats.
        """
        if n_reads_str is None:
            return None
        try:
            val = int(n_reads_str)
            fun = self.OPERATORS.get(">=")
        except:
            match = re.match(r'([<>]=?|==)(\d+)', n_reads_str)
            if match:
                value = int(match.group(2))
                op = match.group(1)
                fun = self.OPERATORS.get(op)
            else:
                raise Exception(f"Could not extract information from given string '{n_reads_str}'")
        return value, fun
    
    def get_bases(self, base_str: str) -> List[str]:
        """
        Extracts the bases from the given string.

        Args:
            base_str (str): The string representing the bases.

        Returns:
            List[str]: A list of unique bases extracted from the string.
            None if given string is None

        Raises:
            Exception: If the string contains unexpected bases.
        """
        if base_str is None:
            return None
        
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
            None if given string is None
                
        Raises:
            Exception: If the string does not match any of the supported formats.
        """
        if string is None:
            return None

        try:
            val = float(string)
            fun = self.OPERATORS.get(">=")
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
        """
        Filters the TSV file based on the specified criteria.
        """
        out = None if self.to_stdout else open(self.output_path, "w")

        with open(self.input_path, "r") as file:
            header = next(file)
            self.output_line(header, out)
            for row in file:
                if self.passes_filter(row):
                    self.output_line(row, out)
        
        if not self.to_stdout:
            out.close()

    def output_line(self, line: str, output: io.TextIOWrapper = None) -> None:
        """
        Writes a line to the output file or stdout.

        Args:
            line (str): The line to be written.
            output (io.TextIOWrapper, optional): The output file. Defaults to None.
        """
        if output:
            output.write(line)
        else:
            sys.stdout.write(line)

    def passes_filter(self, row: str) -> bool:
        """
        Checks if a row passes the specified filter criteria.

        Args:
            row (str): A row from the TSV file.

        Returns:
            bool: True if the row passes the filter, False otherwise.
        """

        def check_func(val, tup: Tuple[int | float, Callable[[int | float, int | float], bool]]) -> bool:
            val_ref = tup[0]
            fun = tup[1]

            return fun(val, val_ref)

        row = row.strip("\n").split("\t")
        chr = row[0]
        site = int(row[1])
        n_reads = int(row[2])
        ref_base = row[3]
        maj_base = row[4]
        perc_mismatched = float(row[17])
        motif = row[18]
        q_score = float(row[19])

        if self.chr:
            if chr != self.chr: return False
        if self.site: 
            if site not in self.site: return False
        if self.n_reads:
            if not check_func(n_reads, self.n_reads): return False
        if self.base:
            if ref_base not in self.base: return False
        if self.mismatched:
            if maj_base == ref_base: return False
        if self.perc_mismatched:
            if not check_func(perc_mismatched, self.perc_mismatched): return False
        if self.motif:
            if motif != self.motif: return False
        if self.q_score:
            if not check_func(q_score, self.q_score): return False
        return True

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
    parser.add_argument("-m", "--mismatched", action="store_true", required=False,
                        help="Filter mismatched positions.")
    parser.add_argument("-p", "--percent_mismatched", type=str, required=False,
                        help="Filter by percent of mismatched reads. To filter perc_mismatched >= x: 'x'; perc_mismatched <= x: '<=x'")
    parser.add_argument("-f", "--motif", type=str, required=False,
                        help="Filter by motif around position.")
    parser.add_argument("-q", "--q_score", type=str, required=False,
                        help="Filter by mean quality. To filter q_mean >= x: 'x'; q_mean <= x: '<=x'")

    args = parser.parse_args()

    filter = Filter(input_path=args.input,
                    output_path=args.output,
                    chr=args.chromosome,
                    site=args.site,
                    n_reads=args.n_reads,
                    base=args.base,
                    mismatched=args.mismatched,
                    perc_mismatched=args.percent_mismatched,
                    motif=args.motif,
                    q_score=args.q_score)
    filter.filter_tsv()