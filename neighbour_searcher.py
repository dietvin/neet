from typing import List, Tuple, Any
import argparse, sys, warnings, os, io
from tqdm import tqdm
from pyfiglet import Figlet
from itertools import takewhile, repeat
from helper_functions import float_between_zero_and_one

class NeighbourSearcher:
    """
    From a PileupExtractor output file, search for errors in surrounding positions. Does this by using
    a sliding window of size 2*window_size+1 across each line of a given file. For a given central 
    line of a window, check for each neighbouring line if it is also a neighbour on the reference sequence
    and optionally if it is considered an error based on the error_threshold.
    The output file contains the same information as the input file, with two added columns containing 
    (1.) if a error is found in the neighbourhood of a given position and (2.) where exactly these errors occur,
    denoted by comma-separated relative positions.

    Attributes
    ----------
    input_path : str
        path to an PileupExtractor output file
    output_path : str
        path to newly created output file
    window_size : int
        size of the window used to find neighbours
    error_threshold : float
        error threshold to filter by. Neighbours are only considered if the error percentage is 
        larger than the threshold (i.e. if a neighbour is considered an error).

    Methods
    -------
    - sort_extractor_file() -> None
    - read_lines_sliding_window() -> None
    - process_neighbourhood(neighbourhood: List[str]) -> str
    - get_num_lines(path: str) -> int:
    """

    input_path: str
    output_path: str
    window_size: int
    error_threshold: float

    def __init__(self, window_size: int, in_path: str, out_path: str, err_thresh: float = None) -> None:
        self.input_path = self.check_get_in_path(in_path) if in_path else None
        self.output_path = self.check_get_out_path(out_path, in_path) if out_path else None
        self.window_size = window_size
        self.error_threshold = err_thresh if err_thresh is not None else 0

    def __str__(self) -> str:
        o = f"NeighbourhoodSearcher instance information:\n\n"
        o += f" - input file path: {self.input_path}\n"
        o += f" - writing output file to {self.output_path}\n"
        o += f" - searching neighbours {self.window_size} bases up- and downstream from a given position\n"
        o += f" - neighbouring positions must have an error percentage of >={self.error_threshold} to be regarded as errors\n"
        return o
    
    def sort_extractor_file(self) -> None:
        """
        Read and sort a given tsv file from the pileup_extractor.
        Overwrite the existing (unsorted) file.
        """
        pass # for now I assume that the output files are sorted by default

    def check_get_in_path(self, in_path: str) -> str:
        """
        Check if the given input path is valid and of the expected file type (.tsv).

        Parameters
        ----------
        in_path : str
            Path to the input file given by the user.

        Returns
        -------
        str
            Valid input file path.

        Raises
        ------
        FileNotFoundError
            If the input file does not exist.

        Warns
        -----
        UserWarning
            If the input file is not of the expected file type (.tsv).
            Warns the user to ensure it is a tab-separated file as produced by the 
            PileupExtractor.
        """
        if not os.path.exists(in_path): # does file exist?
            raise FileNotFoundError(f"Input file not found. File '{in_path}' does not exist.")
        file_type = os.path.splitext(in_path)[1]
        if file_type != ".tsv": # is file likely in pileup format?
            warnings.warn(f"Input file of type {file_type}. Make sure that this is a tab-separated file (.tsv)", Warning)
        
        return in_path

    def check_get_out_path(self, out_path: str, in_path: str) -> str:
        """
        Check if the given out_put path is valid. Can be either a filename or directory.
        If a directory is given, output path will be '[DIR-path]/[INPUT-FILE-BASENAME]_neighbours.tsv'.

        Parameters
        ----------
        out_path : str
            Output path given by the user. Either path to a (non-existing) file or a directory
        in_path : str
            Path to input file given by the user

        Returns
        -------
        str
            Valid path to output file 

        Raises
        ------
        FileNotFoundError
            If the output directory or path to the output file does not exist.
        """
        # check if directory/file exists
        if os.path.isdir(out_path):
            if not os.path.exists(out_path):
                raise FileNotFoundError(f"Output directory not found. '{out_path}' does not exist.")

            in_basename = os.path.splitext(os.path.basename(in_path))[0]
            if not out_path.endswith("/"):
                out_path += "/"

            return os.path.join(out_path, in_basename + "_neighbours.tsv")
        
        else:
            dirname = os.path.dirname(out_path)
            if not os.path.exists(dirname):
                raise FileNotFoundError(f"Path to output file not found. '{dirname}' does not exist.")

            file_extension = os.path.splitext(out_path)[1]
            if file_extension != ".tsv":
                warnings.warn(f"Given output file has extension '{file_extension}'. Note that the output file will be of type '.tsv'.")

            return out_path


    def read_lines_sliding_window(self) -> None:
        """
        Read and process each line of a tsv file from the pileup extractor
        """

        def output_line(line, output: io.TextIOWrapper = None) -> None:
            if output:
                output.write(line)
            else:
                sys.stdout.write(line)

        def write_edge_lines(neighbourhood: List[str], outfile: io.TextIOWrapper, start: bool = True):
            k = self.window_size
            r = range(k) if start else range(k+1, 2*k+1)
            for current_pos in r:
                outline = self.process_edge(current_pos, neighbourhood, start)
                output_line(outline, outfile)

        def write(file_input, n_lines: int = None):
            window_size = 1 + 2 * self.window_size
            lines = []      
            first = True

            o = None if to_stdout else open(self.output_path, "w")

            header = next(file_input)
            header = header.strip("\n")+"\thas_neighbour_error\tneighbour_error_pos\n" 
            output_line(header, o)

            if not to_stdout:
                progress_bar = tqdm() if from_stdin else tqdm(total=self.get_num_lines(self.input_path) - 1)

            for line in file_input:
                lines.append(line)

                if len(lines) > window_size:  
                    lines.pop(0)
            
                if len(lines) == window_size:
                    # once lines hits the max. size for the first time, get and write the first k lines
                    if first:
                        first = False
                        write_edge_lines(lines, o, start=True)

                    outline = self.process_neighbourhood(lines)
                    output_line(outline, o)

                if not to_stdout:
                    progress_bar.update()

            # after the last center line was added to lines, process the last k lines
            write_edge_lines(lines, o, start=False)

            if not to_stdout:
                o.close()
                progress_bar.close()

        from_stdin = False if self.input_path else True
        to_stdout = False if self.output_path else True

        f = Figlet(font="slant")
        if not to_stdout:
            print(f.renderText("Neet - neighbourhood searcher"))
            print(str(self))

        # Reading data from STDIN
        if from_stdin:
            write(sys.stdin)
        # Reading data from file
        else:
            with open(self.input_path, 'r') as file:
                write(file)

    def process_edge(self, current_pos: int, neighbourhood: List[str], start: bool = True) -> str:
        k = self.window_size

        ref_str = neighbourhood[current_pos].strip("\n")
        nb = neighbourhood.copy()
        ref = nb[current_pos].strip("\n").split("\t")

        if start:
            nb = nb[:current_pos+k+1]
            del(nb[current_pos])
        else:       
            del(nb[current_pos])
            nb = nb[current_pos-k:]

        has_nb, nb_info = self.get_neighbour_info(ref, nb)
        ref_str += f"\t{has_nb}\t{nb_info}\n"
        return ref_str


    def process_neighbourhood(self, neighbourhood: List[str]) -> str:
        """
        Get 2*window_size+1 rows ([row i-k, ..., row i, ..., row i+k]) that 
        are next to each other in the tsv file and compare the row i to all 
        remaining ones. Check if the other rows are on the same chromosome 
        and if the relative distance btw them is smaller or equal to k.
        Create summary string that indicates the error position relative to 
        the center position.
        Add new information to row and return

        Parameters
        ----------
        neighbourhood : List[str]
            List of k number of rows extracted from a tsv file

        Returns
        -------
        str
            New line containing the neighbourhood information for a given
            center position
        """
        k = self.window_size

        ref_str = neighbourhood[k].strip("\n")
        nb = neighbourhood.copy()

        ref = nb[k].strip("\n").split("\t")
        del nb[k]

        has_nb, nb_info = self.get_neighbour_info(ref, nb)
        ref_str += f"\t{has_nb}\t{nb_info}\n"
        return ref_str

    def get_neighbour_info(self, ref: List[str], neighbourhood: List[str]) -> Tuple[bool, str]:
        """
        From a range of neighbouring lines in a file (neighbourhood), check if positions in these
        lines are neighbours on the reference genome (based on chromosome and site) and if close
        genomic positions can be regarded as errors based on the given error threshold.

        Parameters
        ----------
        ref : List[str]
            line of the central position for which neighbours are searched
        neighbourhood : List[str]
            list containing the lines surrounding the central position

        Returns
        -------
        Tuple[bool, str]
            - boolean indicating if any neighbouring error was found
            - string giving the relative distance to all neighbouring errors to the central position,
              if any were found 
        """
        k = self.window_size
        ref_chr = ref[0]
        ref_site = int(ref[1])

        has_nb = False
        nb_info = ""

        # for each neighbour check if they are 1.on the same chr and 2.neighbours
        for pos in neighbourhood:
            pos = pos.strip("\n").split("\t")
            chr = pos[0]
            perc_error = float(pos[18])

            if (chr == ref_chr) & (perc_error >= self.error_threshold): # check if same chromosome & if pos is error
                site = int(pos[1])
                relative_pos = site - ref_site

                if (abs(relative_pos) <= k): # check if pos are close to each other
                    has_nb = True
                    nb_info += str(relative_pos)+","
        return has_nb, nb_info


    def get_num_lines(self, path: str) -> int:
        """
        Calculate the number of lines in a given file. Function taken from
        https://stackoverflow.com/questions/845058/how-to-get-line-count-of-a-large-file-cheaply-in-python
        
        Parameters
        ----------
        path : str
            Path to a file

        Returns
        -------
        int
            Number of lines in the given file
        """
        f = open(path, 'rb')
        bufgen = takewhile(lambda x: x, (f.raw.read(1024*1024) for _ in repeat(None)))
        return sum( buf.count(b'\n') for buf in bufgen )



if __name__ == "__main__":
    parser = argparse.ArgumentParser(prog="Neighbourhood searcher",
                                        description="Adds information about neighbouring error positions.")
    parser.add_argument('-w', '--window_size', type=int, required=True, default=3,
                        help='Size of the sliding window = 2*w+1')
    parser.add_argument('-i', '--input', type=str, required=False,
                        help='Path to the input file')
    parser.add_argument('-o', '--output', type=str, required=False,
                        help='Path to the output file')
    parser.add_argument('-t', '--error_threshold', type=float_between_zero_and_one, required=False,
                        help='Filter by minimum fraction of mismatched/deleted/inserted bases')

    args = parser.parse_args()

    neighbour_searcher = NeighbourSearcher(args.window_size, args.input, args.output, args.error_threshold)
    
    neighbour_searcher.read_lines_sliding_window()
