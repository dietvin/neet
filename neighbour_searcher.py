from typing import List
import argparse
from tqdm import tqdm
from pyfiglet import Figlet
from itertools import takewhile, repeat

class NeighbourSearcher:

    input_path: str
    output_path: str
    window_size: int


    def __init__(self, in_path: str, out_path: str, window_size: int) -> None:
        self.input_path = in_path
        self.output_path = out_path
        self.window_size = window_size

    def sort_extractor_file(self) -> None:
        """
        Read and sort a given tsv file from the pileup_extractor.
        Overwrite the existing (unsorted) file.
        """
        pass # for now I assume that the output files are sorted by default

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
        ref_chr = ref[0]
        ref_site = int(ref[1])
        del nb[k]

        # for each neighbour check if they are 1.on the same chr and 2.neighbours
        has_nb = False
        nb_info = ""

        for pos in nb:
            pos = pos.strip("\n").split("\t")
            if pos[0] == ref_chr: # check if same chromosome
                site = int(pos[1])
                relative_pos = site - ref_site
                if abs(relative_pos) <= k: # check if pos are close to each other
                    has_nb = True
                    nb_info += str(relative_pos)+","

        ref_str += f"\t{has_nb}\t{nb_info}\n"
        return ref_str

    def read_lines_sliding_window(self):
        """
        Read and process each line of a tsv file from the pileup extractor
        """
        window_size = 1 + 2 * self.window_size

        n_lines = self.get_num_lines(self.input_path)

        with open(self.input_path, 'r') as file, open(self.output_path, "w") as o:
            progress_bar = tqdm(total=n_lines)
            lines = []
            for line in file:
                lines.append(line)

                if len(lines) > window_size:
                    lines.pop(0)
                
                if len(lines) == window_size:
                    outline = self.process_neighbourhood(lines)
                    o.write(outline)

                progress_bar.update()
            
            progress_bar.close()

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
    parser = argparse.ArgumentParser(prog="Pileup feature extractor",
                                        description="Extracs different characteristics from a\
                                        given pileup file.")
    parser.add_argument('-i', '--input', type=str, required=True,
                        help='Path to the input file')

    parser.add_argument('-o', '--output', type=str, required=True,
                        help='Path to the output file')

    parser.add_argument('-w', '--window_size', type=int, required=True,
                        help='Size of the sliding window = 2*w+1')
    parser.add_argument('-t', '--num_processes', type=int, required=False,
                        help='Number of threads to use for processing.')

    args = parser.parse_args()

    neighbour_searcher = NeighbourSearcher(args.input, args.output, args.window_size)
    
    neighbour_searcher.read_lines_sliding_window()
