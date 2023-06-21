from pileup_extractor import FeatureExtractor
import csv, os, warnings, argparse
from tqdm import tqdm
from pyfiglet import Figlet

class TwoSampleExtractor(FeatureExtractor):
    """
    A class for extracting and merging two TSV files based on position.

    Attributes
    ----------
    input_path_1 : str
        Path to the first input TSV file.
    input_path_2 : str
        Path to the second input TSV file.
    output_path : str
        Path to the output merged TSV file.
    basename_1 : str, optional
        Basename of the first input file (default is 'a').
    basename_2 : str, optional
        Basename of the second input file (default is 'b').

    Methods
    -------
    __init__(in_path_1: str, in_path_2: str, out_path: str)
        Initializes the TwoSampleExtractor instance.
    get_basename(in_path: str) -> str
        Extracts the basename of a file from the given path.
    merge_tsv_files()
        Merges two large TSV files based on position.
    """
    input_path_1: str
    input_path_2: str
    output_path: str
    basename_1: str = "a"
    basename_2: str = "b"

    def __init__(self, in_path_1: str, in_path_2: str, out_path: str) -> None:
        self.input_path_1 = self.check_get_in_path(in_path_1)
        self.input_path_2 = self.check_get_in_path(in_path_2)
        self.output_path = self.check_get_out_path(out_path, self.input_path_1, "_merged.tsv")
        self.basename_1 = self.get_basename(in_path_1)
        self.basename_2 = self.get_basename(in_path_2)
        
    def get_basename(self, in_path: str) -> str:
        """
        Extracts the basename of a file from the given path.

        Parameters
        ----------
        in_path : str
            The input file path.

        Returns
        -------
        str 
            The basename of the file.

        Warnings
        --------
            If the file extension is not ".tsv", a warning is issued.

        Example
        -------
            >>> file_path = "/path/to/file.tsv"
            >>> basename = get_basename(file_path)
            >>> print(basename)
            "file"
            
        """
        file_name = os.path.basename(in_path)
        base_name, file_extension = os.path.splitext(file_name)

        if file_extension != ".tsv":
            warnings.warn(f"Given input file '{file_name}' has extension '{file_extension}'. Note that the input expects '.tsv' format.")

        return base_name

    def merge_tsv_files(self):
        """
        Merge two large TSV files based on the position given in the "chr" and "site" columns.
        The merged file will contain the columns from both input files, with suffixes "_a" and "_b"
        added to indicate the source file.

        Parameters
        ----------
        file1 : str
            Path to the first input TSV file.
        file2 : str
            Path to the second input TSV file.
        output_file : str
            Path to the output merged TSV file.

        Returns
        -------
        None 
            The merged TSV file is created at the specified output path.
        """
        file1 = self.input_path_1
        file2 = self.input_path_2

        f = Figlet(font="slant")
        print(f.renderText("Neet - file merger"))

        progress_bar = tqdm(total=self.get_num_lines(file1)+self.get_num_lines(file2))

        output_file = self.output_path
        # Open input files in read mode, open output file in write mode
        with open(file1, "r") as file1_handle, open(file2, "r") as file2_handle, open(output_file, "w") as output_handle:
            # Create CSV reader and writer objects
            reader1 = csv.reader(file1_handle, delimiter="\t")
            reader2 = csv.reader(file2_handle, delimiter="\t")
            writer = csv.writer(output_handle, delimiter="\t")

            # Read the header row from both input files
            header1 = next(reader1)
            header2 = next(reader2)

            # Find the column indexes of "chr" and "site" columns in both input files
            chr_index1 = header1.index("chr")
            site_index1 = header1.index("site")
            chr_index2 = header2.index("chr")
            site_index2 = header2.index("site")

            # Write the merged header row to the output file
            output_header = header1 + [col + "_b" for col in header2 if col != "chr" and col != "site"]
            writer.writerow(output_header)

            # Initialize current_row1 and current_row2 with the first rows from both input files
            current_row1 = next(reader1, None)
            progress_bar.update()
            current_row2 = next(reader2, None)
            progress_bar.update()
            
            # Merge the rows from both input files based on "chr" and "site" columns
            while current_row1 and current_row2:
                chr1 = current_row1[chr_index1]
                chr2 = current_row2[chr_index2]
                site1 = int(current_row1[site_index1])
                site2 = int(current_row2[site_index2])

                if chr1 == chr2 and site1 == site2:
                    # Rows have matching positions, write the merged row to the output file
                    output_row = current_row1 + [col for i, col in enumerate(current_row2) if i != chr_index2 and i != site_index2]
                    writer.writerow(output_row)
                    # Move to the next row in both input files
                    current_row1 = next(reader1, None)
                    progress_bar.update()
                    current_row2 = next(reader2, None)
                    progress_bar.update()                    


                elif chr1 < chr2 or (chr1 == chr2 and site1 < site2):
                    # Row 2 has a later position (i.e. position is missing in file 1)
                    writer.writerow(current_row1 + [""] * (len(header2) - 2))
                    # Move to the next row in the first input file
                    current_row1 = next(reader1, None)
                    progress_bar.update()

                else:
                    # Row 1 has a later position (i.e. position is missing in file 2)
                    output_row = [chr2, site2] + [""] * (len(header1) - 2) + [col for i, col in enumerate(current_row2) if i != chr_index2 and i != site_index2]
                    writer.writerow(output_row)
                    # Move to the next row in the second input file
                    current_row2 = next(reader2, None)
                    progress_bar.update()

            # Write remaining rows from file1, if any
            while current_row1:
                writer.writerow(current_row1 + [""] * (len(header2) - 2))
                current_row1 = next(reader1, None)
                progress_bar.update()


            # Write remaining rows from file2, if any
            while current_row2:
                chr2 = current_row2[chr_index2]
                site2 = int(current_row2[site_index2])

                output_row = [chr2, site2] + [""] * (len(header1) - 2) + [col for i, col in enumerate(current_row2) if i != chr_index2 and i != site_index2]
                writer.writerow(output_row)
                current_row2 = next(reader2, None)
                progress_bar.update()

            progress_bar.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(prog="Pileup feature extractor",
                                        description="Merges two files containing extracted features as given by PileupExtractor or NeighbourhoodSearcher.")
    parser.add_argument('-a', '--file_a', type=str, required=True,
                        help='Path to the first file from PileupExtractor or NeighbourhoodSearcher')
    parser.add_argument('-b', '--file_b', type=str, required=True,
                        help='Path to the second file from PileupExtractor or NeighbourhoodSearcher')
    parser.add_argument('-o', '--output', type=str, required=False,
                        help='Path to the output file')
    args = parser.parse_args()

    feature_extractor = TwoSampleExtractor(args.file_a, args.file_b, args.output)
    feature_extractor.merge_tsv_files()
