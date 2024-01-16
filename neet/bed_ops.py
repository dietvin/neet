import argparse, os, csv
from typing import List, Tuple, Dict
from . import helper_functions as hs
from tqdm import tqdm

def tsv_to_bed(in_path: str, out_path: str) -> None:
    """
    Converts a TSV (Tab-Separated Values) file to a BED (Browser Extensible Data) file format.

    Args:
        in_path (str): Path to the input TSV file.
        out_path (str): Path to the output BED file.

    Returns:
        None: This function does not return any value. The converted data is written to the output BED file.
    """
    hs.print_update(f"Reading from '{in_path}'.")
    hs.check_input_path(in_path, [".tsv"])
    hs.check_output_path(out_path, [".bed"])

    n_lines = hs.get_num_lines(in_path)

    with open(in_path, "r") as i, open(out_path, "w") as o:
        next(i)
        for line in tqdm(i, total=n_lines, desc="Transforming TSV to BED"):
            line = line.strip("\n").split("\t")
            o.write(f"{line[0]}\t{int(line[1])-1}\t{int(line[1])}\n")
    hs.print_update(f"Finished. Wrote {n_lines} to '{out_path}'")

def intersect_beds(file_a: str,
                   file_b: str,
                   out_path: str,
                   label_a: str, 
                   label_b: str) -> None:
    """
    Intersects two BED files and writes the results to an output file.

    Args:
        file_a (str): Path to the first input BED file.
        file_b (str): Path to the second input BED file.
        out_path (str): Path to the output file.
        label_a (str, optional): Label for file A. If not provided, the filename of file A without extension will be used. Defaults to None.
        label_b (str, optional): Label for file B. If not provided, the filename of file B without extension will be used. Defaults to None.

    Returns:
        None: This function does not return any value. The results are written to the output file.
    """
    hs.print_update(f"Reading from '{file_a}' and '{file_b}'.")
    hs.check_input_path(file_a, [".bed"])
    hs.check_input_path(file_b, [".bed"])
    hs.check_output_path(out_path, [".bed"])

    label_a = label_a
    label_b = label_b

    with open(file_a, 'r') as file1:
        hs.print_update(f"Reading coordinates from file '{file_a}'... ", line_break=False)
        file_a_coordinates = set(tuple(line.strip().split('\t')[:3]) for line in file1)
        hs.print_update("Done.", with_time=False)
    # Read coordinates from file 2
    with open(file_b, 'r') as file2:
        hs.print_update(f"Reading coordinates from file '{file_b}'... ", line_break=False)
        file_b_coordinates = set(tuple(line.strip().split('\t')[:3]) for line in file2)
        hs.print_update("Done.", with_time=False)

    # Find shared coordinates
    shared_coordinates = file_a_coordinates.intersection(file_b_coordinates)
    n_shared = len(shared_coordinates)
    # Find coordinates exclusive to file 1
    exclusive_file_a = file_a_coordinates.difference(file_b_coordinates)
    n_a_excl = len(exclusive_file_a)
    # Find coordinates exclusive to file 2
    exclusive_file_b = file_b_coordinates.difference(file_a_coordinates)
    n_b_excl = len(exclusive_file_b)

    hs.print_update(f"Found {n_shared} shared positions, {n_a_excl} exclusive to A and {n_b_excl} exclusive to B. Writing to output...")

    with open(out_path, "w") as out:
        for line in shared_coordinates:
            out.write("\t".join(line)+f"\t{label_a}+{label_b}\n")
        for line in exclusive_file_a:
            out.write("\t".join(line)+f"\t{label_a}\n")
        for line in exclusive_file_b:
            out.write("\t".join(line)+f"\t{label_b}\n")

    hs.print_update(f"Finished. Wrote {n_shared+n_a_excl+n_b_excl} coordinates to '{out_path}'.")

def get_coord_names(path: str, keep_name_pos: bool = True) -> Tuple[set[Tuple[str, int]], Dict[Tuple[str, int], str]]:
    """
    Extracts coordinate names from a file and returns a tuple containing two elements:
    
    Parameters:
    - path (str): The path to the file containing coordinate information.
    - keep_name_pos (bool, optional): If True, includes coordinate names in the result; 
      if False, only includes the coordinate positions. Default is True.

    Returns:
    Tuple[Set[Tuple[str, int]], Dict[Tuple[str, int], str]]: A tuple containing two elements.
        - The first element is a set of coordinate positions represented as tuples (chromosome, position).
        - The second element is a dictionary mapping each coordinate position to its corresponding name,
          if `keep_name_pos` is True. If `keep_name_pos` is False, the dictionary is empty.
    """
    with open(path, "r") as file:
        file_pos = []
        name_pos = {}
        for line in file:
            line = line.strip().split("\t")
            chromosome = line[0]
            start = int(line[1])
            end = int(line[2])

            if keep_name_pos:
                for i in range(start, end):
                    file_pos.append((chromosome, i))

                    if len(line) >= 4:
                        name_pos[(chromosome, i)] = line[3]
            else:
                for i in range(start, end):
                    file_pos.append((chromosome, i))

        if keep_name_pos:
            return set(file_pos), name_pos
        else:
            return set(file_pos), {}

def add_bed_info(tsv_file: str, bed_file: str, out_file: str) -> None:
    """
    Reads data from a TSV file and a BED file, combines the information, and writes the updated data to an output file.

    Args:
        tsv_file (str): Path to the TSV file.
        bed_file (str): Path to the BED file.
        out_file (str): Path to the output file.

    Returns:
        None
    """
    hs.print_update(f"Reading from '{tsv_file}' and '{bed_file}'")
    hs.check_input_path(tsv_file, extensions=[".tsv"])
    hs.check_input_path(bed_file, extensions=[".bed"])
    hs.check_output_path(out_file, extensions=[".tsv"])

    bed_data = {}
    with open(bed_file, "r") as bed:
        data = []
        for row in bed:
            row = row.strip().split("\t")
            chromosome = row[0]
            start = int(row[1]) #+ 1 # get it to 1-indexed to fit tsv file
            data = row[3:]
            bed_data[(chromosome, start)] = data
        n_cols = len(data)

    with open(tsv_file, "r") as tsv, open(out_file, "w") as out:
        tsv_reader = csv.reader(tsv, delimiter="\t")
        out_writer = csv.writer(out, delimiter="\t")
        header = next(tsv_reader)
        new_cols = [f"bed_col_{i}" for i in range(1,n_cols+1)]
        header += new_cols

        n_lines = hs.get_num_lines(tsv_file) - 1
        count = 0

        out_writer.writerow(header)
        for row in tqdm(tsv_reader, total=n_lines):
            count += 1
            chromosome = row[0]
            site = int(row[1])
            bed_entry = bed_data.get((chromosome, site), None)
            if bed_entry:
                row = row + bed_entry
            out_writer.writerow(row)

    hs.print_update(f"Finished. Wrote {count} lines to '{out_file}'")

def merge(file_paths: str, output_file_path: str):
    """
    Merges multiple BED files into a single BED file.

    Parameters:
    - file_paths (str): Comma-separated list of paths to the input BED files to be merged.
    - output_file_path (str): Path to the output BED file where the merged data will be written.

    Returns:
    None

    This function reads multiple BED files, merges their data based on overlapping coordinates, 
    and writes the result to a new BED file. If the input BED files have names associated with the 
    coordinates, the output file will include a comma-separated list of names.

    Note: The function assumes the input BED files are tab-delimited and have the following format:
    `chromosome   start   end   name`
    """
    merged_data = {}
    file_path_list = file_paths.split(",")
    hs.print_update(f"Merging {len(file_path_list)} bed files.")

    for path in file_path_list:
        hs.check_input_path(path, [".bed"])

    hs.check_output_path(output_file_path, [".bed"])

    for file_path in file_path_list:
        n_lines = hs.get_num_lines(file_path)
        with open(file_path, "r") as bed:
            for line in tqdm(bed, total=n_lines, desc=f"Processing file '{file_path}'"):
                line = line.strip().split("\t")
                if len(line) >= 3:
                    chromosome = line[0]
                    start = int(line[1])
                    end = int(line[2])
                    name_value = line[3] if len(line) >= 4 else ""

                    for i in range(start, end):
                        pos_key = (chromosome, i)
                        if pos_key not in merged_data:
                            merged_data[pos_key] = [chromosome, i, i+1, name_value]
                        elif len(name_value) > 0:
                            merged_data[pos_key][3] += f",{name_value}"
    
    sorted_data = sorted(merged_data.values(), key=lambda x: (x[0], x[1]))
    with open(output_file_path, 'w') as output_file:
        for entry in sorted_data:
            output_file.write(f"{entry[0]}\t{entry[1]}\t{entry[2]}\t{entry[3]}\n")

    hs.print_update(f"Finished. Wrote {len(sorted_data)} lines to '{output_file_path}'")
    
def difference(bed1: str, bed2: str, out: str):
    """
    Finds the exclusive coordinates in the first BED file that do not overlap with the second BED file
    and writes the result to a new BED file.

    Parameters:
    - bed1 (str): Path to the first BED file.
    - bed2 (str): Path to the second BED file.
    - out (str): Path to the output BED file.

    Returns:
    None

    This function reads two BED files, finds the coordinates exclusive to the first file, and writes
    the result to a new BED file. If the input BED files have names associated with the coordinates,
    the output file will include those names.

    Note: The function uses the `get_coord_names` function to extract coordinate positions and names
    from the input BED files.
    """
    hs.print_update(f"Reading from '{bed1}' and '{bed2}'.")
    hs.check_input_path(bed1, [".bed"])
    hs.check_input_path(bed2, [".bed"])
    hs.check_output_path(out, [".bed"])
    file1_pos, name_pos = get_coord_names(bed1)
    file2_pos, _ = get_coord_names(bed2, keep_name_pos=False)
    
    exclusive_file1 = file1_pos.intersection(file2_pos)
    del(file1_pos)
    del(file2_pos)

    with open(out, "w") as out_file:
        count = 0
        for coordinate in tqdm(exclusive_file1):
            count += 1
            has_name = coordinate in name_pos
            name = name_pos[coordinate] if has_name else ""
            out_file.write(f"{coordinate[0]}\t{coordinate[1]}\t{coordinate[1]+1}\t{name}\n")
    hs.print_update(f"Finished. Wrote {count} lines to '{out}'")