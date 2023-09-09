import argparse, os, csv
from typing import List, Tuple, Dict

def tsv_to_bed(in_path: str, out_path: str) -> None:
    """
    Converts a TSV (Tab-Separated Values) file to a BED (Browser Extensible Data) file format.

    Args:
        in_path (str): Path to the input TSV file.
        out_path (str): Path to the output BED file.

    Returns:
        None: This function does not return any value. The converted data is written to the output BED file.
    """
    print(f"Reading from '{in_path}'. Writing to '{out_path}'")
    with open(in_path, "r") as i, open(out_path, "w") as o:
        next(i)
        for line in i:
            line = line.strip("\n").split("\t")
            o.write(f"{line[0]}\t{int(line[1])-1}\t{int(line[1])}\n")

def intersect_beds(file_a: str,
                   file_b: str,
                   out_path: str,
                   label_a: str|None = None, 
                   label_b: str|None = None) -> None:
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
    label_a = label_a if label_a else os.path.splitext(os.path.basename(file_a))[0]
    label_b = label_b if label_b else os.path.splitext(os.path.basename(file_b))[0]

    with open(file_a, 'r') as file1:
        file_a_coordinates = set(tuple(line.strip().split('\t')[:3]) for line in file1)

    # Read coordinates from file 2
    with open(file_b, 'r') as file2:
        file_b_coordinates = set(tuple(line.strip().split('\t')[:3]) for line in file2)

    # Find shared coordinates
    shared_coordinates = file_a_coordinates.intersection(file_b_coordinates)

    # Find coordinates exclusive to file 1
    exclusive_file_a = file_a_coordinates.difference(file_b_coordinates)

    # Find coordinates exclusive to file 2
    exclusive_file_b = file_b_coordinates.difference(file_a_coordinates)

    with open(out_path, "w") as out:
        for line in shared_coordinates:
            out.write("\t".join(line)+f"\t{label_a}+{label_b}\n")
        for line in exclusive_file_a:
            out.write("\t".join(line)+f"\t{label_a}\n")
        for line in exclusive_file_b:
            out.write("\t".join(line)+f"\t{label_b}\n")

def get_coord_names(path: str, keep_name_pos: bool = True) -> Tuple[set[Tuple[str, int]], Dict[Tuple[str, int], str]]:
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

def intersect(bed1: str, bed2: str, out: str):
    file1_pos, name1_pos = get_coord_names(bed1)
    file2_pos, name2_pos = get_coord_names(bed2)
    
    shared_coordinates = file1_pos.intersection(file2_pos)
    del(file1_pos)
    del(file2_pos)

    with open(out, "w") as out_file:
        for coordinate in shared_coordinates:
            has_name_in1 = coordinate in name1_pos
            has_name_in2 = coordinate in name2_pos

            if has_name_in1 & has_name_in2:
                name = f"{name1_pos[coordinate]},{name2_pos[coordinate]}"
            elif has_name_in1:
                name = {name1_pos[coordinate]}
            elif has_name_in2:
                name = {name2_pos[coordinate]}
            else:
                name = ""
            out_file.write(f"{coordinate[0]}\t{coordinate[1]}\t{coordinate[1]+1}\t{name}\n")

def difference(bed1: str, bed2: str, out: str):
    file1_pos, name_pos = get_coord_names(bed1)
    file2_pos, _ = get_coord_names(bed2, keep_name_pos=False)
    
    exclusive_file1 = file1_pos.intersection(file2_pos)
    del(file1_pos)
    del(file2_pos)

    with open(out, "w") as out_file:
        for coordinate in exclusive_file1:
            has_name = coordinate in name_pos
            name = name_pos[coordinate] if has_name else ""
            out_file.write(f"{coordinate[0]}\t{coordinate[1]}\t{coordinate[1]+1}\t{name}\n")


def merge(file_paths: str, output_file_path: str):
    merged_data = {}
    file_path_list = file_paths.split(",")

    for file_path in file_path_list:
        with open(file_path, "r") as bed:
            for line in bed:
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
        for entries in sorted_data:
            for entry in entries:
                output_file.write(f"{entry[0]}\t{entry[1]}\t{entry[2]}\t{entry[3]}\n")


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

        out_writer.writerow(header)
        for row in tsv_reader:
            chromosome = row[0]
            site = int(row[1])
            bed_entry = bed_data.get((chromosome, site), None)
            if bed_entry:
                row = row + bed_entry
            out_writer.writerow(row)




if __name__=="__main__":
    parser = argparse.ArgumentParser(prog="Bed operations",
                                     description="Provides different BED functionalities")
    
    subparsers = parser.add_subparsers(dest="subcommand")

    tsv_bed_parser = subparsers.add_parser("tsv_to_bed", prog="TSV to BED",
                                           description="Transform TSV output to BED file to be loaded into IGV.")
    tsv_bed_parser.add_argument('-i', '--input', type=str, required=True,
                                help='Path to the input file')
    tsv_bed_parser.add_argument('-o', '--output', type=str, required=True,
                                help='Path to the output file.')

    intersect_parser = subparsers.add_parser("intersect_bed", prog="Intersect BED files", 
                                           description="Extracts shared and exclusive positions from two bed files.")
    intersect_parser.add_argument("-a", "--file_a", type=str, required=True,
                                  help="First BED file")
    intersect_parser.add_argument("-b", "--file_b", type=str, required=True,
                                  help="Second BED file")
    intersect_parser.add_argument('-o', '--output', type=str, required=True,
                                  help='Path to the output file.')
    intersect_parser.add_argument("--label_a", type=str, required=False,
                                  help="Label given to the output for file a")
    intersect_parser.add_argument("--label_b", type=str, required=False,
                                  help="Label given to the output for file b")

    add_bed_parser = subparsers.add_parser("add_bed_info", prog="Add bed info to TSV file from PileupExtractor.",
                                           description="Add the name column from a BED file to the respective line in the TSV file.")
    add_bed_parser.add_argument("-i", "--input", type=str, required=True,
                                help="TSV output from PileupExtractor/NeighbourhoodSearcher")
    add_bed_parser.add_argument("-o", "--output", type=str, required=True,
                                help="Path to the output file")
    add_bed_parser.add_argument("-b", "--bed", type=str, required=True,
                                help="BED file containing additional information in the 'name' column")

    merge_parser = subparsers.add_parser("merge", prog="Merge bed files", description="Merges the position contained in multiple bed files.")
    merge_parser.add_argument("-i", "--input", type=str, required=True,
                                help="Paths to multiple bed files, comma separated. For example: '/path/to/file1.bed,/path/to/file2.bed,...'")
    merge_parser.add_argument("-o", "--output", type=str, required=True,
                                help="Path to the output file")

    merge_parser = subparsers.add_parser("difference", prog="Get difference bed files", 
                                         description="Extract the positons present in file 1 and not present in file 2.")
    merge_parser.add_argument("-i1", "--input1", type=str, required=True,
                                help="Path to bed file 1")
    merge_parser.add_argument("-i2", "--input2", type=str, required=True,
                                help="Path to bed file 2")
    merge_parser.add_argument("-o", "--output", type=str, required=True,
                                help="Path to the output file")


    args = parser.parse_args()

    if args.subcommand == "tsv_to_bed":
        tsv_to_bed(args.input, args.output)
    elif args.subcommand == "intersect_bed":
        intersect_beds(args.file_a, args.file_b, args.output, args.label_a, args.label_b)
    elif args.subcommand == "add_bed_info":
        add_bed_info(args.input, args.bed, args.output)
    elif args.subcommand == "merge":
        merge(args.input, args.output)
    elif args.subcommand == "difference":
        difference(bed1=args.input1, bed2=args.input1, out=args.output)
    else:
        parser.print_help()