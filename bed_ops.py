import argparse, os

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
                   label_a: str = None, 
                   label_b: str = None) -> None:
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

    args = parser.parse_args()

    if args.subcommand == "tsv_to_bed":
        tsv_to_bed(args.input, args.output)
    elif args.subcommand == "intersect_bed":
        intersect_beds(args.file_a, args.file_b, args.output, args.label_a, args.label_b)
    else:
        parser.print_help()