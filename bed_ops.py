import argparse, os, csv

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


def filter_with_bed(tsv_file: str, bed_file: str, out_file: str) -> None:
    bed_data = []
    with open(bed_file, "r") as bed:
        for row in bed:
            row = row.strip().split("\t")
            chromosome = row[0]
            start = int(row[1]) #+ 1 # get it to 1-indexed to fit tsv file
            bed_data.append((chromosome, start))
    bed_data = set(bed_data)
    
    with open(tsv_file, "r") as tsv, open(out_file, "w") as out:
        header = next(tsv)
        out.write(header)

        for row_str in tsv:
            row = row_str.strip().split("\t")
            chromosome = row[0]
            start = int(row[1])
            if (chromosome, start) in bed_data:
                out.write(row_str)





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


    filter_parser = subparsers.add_parser("filter_tsv", prog="Filter TSV by BED", description="Filter existing TSV file based on positions in BED file.")
    filter_parser.add_argument("-i", "--input", type=str, required=True,
                                help="TSV output from PileupExtractor/NeighbourhoodSearcher")
    filter_parser.add_argument("-o", "--output", type=str, required=True,
                                help="Path to the output file")
    filter_parser.add_argument("-b", "--bed", type=str, required=True,
                                help="BED file containing single positions")


    args = parser.parse_args()

    if args.subcommand == "tsv_to_bed":
        tsv_to_bed(args.input, args.output)
    elif args.subcommand == "intersect_bed":
        intersect_beds(args.file_a, args.file_b, args.output, args.label_a, args.label_b)
    elif args.subcommand == "add_bed_info":
        add_bed_info(args.input, args.bed, args.output)
    elif args.subcommand == "filter_tsv":
        filter_with_bed(args.input, args.bed, args.output)
    else:
        parser.print_help()