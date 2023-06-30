import argparse

def tsv_to_bed(in_path: str, out_path: str) -> None:
    print(f"Reading from '{in_path}'. Writing to '{out_path}'")
    with open(in_path, "r") as i, open(out_path, "w") as o:
        next(i)
        for line in i:
            line = line.strip("\n").split("\t")
            o.write(f"{line[0]}\t{int(line[1])-1}\t{int(line[1])}\n")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(prog="TSV to BED",
                                     description="Transform TSV output to BED file to be loaded into IGV.")
    parser.add_argument('-i', '--input', type=str, required=True,
                        help='Path to the input file')
    parser.add_argument('-o', '--output', type=str, required=True,
                        help='Path to the output file.')
    args = parser.parse_args()
    
    tsv_to_bed(args.input, args.output)