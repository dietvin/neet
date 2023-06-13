import re
from typing import Dict, List, Tuple, Union
import numpy as np
import argparse

class FeatureExtractor:
    input_path : str
    output_path : str
    ref_sequences : Dict[str, str]

    def __init__(self, in_path: str, out_path: str, ref_path: str) -> None:
        self.input_path = in_path
        self.output_path = out_path
        self.ref_sequences = self.get_references(ref_path)

    def get_references(self, path: str) -> Dict[str, str]:
        """
        Reads a fasta file and stores the sequences in a dictionary (values) with the 
        corresponding chromosome names (keys).

        Parameters
        ----------
        path : str
            filepath to a fasta file

        Returns
        -------
        dict[str]
            Dictionary where the key is the chromosome name and the value is the sequence
        """
        with open(path, "r") as ref:
            lines = ref.readlines()
            i = 0
            refs = {}
            for i in range(len(lines)):
                line = lines[i]
                if line.startswith(">"):
                    refs[line[1:].strip()] = lines[i+1].strip()
        
        return refs
    
    def process_file(self):
        """
        Reads .pileup file and processes it, writing the results to a new file.

        Parameters
        ----------
        infile : str
            path to the input .pileup file
        outfile : str
            path to the output tsv file
        ref : str
            path to the reference fasta
            
        Returns
        -------
        None
        """
        with open(self.input_path, "r") as i:
            lines = i.readlines()

            with open(self.output_path, "w") as o:
                header = f"chr\tsite\tn_reads\tref_base\tmajority_base\tn_a\tn_c\tn_g\tn_t\tn_a_rel\tn_c_rel\tn_g_rel\tn_t_rel\tmotif\tperc_mismatched\tq_mean\tq_std\n"
                o.write(header)
                for line in lines:
                    outline = self.process_position(line.split("\t"))
                    o.write(outline) 

    def process_position(self, line: List[str]) -> str:
        """
        Takes a line from a pileup file and processes it.

        Parameters
        ----------
        line : list[str]
            list containing each element from the pileup line.

        Returns
        -------
        str
            New line derived from the initial one. Can be written to a new file in consequent
            steps.
        """
        # extract elements from list
        chr, site, ref_base, n_reads, read_bases, read_qualities = line[0], int(line[1]), line[2], int(line[3]), line[4], line[5]
        
        # get reference sequence 
        ref = self.ref_sequences[chr]
        # get absolute number of A, C, G, T, ins, del
        count = self.parse_pileup_string(read_bases, ref_base)

        # get relative number of A, C, G and T counts
        count = self.get_relative_count(count, n_reads)

        # get majority base
        majority_base = self.get_majority_base(count)

        # get 11b motif
        motif = self.get_motif(chr, site, ref, k=5)

        # get allele fraction
        allele_fraction = self.get_allele_fraction(count, ref_base)

        # get qualitiy measures
        quality_mean, quality_std = self.get_read_quality(read_qualities)

        out = f'{chr}\t{site}\t{n_reads}\t{ref_base}\t{majority_base}\t{count["a"]}\t{count["c"]}\t{count["g"]}\t{count["t"]}\t{count["a_rel"]}\t{count["c_rel"]}\t{count["g_rel"]}\t{count["t_rel"]}\t{motif}\t{allele_fraction}\t{quality_mean}\t{quality_std}\n'
        return out

    def remove_indels(self, pileup_string: str) -> str:
        """
        Takes a pileup string and removes all occurences of the following patterns:
        '\+[0-9]+[ACGTNacgtn]+' for insertions
        '\-[0-9]+[ACGTNacgtn]+' for deletions

        Parameters
        ----------
        pileup_string : str
            Pileup string extracted from the fifth column of a pileup file

        Returns
        -------
        str
            Pileup strings with all occurences of the patterns above removed
        """
        pattern = "(\+|\-)[0-9]+[ACGTNacgtn]+"
        
        # get the start and end indices of all found patterns 
        coords = []
        for m in re.finditer(pattern, pileup_string):
            str_len = int(pileup_string[m.start(0)+1]) + 1
            coords.append((m.start(0), m.start(0)+1+str_len))
            
        # remove the patterns by the indices
        for start, end in reversed(coords): # reverse list as to not shift the index downstream
            pileup_string = pileup_string[:start] + pileup_string[end:]

        return pileup_string

    def parse_pileup_string(self, pileup_string: str, ref_base: str) -> Dict[str, Union[str, int]]:
        """
        Extracts the number of each base called at a given position, as well as the number
        of insertions and deletions. Information is extracted from a pileup string (fifth
        column in a pileup file).

        Parameters
        ----------
        pileup_string : str
            Pileup string extracted from the fifth column of a pileup file
        ref_base : str
            reference base at the position corresponding to the pileup string

        Returns
        -------
        dict
            Dictionary containing the number of A, T, C, G, 
            insertions and deletions.
        """
        pileup_string = pileup_string.lower()
        # remove all occurences of a caret and the following letter (could contain a,c,g,t)
        pileup_string = re.sub(r'\^.', '', pileup_string)

        ref_base = ref_base.lower()
        count_dict = {"a": 0, "t": 0, "c": 0, "g": 0, "ins": 0, "del": 0}

        # get number of insertions
        count_dict["ins"] = len(re.findall(r'\+[0-9]+[ACGTNacgtn]+', pileup_string))

        # get number of deletions
        count_dict["del"] = len(re.findall(r'\-[0-9]+[ACGTNacgtn]*|\*', pileup_string))

        # remove indel patterns to count the number of mismatches correctly
        pileup_string = self.remove_indels(pileup_string)

        # get number of mismatches (i.e. [ACGT])
        count_dict["a"] = pileup_string.count("a")
        count_dict["t"] = pileup_string.count("t")
        count_dict["c"] = pileup_string.count("c")
        count_dict["g"] = pileup_string.count("g")

        # get number of matches (determine where to count matches bases on ref_base)
        n_matches = pileup_string.count('.') + pileup_string.count(',')
        count_dict[ref_base] = n_matches

        return count_dict

    def get_relative_count(self, count_dict: Dict[str, Union[str, int]], n_reads: int) -> Dict[str, Union[str, int, float]]:
        """
        Gets a dictionary containing the absolute counts for A, C, G and T 
        and calculates the relative proportions

        Parameters
        ----------
        count_dict : dict[int]
            Dictionary containing the absolute counts for A, C, G and T
        n_reads : int
            Number of reads at the given position

        Returns
        -------
        dict[float]
            Dictionary containing the relative counts for A, C, G and T
        """
        try:
            count_dict["a_rel"] = count_dict["a"] / n_reads
            count_dict["c_rel"] = count_dict["c"] / n_reads
            count_dict["g_rel"] = count_dict["g"] / n_reads
            count_dict["t_rel"] = count_dict["t"] / n_reads
        except ZeroDivisionError:
            count_dict["a_rel"] = 0
            count_dict["c_rel"] = 0
            count_dict["g_rel"] = 0
            count_dict["t_rel"] = 0

        return count_dict

    def get_majority_base(self, count_dict: Dict[str, Union[str, int, float]]) -> str:
        """
        Gets a dictionary containing the absolute counts for A, C, G and T and returns the
        key of the one with the highest count.

        Parameters
        ----------
        count_dict : dict
            dictionary containing the absolute counts for A, C, G and T

        Returns
        -------
        str
            Key from the dictionary corresponding to the largest value
        """
        dict_subset = dict((k, count_dict[k]) for k in ("a", "c", "g", "t"))
        return max(dict_subset, key = dict_subset.get).upper()

    def get_motif(self, chr: str, site: int, ref: str, k: int) -> str:
        """
        Extracts the motif of k bases up- and downstream from a given chromosomal site.
        Around the start and end of a refernce sequence the missing bases are filled with
        Ns.

        Parameters
        ----------
        chr : str
            name of the chromosome
        site : int
            position on the chromosome (1-indexed)
        ref : str
            reference sequence for the given chromosome 
        k : int
            number of bases to be regarded in both up- and downstream direction 
            
        Returns
        -------
        str
            sequence of k bases around the center site
        """ 
        idx = site-1
        n_ref = len(ref)

        if idx >= 0 and idx < n_ref:
            idx_l = idx-k
            idx_r = idx+k+1
            # left overhang
            if idx_l < 0:
                len_overhang = abs(idx_l)
                overhang = "N" * len_overhang
                motif = overhang + ref[:idx_r]
            # right overhang
            elif idx_r > n_ref:
                len_overhang = idx_r - n_ref
                overhang = "N" * len_overhang
                motif = ref[idx_l:] + overhang
            # no overhang
            else:
                motif = ref[idx_l:idx_r]

            return motif
        
    def get_allele_fraction(self, count_dict: Dict[str, Union[str, int, float]], ref_base: str) -> int:
        """
        Calculates the number of mismatched reads

        Parameters
        ----------
        count_dict : dict
            Dictionary containing the number of occurences of A,C,G,T for a given position
        ref_base : str
            reference base at the given position

        Returns
        -------
        int
            Number of mismatched reads a the given position
        """
        mismatch_count_sum = 0
        for b in ["a", "c", "g", "t"]:
            if b != ref_base.lower():
                mismatch_count_sum += count_dict[b+"_rel"]

        return mismatch_count_sum

    def get_read_quality(self, read_qualities: str) -> Tuple[float, float]:
        """
        Calculates the mean and std from the read qualities given in the sixth row
        of a pileup file.

        Parameters
        ----------
        read_qualities : str
            Read quality string from pileup file

        Returns
        -------
        tuple[float, float]
            Mean and standard deviation of read qualities
        """
        # transform string to list of corresponding phred numeric values
        vals = [code - 33 for code in read_qualities.encode("ascii")]

        mean = sum(vals)/len(vals)
        std = np.std(vals)

        return mean, std 
    

if __name__ == "__main__":
    parser = argparse.ArgumentParser(prog="Pileup feature extractor",
                                        description="Extracs different characteristics from a\
                                        given pileup file.")
    parser.add_argument('-i', '--input', type=str, required=True,
                        help='Path to the input file')

    parser.add_argument('-o', '--output', type=str, required=True,
                        help='Path to the output file')

    parser.add_argument('-r', '--reference', type=str, required=True,
                        help='Path to the reference file')
    args = parser.parse_args()

    feature_extractor = FeatureExtractor(args.input, args.output, args.reference)
    feature_extractor.process_file()       