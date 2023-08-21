from typing import Any, List
import argparse, os, warnings
from itertools import takewhile, repeat
import datetime


def print_update(message: str) -> None:
    time = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    print(f"{time}  |  {message}")


def get_num_lines(path: str) -> int:
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


def check_path(path: str, extensions: List[str]) -> None:
    """
    Check if the specified file path exists and has the expected file extension.

    This function verifies whether the file specified by the given path exists and has a valid extension.
    If the file does not exist, it raises a FileNotFoundError with a detailed error message.
    If the file extension does not match any of the expected extensions, it raises a Warning.

    Parameters:
        path (str): The file path to be checked.
        extensions (List[str]): A list of expected file extensions (e.g., ['.txt', '.csv']).

    Raises:
        FileNotFoundError: If the specified file path does not exist.
        Warning: If the file extension is not among the expected extensions.
    """
    if not os.path.exists(path): # does file exist?
        raise FileNotFoundError(f"Input file not found. File '{path}' does not exist.")
    file_type = os.path.splitext(path)[1]
    if not file_type in extensions:
        warnings.warn(f"Found file extension {file_type}. Expected file extension to be one of: {extensions}. If this is deliberate, ignore warning.", Warning)


def check_get_in_path(in_path: str, 
                      exp_extensions: List[str] = [".msf", ".pup", ".pileup"],
                      warn_expected_text: str = "Expected .pileup file or similar.") -> str:
    """
    Check if the given input path is valid and of the expected file type.

    Parameters
    ----------
    in_path : str
        Path to the input file given by the user.
    exp_extensions : List[str]
        File extensions of the given input format
    warn_expected_text : str
        Text to be displayed in the warning if another format is given
        
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
        If the input file is not of the expected file type (.msf, .pup, .pileup).
        Warns the user to ensure it is a pileup file as produced by Samtools' mpileup function.
    """
    if not os.path.exists(in_path): # does file exist?
        raise FileNotFoundError(f"Input file not found. File '{in_path}' does not exist.")
    file_type = os.path.splitext(in_path)[1]
    if not file_type in exp_extensions: # is file likely in pileup format?
        warnings.warn(f"Input file of type {file_type}. {warn_expected_text}", Warning)
    
    return in_path

def check_get_out_path(out_path: str, in_path: str, suffix: str = "_extracted.tsv") -> str:
    """
    Check if the given out_put path is valid. Can be either a filename or directory.
    If a directory is given, output path will be '[DIR-path]/[INPUT-FILE-BASENAME]_extracted.tsv'.

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
    if os.path.isdir(out_path):
        if not os.path.exists(out_path):
            raise FileNotFoundError(f"Output directory not found. '{out_path}' does not exist.")

        in_basename = os.path.splitext(os.path.basename(in_path))[0]
        if not out_path.endswith("/"):
            out_path += "/"

        return os.path.join(out_path, in_basename + suffix)
    
    else:
        dirname = os.path.dirname(out_path)
        if not os.path.exists(dirname):
            raise FileNotFoundError(f"Path to output file not found. '{dirname}' does not exist.")

        file_extension = os.path.splitext(out_path)[1]
        if file_extension != ".tsv":
            warnings.warn(f"Given output file has extension '{file_extension}'. Note that the output file will be of type '.tsv'.")

        return out_path


def positive_int(value: Any) -> int:
    """
    Convert the given value to an integer and validate that it is positive.

    Parameters
    ----------
    value : int
        Value given in the command line when filtering by number of reads
    
    Returns
    -------
    int
        Same number given, but only if it is a positive integer
    """
    ival = int(value)
    if ival < 0:
        raise argparse.ArgumentTypeError(f"{value} is not a positive integer.")
    return ival

def positive_float(value: Any) -> int:
    """
    Convert the given value to a float and validate that it is positive.

    Parameters
    ----------
    value : Any
        The value to be converted and validated.

    Returns
    -------
    float
        The converted and validated float value.

    Raises
    ------
    argparse.ArgumentTypeError
        If the value is not a positive float 

    """
    fvalue = float(value)
    if fvalue < 0:
        raise argparse.ArgumentTypeError(f"{value} is not a positive float.")
    return fvalue

def float_between_zero_and_one(value: Any) -> float:
    """
    Convert the given value to a float and validate that it is between 0 and 1 (inclusive).

    Parameters
    ----------
    value : Any
        The value to be converted and validated.

    Returns
    -------
    float
        The converted and validated float value.

    Raises
    ------
    argparse.ArgumentTypeError
        If the value is not a float between 0 and 1.

    """
    fvalue = float(value)
    if not 0 <= fvalue <= 1:
        raise argparse.ArgumentTypeError(f"{value} is not a float between 0 and 1")
    return fvalue


if __name__ == "__main__":
    pass