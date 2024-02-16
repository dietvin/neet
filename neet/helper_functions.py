from typing import Any, List, Tuple, Dict
import argparse, os, warnings, sys
from itertools import takewhile, repeat
import datetime

import pkg_resources
import importlib.resources as impresources
from . import summary_style


def print_update(message: str, line_break: bool = True, with_time: bool = True) -> None:
    time = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    if line_break:
        message += "\n"
    if with_time:
        sys.stdout.write(f"{time}  |  {message}")
    else:
        sys.stdout.write(f"{message}")


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


def check_input_path(path: str, extensions: List[str]) -> None:
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

def check_output_path(path: str, extensions: List[str]) -> None:
    """
    Check if the specified base directory to the file exists and if the file has the expected extension.
    It created the directory if it does not exist.

    Parameters:
        path (str): The file path to be checked.
        extensions (List[str]): A list of expected file extensions (e.g., ['.txt', '.csv']).

    Raises:
        FileNotFoundError: If the specified file path does not exist.
        Warning: If the file extension is not among the expected extensions.
    """
    basedir = os.path.dirname(path)
    check_create_dir(basedir)

    file_type = os.path.splitext(path)[1]
    if not file_type in extensions:
        warnings.warn(f"Found file extension {file_type}. Expected file extension to be one of: {extensions}. If this is deliberate, ignore warning.", Warning)


def check_create_dir(dirname: str):
    """
    Check if the directory exists, and if not, attempt to create it.

    Parameters:
    - dirname (str): The directory path to check and create.

    Raises:
    - Exception: If the directory creation fails.
    """
    if not os.path.isdir(dirname): # does directory of the given file exist? If not try to create it
        print_update(f"Creating directory {dirname}")
        try: 
            os.makedirs(dirname)
        except Exception as e:
            raise Exception(f"Could not create directory '{dirname}'. Error: {e}")
 
def is_directory(path: str) -> bool:
    """
    Checks if path corresponds to a file (with an extension) or a directory.
    Note that it does not check if the path exists, only if the format is 
    fitting to either.
    
    Returns:
        bool: True if the path is a directory, False if it is a file
    """
    _, extension = os.path.splitext(path)
    if not extension or path.endswith(os.path.sep):
        return True
    return False

def process_outpath(out: str, filename: str|None, ext: List[str]) -> str:
    """
    Process the output path and return a valid output file path.

    Parameters:
    - out (str): The specified output path.

    Raises:
    - FileNotFoundError: If the specified output directory does not exist.
    - Warning: If the output file has an extension other than '.html'. A warning is issued, but the function continues execution.

    Returns:
    - str: The processed output file path.
    """
    if is_directory(out):
        check_create_dir(out)
        return os.path.join(out, filename)
    else: 
        check_output_path(out, ext)
        return out


def get_references(path: str, give_update: bool = True) -> Dict[str, str]:
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

    def stdout_progress(n: int):
        if give_update:
            sys.stdout.write(f"\rSequences found: {n}")
            sys.stdout.flush()

    if give_update:
        print_update(f"Processing reference genome from file '{path}'.", line_break=False)
    with open(path, "r") as ref:
        refs = {}
        line = next(ref)
        if not line.startswith(">"):
            raise Exception(f"Fasta format error. The first line of fasta file '{path}' does not contain a header (starting with '>').")
        
        chr_name = line[1:].strip().split(" ")[0]
        seq = ""

        chr_found = 1
        stdout_progress(chr_found)

        for line in ref:
            if line.startswith(">"):
                refs[chr_name] = seq
                chr_name = line[1:].strip().split(" ")[0]
                seq = ""

                chr_found += 1
                stdout_progress(chr_found)

            else:
                seq += line.strip().upper()
                
        refs[chr_name] = seq # add the last dict entry
        sys.stdout.write("\n")
    return refs

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


def load_html_template_str() -> Tuple[str, str]:
    """
    Load static files style.css and plotly.js from folder summary_style. Returns both as strings.
    """
    # https://stackoverflow.com/questions/6028000/how-to-read-a-static-file-from-inside-a-python-package
    css_style = (impresources.files(summary_style) / "style.css")
    plotly_js = (impresources.files(summary_style) / "plotly.js")

    with css_style.open("r") as css:
        css_string = css.read()
    with plotly_js.open("r") as plotly_js:
        plotly_js_string = plotly_js.read()

    return css_string, plotly_js_string