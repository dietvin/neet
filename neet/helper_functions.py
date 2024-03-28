from typing import Any, List, Tuple, Dict
import argparse, os, warnings, sys
from itertools import takewhile, repeat
import datetime

import pkg_resources
import importlib.resources as impresources
import summary_style


def print_update(message: str, line_break: bool = True, with_time: bool = True) -> None:
    time = get_time()
    if line_break:
        message += "\n"
    if with_time:
        sys.stdout.write(f"{time} | {message}")
    else:
        sys.stdout.write(f"{message}")

def get_time() -> str:
    """
    Returns the current time in the format 'YYYY-mm-DD HH:MM:SS'
    """
    return datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')

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
    time = get_time()
    def stdout_progress(n: int):
        if give_update:
            sys.stdout.write(f"\r{time} | Processing reference genome from file '{path}'. Sequences found: {n}")
            sys.stdout.flush()
            
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




from threading import Thread, Lock
from queue import Queue
from tqdm import tqdm
import tempfile, logging

class SortedFileMerge:
    """
    Class to merge sorted temporary files. Taken from a list of paths to these temporary files,
    each thread takes two files, combines these and adds the newly created temporary file back
    into the pool. This way each iteration the total number of temporary files decreases by 1.
    Once two files are processed, they get deleted. In the end one final temporary file remains 
    containing all sites sorted by chromosome and coordinate.

    Attributes:
        tmpfiles (List[str]): A list of paths to sorted temporary files.
        remaining (int): Number of temporary files remaining to be merged.
        remaining_lock (Lock): Lock for synchronizing access to `remaining`.
        n_threads (int): Number of threads to use for merging.
    """
    tmpfiles: List[str]
    remaining: int
    remaining_lock: Lock
    n_threads: int

    logger: logging.Logger

    def __init__(self, tmpfiles: List[str], n_threads: int = 2) -> None:
        """
        Initializes SortedFileMerge with a list of temporary files and number of threads.

        Args:
            tmpfiles (List[str]): A list of paths to sorted temporary files.
            n_threads (int, optional): Number of threads for merging. Defaults to 2.
        """
        self.tmpfiles = tmpfiles
        self.remaining = len(tmpfiles)
        self.remaining_lock = Lock()
        self.n_threads = n_threads

        self.logger = logging.getLogger("SortedFileMerge")

    def start(self) -> str:
        """
        Starts the merging process and returns the path of the final merged file.

        Returns:
            str: Path to the final merged file.
        """
        tmp_queue = Queue()
        for i in self.tmpfiles:
            tmp_queue.put(i)

        log_queue = Queue()
        log_thread = Thread(target=self.log, args=(log_queue,))
        log_thread.start()

        progress_queue = Queue()
        progress_thread = Thread(target=self.progress, args=(progress_queue, len(self.tmpfiles)))
        progress_thread.start()

        threads = []
        for i in range(self.n_threads):
            t = Thread(target=self.combine_two_tmps, args=(tmp_queue, progress_queue, log_queue))
            t.start()
            threads.append(t)     

        for t in threads:
            t.join()

        progress_queue.put(False)
        progress_thread.join()

        log_queue.put(None)
        log_thread.join()

        return tmp_queue.get()

    def progress(self, progress_queue: Queue, n_files: int) -> None:
        """
        Tracks the progress of merging temporary files.

        Args:
            progress_queue (Queue): Queue to communicate progress.
            n_files (int): Total number of files to merge.
        """
        with tqdm(desc=f"{get_time()} | Merging temporary files", total=n_files-1) as progress:
            increase = progress_queue.get()        
            while increase:
                progress.update()
                increase = progress_queue.get()

    def log(self, log_queue: Queue) -> None:
        msg = log_queue.get()
        while msg:
            self.logger.debug(msg)
            msg = log_queue.get()            

    def combine_two_tmps(self, tmpfiles_queue: Queue, progress_queue: Queue, log_queue: Queue) -> None:
        """
        Combines two temporary files into one until all files are merged.

        Args:
            tmpfiles_queue (Queue): Queue containing paths to temporary files.
            progress_queue (Queue): Queue to communicate progress.
        """
        while self.remaining > 1:
            with self.remaining_lock:
                self.remaining -= 1

            tmp_path1 = tmpfiles_queue.get()
            tmp_path2 = tmpfiles_queue.get()
            
            with open(tmp_path1, "r") as f1, open(tmp_path2, "r") as f2, tempfile.NamedTemporaryFile(mode="w", prefix="neet_", delete=False) as new_tmp:
                line1 = f1.readline()
                line2 = f2.readline()
                while line1 != "" and line2 != "":
                    idx1, idx2 = int(line1.split("|")[0]), int(line2.split("|")[0])
                    # check which row is earlier and write the earlier line to file, then update the written file
                    if idx1 < idx2:
                        new_tmp.write(line1)
                        line1 = f1.readline()
                    else:
                        new_tmp.write(line2)
                        line2 = f2.readline()
                        
                # write the remaining lines to the new file
                if line1 == "":
                    while line2 != "":
                        new_tmp.write(line2)
                        line2 = f2.readline()
                else:
                    while line1 != "":
                        new_tmp.write(line1)
                        line1 = f1.readline()

                os.unlink(tmp_path1)
                os.unlink(tmp_path2)

                tmpfiles_queue.put(new_tmp.name)
                progress_queue.put(True)
                
                log_queue.put(f"Combined '{tmp_path1}' and '{tmp_path2}' into '{new_tmp.name}'")