import os
import numpy as np


def read_file_folder_cell_bin_categories(filepath, file):
    """
    Read cell bin categories from file

    Parameters:
    -----------
    filepath : str
        Full path to the output folder
    file : str
        Filename to read

    Returns:
    --------
    bin_num : numpy array
        Array of bin numbers
    """
    original_dir = os.getcwd()

    try:
        os.chdir(filepath)

        # Read the data using numpy, skipping comment lines starting with '#'
        # Data is space-separated
        bin_num = np.loadtxt(file, comments='#', dtype=int)

    finally:
        # Change back to original directory
        os.chdir(original_dir)

    return bin_num
