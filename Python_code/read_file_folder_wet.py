import os
import numpy as np


def read_file_folder_wet(filepath, file):
    """
    Read file from folder

    Parameters:
    -----------
    filepath : str
        Full path to the output folder
    file : str
        Filename to read

    Returns:
    --------
    A, B : numpy arrays
        Two columns of data from the file
    """
    original_dir = os.getcwd()

    try:
        os.chdir(filepath)

        # Read the file, skipping comment lines starting with '#'
        # Data is space-separated
        data = np.loadtxt(file, comments='#')

        # Extract columns
        if data.ndim == 1:
            # Single row case
            A = np.array([data[0]])
            B = np.array([data[1]])
        else:
            A = data[:, 0]
            B = data[:, 1]

    finally:
        # Change back to original directory
        os.chdir(original_dir)

    return A, B
