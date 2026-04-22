import os
import numpy as np


def read_file_folder_vel(filepath, time):
    """
    Read cell velocity/force data from file

    Parameters:
    -----------
    filepath : str
        Full path to the output folder
    time : int
        Time step for filename

    Returns:
    --------
    cForces : numpy array
        Array of cell forces (columns 4:5 from file, 0-indexed as 3:5)
    """
    original_dir = os.getcwd()

    try:
        os.chdir(filepath)

        # Read the matrix data
        data = np.loadtxt(f"c_vel_{time}.dat")

        # Extract columns 4:5 (0-indexed: 3:5)
        cForces = data[:, 3:5]

    finally:
        # Change back to original directory
        os.chdir(original_dir)

    return cForces
