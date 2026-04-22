import os
import numpy as np


def read_file_folder_v_vel(filepath, time):
    """
    Read vertex velocity/force data from file

    Parameters:
    -----------
    filepath : str
        Full path to the output folder
    time : int
        Time step for filename

    Returns:
    --------
    cForces : numpy array
        Array of vertex forces (columns 3:4 from file, 0-indexed as 2:4)
    """
    original_dir = os.getcwd()

    try:
        os.chdir(filepath)

        # Read the matrix data
        data = np.loadtxt(f"v_vel_{time}.dat")

        # Extract columns 3:4 (0-indexed: 2:4)
        cForces = data[:, 2:4]

    finally:
        # Change back to original directory
        os.chdir(original_dir)

    return cForces
