import os


def read_file_folder_cVertices(filepath, file):
    """
    Read cell vertices data from file

    Parameters:
    -----------
    filepath : str
        Full path to the output folder
    file : str
        Filename to read

    Returns:
    --------
    data : list of lists
        Each element is a list containing row data for that cell
    """
    original_dir = os.getcwd()

    try:
        os.chdir(filepath)

        # Initialize empty list to store the data
        data = []

        # Read the file row-wise
        with open(file, 'r') as fid:
            for line in fid:
                line = line.strip()
                if line:
                    # Convert the line to a list of numbers
                    row_data = [float(x) for x in line.split()]
                    # Store the row data in the list
                    data.append(row_data)

    finally:
        # Change back to original directory
        os.chdir(original_dir)

    return data
