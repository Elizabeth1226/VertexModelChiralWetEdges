import os


def read_file_folder_X(filepath, filename):
    """
    Read .dat file where each row starts with count of elements

    Parameters:
    -----------
    filepath : str
        Full path to the output folder
    filename : str
        Filename to read

    Returns:
    --------
    row_data : list of lists
        Each element is a list containing the numeric vector of elements for that row
    """
    original_dir = os.getcwd()

    try:
        os.chdir(filepath)

        # Initialize list to hold each row's data
        row_data = []

        # Open and read file
        with open(filename, 'r') as fid:
            for line in fid:
                line = line.strip()
                if line and not line.startswith('#'):
                    # Convert string to numeric array
                    nums = [float(x) for x in line.split()]
                    if nums:
                        # First element is the count
                        n = int(nums[0])
                        # Store remaining elements in list
                        row_data.append(nums[1:n+1])

    finally:
        # Change back to original directory
        os.chdir(original_dir)

    return row_data
