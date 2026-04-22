import os
import numpy as np


def read_file_folder_cell_spin(Nx, kP, P0, C1, dw, di, disorder, tMAX, seed, file):
    """
    Read cell spin data from file

    Parameters:
    -----------
    Nx, kP, P0, C1, dw, di, disorder, tMAX, seed : numeric
        Parameters for folder path construction
    file : str
        Filename to read

    Returns:
    --------
    Cspin : numpy array
        Array of cell spin values
    """
    # Determine folder path based on P0 value
    if P0 in [1, 2, 3, 4, 5]:
        if disorder > 0 and C1 > 0:
            folder_path = f'/Volumes/TOSHIBA_EXT/output/out_Nx_{Nx}_kA_0.5_kP_{kP:.2f}_P0_{int(P0)}_disorder_{disorder:.2f}_tMAX_{tMAX}_dw_{dw}_di_{di}_ca1_{C1:.2f}_ca2_0_iC_0_iP_0.25_iS_1_iT_{seed}/'
        else:
            if kP in [0.1, 0.5]:
                folder_path = f'/Users/wangqianyu/Desktop/VertexModel/vertex_model/chiral_code/output/out_Nx_{Nx}_kA_0.5_kP_{kP:.1f}_P0_{int(P0)}_disorder_0_tMAX_{tMAX}_dw_{dw}_di_{di}_ca1_0_ca2_0_iC_0_iP_0.25_iS_1_iT_{seed}/'
            else:
                folder_path = f'/Users/wangqianyu/Desktop/VertexModel/vertex_model/chiral_code/output/out_Nx_{Nx}_kA_0.5_kP_{kP:.2f}_P0_{int(P0)}_disorder_0_tMAX_{tMAX}_dw_{dw}_di_{di}_ca1_0_ca2_0_iC_0_iP_0.25_iS_1_iT_{seed}/'
    else:
        if disorder > 0 and C1 > 0:
            folder_path = f'/Volumes/TOSHIBA_EXT/output/out_Nx_{Nx}_kA_0.5_kP_{kP:.2f}_P0_{P0:.1f}_disorder_{disorder:.2f}_tMAX_{tMAX}_dw_{dw}_di_{di}_ca1_{C1:.2f}_ca2_0_iC_0_iP_0.25_iS_1_iT_{seed}/'
        else:
            if kP in [0.1, 0.5]:
                folder_path = f'/Users/wangqianyu/Desktop/VertexModel/vertex_model/chiral_code/output/out_Nx_{Nx}_kA_0.5_kP_{kP:.1f}_P0_{P0:.1f}_disorder_0_tMAX_{tMAX}_dw_{dw}_di_{di}_ca1_0_ca2_0_iC_0_iP_0.25_iS_1_iT_{seed}/'
            else:
                folder_path = f'/Users/wangqianyu/Desktop/VertexModel/vertex_model/chiral_code/output/out_Nx_{Nx}_kA_0.5_kP_{kP:.2f}_P0_{P0:.1f}_disorder_0_tMAX_{tMAX}_dw_{dw}_di_{di}_ca1_0_ca2_0_iC_0_iP_0.25_iS_1_iT_{seed}/'

    original_dir = os.getcwd()

    try:
        os.chdir(folder_path)

        # Read the data using numpy, skipping comment lines starting with '#'
        # Data is space-separated
        data = np.loadtxt(file, comments='#', dtype=int)

        # Extract columns
        if data.ndim == 1:
            # Single row case
            Cid = np.array([data[0]])
            Cspin = np.array([data[1]])
        else:
            Cid = data[:, 0]
            Cspin = data[:, 1]

    finally:
        # Change back to original directory
        os.chdir('/Users/wangqianyu/Documents/MATLAB/UROP24')

    return Cspin
