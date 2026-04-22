import os
import numpy as np
from read_file_folder_X import read_file_folder_X
from read_file_folder_cVertices import read_file_folder_cVertices
from read_file_folder_cell_bin_categories import read_file_folder_cell_bin_categories
from read_file_folder import read_file_folder


def Vr_ave_func(filepath, time):
    """
    Calculate radial velocity average function

    Parameters:
    -----------
    filepath : str
        Full path to the output folder
    time : int
        Time step

    Returns:
    --------
    cbins : numpy array
        Bin positions
    vr_list : numpy array
        Radial velocity list
    """
    c_v_xpos_list = read_file_folder_X(filepath, f'X_{time}.dat')
    c_v_ypos_list = read_file_folder_X(filepath, f'Y_{time}.dat')
    cVertices = read_file_folder_cVertices(filepath, f'c_vertices_{time}.dat')

    Nc = len(cVertices)

    original_dir = os.getcwd()

    try:
        os.chdir(filepath)
        data = np.loadtxt(f"c_vel_{time}.dat")
        cVelocities = data[:, 1:3]
    finally:
        os.chdir(original_dir)

    cActivity = read_file_folder_cell_bin_categories(filepath, f'c_activity_{time}.dat')
    _, perioY = read_file_folder(filepath, f'perio_{time}.dat')

    cbins = np.arange(0, perioY + 1 * perioY / np.sqrt(Nc), 1 * perioY / np.sqrt(Nc))
    Nb = len(cbins)

    xpos_center_active = []
    ypos_center_active = []
    xpos_center = np.zeros(Nc)
    ypos_center = np.zeros(Nc)

    for cid in range(Nc):
        cVertices_row = cVertices[cid]
        c_v_xpos_list_current = c_v_xpos_list[cid]
        c_v_ypos_list_current = c_v_ypos_list[cid]
        Nv_current = cVertices_row[1]
        vid_count = 0

        for vid in cVertices_row[2:2 + Nv_current]:
            xpos_center[cid] += c_v_xpos_list_current[vid_count]
            ypos_center[cid] += c_v_ypos_list_current[vid_count]
            vid_count += 1

        xpos_center[cid] /= Nv_current
        ypos_center[cid] /= Nv_current

        if cActivity[cid] == 1:
            xpos_center_active.append(xpos_center[cid])
            ypos_center_active.append(ypos_center[cid])

    chiral_cell_center = np.zeros(2)
    chiral_cell_center[0] = np.sum(xpos_center_active) / len(xpos_center_active)
    chiral_cell_center[1] = np.sum(ypos_center_active) / len(ypos_center_active)

    xpos_center = xpos_center - chiral_cell_center[0]
    ypos_center = ypos_center - chiral_cell_center[1]

    cCategories = np.zeros(Nc, dtype=int)
    for cid in range(Nc):
        dist = np.sqrt(xpos_center[cid]**2 + ypos_center[cid]**2)
        for bin_id in range(len(cbins) - 1):
            if cbins[bin_id] <= dist < cbins[bin_id + 1]:
                cCategories[cid] = bin_id

    theta_list = np.zeros(Nc)
    for cid in range(Nc):
        theta = np.arctan2(ypos_center[cid], xpos_center[cid])
        theta_list[cid] = theta  # Fixed: was theta_list[Nc] in MATLAB

    vr_list = np.zeros(Nb)
    bin_count = np.zeros(Nb)

    for cid in range(Nc):
        vr = cVelocities[cid, 0] * np.cos(theta_list[cid]) + cVelocities[cid, 1] * np.sin(theta_list[cid])
        vr_list[cCategories[cid]] += vr
        bin_count[cCategories[cid]] += 1

    for bin_id in range(Nb):
        if bin_count[bin_id] > 0:
            vr_list[bin_id] /= bin_count[bin_id]

    return cbins, vr_list
