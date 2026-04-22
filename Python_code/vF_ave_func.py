import os
import numpy as np
import matplotlib.pyplot as plt
import numpy as np
from read_file_folder_X import read_file_folder_X
from read_file_folder_cVertices import read_file_folder_cVertices
from read_file_folder_cell_bin_categories import read_file_folder_cell_bin_categories
from read_file_folder import read_file_folder



def read_v_xy_file_full(filepath):
    """
    Read all columns from v_xy_*.dat file
    
    Returns:
    --------
    data : numpy array
        Full data array with all columns
    """
    data = np.loadtxt(filepath, comments='#')
    return data

def vF_ave_func(filepath, time, Lb):
    """
    Calculate velocity magnitude average function for vertices

    parameters:
    -----------
    filepath : str
        Full path to the output folder
    time : int
        Time step
    Lb : float
        Bin size parameter (multiplier for bin width)

    Returns:
    --------
    cbins : numpy array
        Bin positions
    v_fr : numpy array
        radial velocity component average list
    v_fphi : numpy array
        azimuthal velocity component average list
    """

    data = read_v_xy_file_full(os.path.join(filepath, f'v_xy_{time}.dat'))
    # print(f"Data shape: {data.shape}")
    # print(f"Columns: {data.shape[1] if data.ndim > 1 else 1}")
    # print(data[:,0])

    v_fx = data [:,4]
    v_fy = data [:,5]
    v_x = data [:,2]
    v_y = data [:,3]

    Nv = len(v_fx)


    # Read data files
    c_v_xpos_list = read_file_folder_X(filepath, f'X_{time}.dat')
    c_v_ypos_list = read_file_folder_X(filepath, f'Y_{time}.dat')
    cVertices = read_file_folder_cVertices(filepath, f'c_vertices_{time}.dat')
    Nc = len(cVertices)

    cActivity = read_file_folder_cell_bin_categories(filepath, f'c_activity_{time}.dat')
    _, perioY = read_file_folder(filepath, f'perio_{time}.dat')

    # Calculate cell centers
    xpos_center_active = []
    ypos_center_active = []
    xpos_center = np.zeros(Nc)
    ypos_center = np.zeros(Nc)

    for cid in range(Nc):
        cVertices_row = cVertices[cid]
        c_v_xpos_list_current = c_v_xpos_list[cid]
        c_v_ypos_list_current = c_v_ypos_list[cid]
        Nv_current = int(cVertices_row[1])

        # Calculate cell center
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

    v_x_rel = v_x - chiral_cell_center[0]
    v_y_rel = v_y - chiral_cell_center[1]


    cbins = np.arange(0, perioY + Lb * perioY / np.sqrt(Nc), Lb * perioY / np.sqrt(Nc))
    nbins = len(cbins)

    # Bin vertices into categories based on distance from chiral cell center
    vCategories = np.zeros(Nv, dtype=int)
    for v_id in range(Nv):
        dist = np.sqrt(v_x_rel[v_id]**2 + v_y_rel[v_id]**2)
        for bin_id in range(nbins - 1):
            if cbins[bin_id] <= dist < cbins[bin_id + 1]:
                vCategories[v_id] = bin_id

    hist_vCategories, _ = np.histogram(vCategories, bins=cbins)
    # print('Vertex category histogram:', hist_vCategories)

    # plt.figure()
    # plt.hist(vCategories, bins=cbins)

    v_fr = np.zeros(nbins)
    v_fphi = np.zeros(nbins)
    count = np.zeros(nbins)
    for v_id in range(Nv):
        f_r = (v_fx[v_id] * v_x_rel[v_id] + v_fy[v_id] * v_y_rel[v_id]) / np.sqrt(v_x_rel[v_id]**2 + v_y_rel[v_id]**2)
        f_phi = (v_fx[v_id] * (-v_y_rel[v_id]) + v_fy[v_id] * v_x_rel[v_id]) / np.sqrt(v_x_rel[v_id]**2 + v_y_rel[v_id]**2)
        v_fr[vCategories[v_id]] += f_r
        v_fphi[vCategories[v_id]] += f_phi
        count[vCategories[v_id]] += 1


    v_fr = np.divide(v_fr, count, out=np.zeros_like(v_fr), where=count != 0)
    v_fphi = np.divide(v_fphi, count, out=np.zeros_like(v_fphi), where=count != 0)

    return cbins, v_fr, v_fphi