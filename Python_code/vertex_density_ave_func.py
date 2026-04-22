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
    
    
def vertex_density_ave_func(folder_path, time, Lb):
    filepath = os.path.join(folder_path, f'v_xy_{time}.dat')

    data = read_v_xy_file_full(filepath)
    print(data[:,0])

    v_fx = data [:,6]
    v_fy = data [:,7]
    v_x = data [:,2]
    v_y = data [:,3]

    Nv = len(v_fx)


    # Read data files
    c_v_xpos_list = read_file_folder_X(folder_path, f'X_{time}.dat')
    c_v_ypos_list = read_file_folder_X(folder_path, f'Y_{time}.dat')
    cVertices = read_file_folder_cVertices(folder_path, f'c_vertices_{time}.dat')
    Nc = len(cVertices)

    cActivity = read_file_folder_cell_bin_categories(folder_path, f'c_activity_{time}.dat')
    _, perioY = read_file_folder(folder_path, f'perio_{time}.dat')

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

    # Count vertices in each bin
    count = np.zeros(nbins)
    for v_id in range(Nv):
        count[vCategories[v_id]] += 1

    # Bin centers and annular ring areas
    bin_centers = 0.5 * (cbins[:-1] + cbins[1:])
    bin_areas = np.pi * (cbins[1:]**2 - cbins[:-1]**2)

    # Vertex density = count / area; uncertainty via Poisson statistics
    vertex_density = np.divide(count[:-1], bin_areas, out=np.zeros(nbins - 1), where=bin_areas != 0)
    uncertainty = np.divide(np.sqrt(count[:-1]), bin_areas, out=np.zeros(nbins - 1), where=bin_areas != 0)

    return vertex_density, bin_centers, uncertainty