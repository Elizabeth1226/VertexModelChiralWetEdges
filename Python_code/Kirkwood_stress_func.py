import numpy as np
from read_file_folder_X import read_file_folder_X
from read_file_folder_cVertices import read_file_folder_cVertices
from read_file_folder_vel import read_file_folder_vel
from read_file_folder_v_vel import read_file_folder_v_vel
from read_file_folder_cell_bin_categories import read_file_folder_cell_bin_categories
from read_file_folder import read_file_folder
import os
cwd = os.getcwd()

# NOTE: These helper functions are not in the provided MATLAB files
# You'll need to implement them separately
# from helper_functions import get_rotation_angle, tensor_coordinate_transformation, get_Q_angle


def get_rotation_angle(x, y):
    """Helper function to calculate rotation angle"""
    return np.arctan2(y, x)


def tensor_coordinate_transformation(tensor, theta):
    """Helper function to perform tensor coordinate transformation"""
    rotation_matrix = np.array([[np.cos(theta), np.sin(theta)],
                                 [-np.sin(theta), np.cos(theta)]])
    return rotation_matrix @ tensor @ rotation_matrix.T


def get_Q_angle(Q):
    """Helper function to calculate Q tensor angle"""
    # Calculate angle of Q tensor (typically the orientation angle)
    return 0.5 * np.arctan2(2 * Q[0, 1], Q[0, 0] - Q[1, 1])


def Kirkwood_stress_func(kP, P0, time, force_type, rho_0, folder_path, Lb):
    """
    Calculate the Kirkwood stress of all cells at given time

    Parameters:
    -----------
    folder_path : str
        Path to the data folder
    kP : float
        Perimeter stiffness
    P0 : float
        Target perimeter
    time : int
        Time step
    force_type : int
        1: cell force; 2: cell elastic force; 3: cell perimeter force;
        4: cell area force; 5: cell chiral force
    rho_0 : float
        Density parameter
    Lb : float
        Bin size parameter (multiplier for bin width)

    Returns:
    --------
    Multiple numpy arrays containing stress components for total, active, and passive cells

    NOTE: Q_theoretical is only defined for force_type = cell perimeter force (3)
    """

    rho_0 = 1

    # Read data files
    c_v_xpos_list = read_file_folder_X(folder_path, f'X_{time}.dat')
    c_v_ypos_list = read_file_folder_X(folder_path, f'Y_{time}.dat')
    cVertices = read_file_folder_cVertices(folder_path, f'c_vertices_{time}.dat')
    Nc = len(cVertices)

    if force_type == 1:
        cForces = read_file_folder_cVertices(folder_path, f'cell_force_{time}.dat')
    elif force_type == 5:
        cForces = read_file_folder_cVertices(folder_path, f'cell_chiral_force_{time}.dat')
    elif force_type == 2:
        cForces = read_file_folder_cVertices(folder_path, f'cell_elastic_force_{time}.dat')
    elif force_type == 3:
        cForces = read_file_folder_cVertices(folder_path, f'cell_perimeter_force_{time}.dat')
    elif force_type == 4:
        cForces = read_file_folder_cVertices(folder_path, f'cell_area_force_{time}.dat')


    cActivity = read_file_folder_cell_bin_categories(folder_path, f'c_activity_{time}.dat')
    _, perioY = read_file_folder(folder_path, f'perio_{time}.dat')

    cbins = np.arange(0, perioY + Lb * perioY / np.sqrt(Nc), Lb * perioY / np.sqrt(Nc))
    nbins = len(cbins)

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

    xpos_center = xpos_center - chiral_cell_center[0]
    ypos_center = ypos_center - chiral_cell_center[1]

    # Bin cells
    cCategories = np.zeros(Nc, dtype=int)
    for cid in range(Nc):
        dist = np.sqrt(xpos_center[cid]**2 + ypos_center[cid]**2)
        for bin_id in range(nbins - 1):
            if cbins[bin_id] <= dist < cbins[bin_id + 1]:
                cCategories[cid] = bin_id

    # Initialize bin arrays
    trace_bin = np.zeros(nbins)
    trace_theoretical_bin = np.zeros(nbins)
    trace_bin_active = np.zeros(nbins)
    trace_theoretical_bin_active = np.zeros(nbins)
    antisymmetric_part_bin = np.zeros((2 * nbins, 2))
    antisymmetric_part_bin_active = np.zeros((2 * nbins, 2))
    shear_part_bin = np.zeros((2 * nbins, 2))
    shear_part_bin_active = np.zeros((2 * nbins, 2))
    trace_bin_passive = np.zeros(nbins)
    trace_theoretical_bin_passive = np.zeros(nbins)
    antisymmetric_part_bin_passive = np.zeros((2 * nbins, 2))
    shear_part_bin_passive = np.zeros((2 * nbins, 2))
    Q_theoretical_bin = np.zeros((2 * nbins, 2))
    Q_bin = np.zeros((2 * nbins, 2))
    bin_count = np.zeros(nbins)
    bin_count_active = np.zeros(nbins)
    bin_count_passive = np.zeros(nbins)

    # Calculate stress for each cell
    for cid in range(Nc):
        cVertices_row = cVertices[cid]
        c_v_xpos_list_current = np.array(c_v_xpos_list[cid]) - xpos_center[cid] - chiral_cell_center[0]
        c_v_ypos_list_current = np.array(c_v_ypos_list[cid]) - ypos_center[cid] - chiral_cell_center[1]
        Nv_current = int(cVertices_row[1])
        cForces_row = cForces[cid]
        
        # Extract forces (cForces_row: list of [Nv, Fx1, Fy1, Fx2, Fy2, ..., FxNv, FyNv])
        vFx_list_current = np.asarray(cForces_row[1::2])
        vFy_list_current = np.asarray(cForces_row[2::2])


        # Calculate area and perimeter
        area = 0
        perimeter = 0
        for i in range(Nv_current - 1):
            area += abs(0.5 * (c_v_xpos_list_current[i + 1] * c_v_ypos_list_current[i] -
                               c_v_xpos_list_current[i] * c_v_ypos_list_current[i + 1]))
            perimeter += np.sqrt((c_v_xpos_list_current[i + 1] - c_v_xpos_list_current[i])**2 +
                                 (c_v_ypos_list_current[i + 1] - c_v_ypos_list_current[i])**2)

        area += abs(0.5 * (c_v_xpos_list_current[0] * c_v_ypos_list_current[Nv_current - 1] -
                           c_v_xpos_list_current[Nv_current - 1] * c_v_ypos_list_current[0]))
        perimeter += np.sqrt((c_v_xpos_list_current[0] - c_v_xpos_list_current[Nv_current - 1])**2 +
                             (c_v_ypos_list_current[0] - c_v_ypos_list_current[Nv_current - 1])**2)

        # Calculate Q tensor
        l_current = np.zeros(2)
        Q = np.zeros((2, 2))
        for i in range(Nv_current - 1):
            l_current[0] = c_v_xpos_list_current[i + 1] - c_v_xpos_list_current[i]
            l_current[1] = c_v_ypos_list_current[i + 1] - c_v_ypos_list_current[i]
            edgeLength = np.sqrt(l_current @ l_current)
            Q += np.outer(l_current, l_current) / edgeLength / perimeter

        l_current[0] = c_v_xpos_list_current[0] - c_v_xpos_list_current[Nv_current - 1]
        l_current[1] = c_v_ypos_list_current[0] - c_v_ypos_list_current[Nv_current - 1]
        edgeLength = np.sqrt(l_current @ l_current)
        Q += np.outer(l_current, l_current) / edgeLength / perimeter - 0.5 * np.eye(2)

        
        # Calculate stress tensor
        r_current = np.zeros(2)
        v_current = np.zeros(2)
        outerProduct_current = np.zeros((2, 2))

        for i in range(Nv_current):
            r_current[0] = c_v_xpos_list_current[i]
            r_current[1] = c_v_ypos_list_current[i]
            v_current[0] = vFx_list_current[i]
            v_current[1] = vFy_list_current[i]
            outerProduct_current += np.outer(r_current, v_current)

        outerProduct_current = -rho_0 * outerProduct_current / area
        theta = get_rotation_angle(xpos_center[cid], ypos_center[cid])
        sigma_current = tensor_coordinate_transformation(outerProduct_current, theta)
        isotropic_part_current = 0.5 * np.trace(sigma_current) * np.eye(2)
        antisymmetric_part_current = 0.5 * (sigma_current - sigma_current.T)
        shear_part_current = 0.5 * (sigma_current + sigma_current.T - 2 * isotropic_part_current)
        Q_theoretical = area / (2 * kP) / perimeter / (perimeter - P0) * sigma_current - 0.5 * np.eye(2)
        Q = tensor_coordinate_transformation(Q, theta)


        # Store in appropriate bins
        if cActivity[cid] > 0.0:
            trace_theoretical_bin_active[cCategories[cid]] += 2 * (area - 1) + 2 * kP / area * (perimeter - P0) * perimeter
            trace_bin_active[cCategories[cid]] += np.trace(sigma_current)
            antisymmetric_part_bin_active[2 * cCategories[cid]:2 * cCategories[cid] + 2, :] += antisymmetric_part_current
            shear_part_bin_active[2 * cCategories[cid]:2 * cCategories[cid] + 2, :] += shear_part_current
            bin_count_active[cCategories[cid]] += 1
        else:
            trace_theoretical_bin_passive[cCategories[cid]] += 2 * (area - 1) + 2 * kP / area * (perimeter - P0) * perimeter
            trace_bin_passive[cCategories[cid]] += np.trace(sigma_current)
            antisymmetric_part_bin_passive[2 * cCategories[cid]:2 * cCategories[cid] + 2, :] += antisymmetric_part_current
            shear_part_bin_passive[2 * cCategories[cid]:2 * cCategories[cid] + 2, :] += shear_part_current
            bin_count_passive[cCategories[cid]] += 1

        trace_theoretical_bin[cCategories[cid]] += 2 * (area - 1) + 2 * kP / area * (perimeter - P0) * perimeter
        trace_bin[cCategories[cid]] += np.trace(sigma_current)
        antisymmetric_part_bin[2 * cCategories[cid]:2 * cCategories[cid] + 2, :] += antisymmetric_part_current
        shear_part_bin[2 * cCategories[cid]:2 * cCategories[cid] + 2, :] += shear_part_current
        bin_count[cCategories[cid]] += 1
        Q_theoretical_bin[2 * cCategories[cid]:2 * cCategories[cid] + 2, :] += Q_theoretical
        Q_bin[2 * cCategories[cid]:2 * cCategories[cid] + 2, :] += Q

    # Average over bins
    trace_total = []
    trace_active = []
    trace_passive = []
    trace_theoretical_total = []
    trace_theoretical_active = []
    trace_theoretical_passive = []
    antisymmetric_xy_total = []
    antisymmetric_xy_active = []
    antisymmetric_xy_passive = []
    shear_xx_total = []
    shear_xy_total = []
    shear_xx_active = []
    shear_xy_active = []
    shear_xx_passive = []
    shear_xy_passive = []
    shear_det_active = []
    shear_det_passive = []
    shear_angle_active = []
    shear_angle_passive = []
    Q_angle_theoretical = []
    Q_det_theoretical = []
    Q_angle = []
    Q_det = []
    bin_pos_total = []
    bin_pos_active = []
    bin_pos_passive = []

    for bin_id in range(len(bin_count_active)):
        if bin_count_active[bin_id] > 0:
            trace_theoretical_active.append(trace_theoretical_bin_active[bin_id] / bin_count_active[bin_id])
            trace_active.append(trace_bin_active[bin_id] / bin_count_active[bin_id])
            antisymmetric_xy_active.append(antisymmetric_part_bin_active[2 * bin_id, 1] / bin_count_active[bin_id])
            shear_xx_active.append(shear_part_bin_active[2 * bin_id, 0] / bin_count_active[bin_id])
            shear_xy_active.append(shear_part_bin_active[2 * bin_id, 1] / bin_count_active[bin_id])
            shear_det_active.append(np.linalg.det(shear_part_bin_active[2 * bin_id:2 * bin_id + 2, :] / bin_count_active[bin_id]))
            shear_angle_active.append(get_Q_angle(shear_part_bin_active[2 * bin_id:2 * bin_id + 2, :] / bin_count_active[bin_id]))
            bin_pos_active.append(cbins[bin_id])

    for bin_id in range(len(bin_count_passive)):
        if bin_count_passive[bin_id] > 0:
            trace_theoretical_passive.append(trace_theoretical_bin_passive[bin_id] / bin_count_passive[bin_id])
            trace_passive.append(trace_bin_passive[bin_id] / bin_count_passive[bin_id])
            antisymmetric_xy_passive.append(antisymmetric_part_bin_passive[2 * bin_id, 1] / bin_count_passive[bin_id])
            shear_xx_passive.append(shear_part_bin_passive[2 * bin_id, 0] / bin_count_passive[bin_id])
            shear_xy_passive.append(shear_part_bin_passive[2 * bin_id, 1] / bin_count_passive[bin_id])
            shear_det_passive.append(np.linalg.det(shear_part_bin_passive[2 * bin_id:2 * bin_id + 2, :] / bin_count_passive[bin_id]))
            shear_angle_passive.append(get_Q_angle(shear_part_bin_passive[2 * bin_id:2 * bin_id + 2, :] / bin_count_passive[bin_id]))
            bin_pos_passive.append(cbins[bin_id])

    for bin_id in range(len(bin_count)):
        if bin_count[bin_id] > 0:
            trace_theoretical_total.append(trace_theoretical_bin[bin_id] / bin_count[bin_id])
            trace_total.append(trace_bin[bin_id] / bin_count[bin_id])
            antisymmetric_xy_total.append(antisymmetric_part_bin[2 * bin_id, 1] / bin_count[bin_id])
            shear_xx_total.append(shear_part_bin[2 * bin_id, 0] / bin_count[bin_id])
            shear_xy_total.append(shear_part_bin[2 * bin_id, 1] / bin_count[bin_id])
            bin_pos_total.append(cbins[bin_id])
            Q_det_theoretical.append(np.linalg.det(Q_theoretical_bin[2 * bin_id:2 * bin_id + 2, :] / bin_count[bin_id]))
            Q_angle_theoretical.append(get_Q_angle(Q_theoretical_bin[2 * bin_id:2 * bin_id + 2, :] / bin_count[bin_id]))
            Q_det.append(np.linalg.det(Q_bin[2 * bin_id:2 * bin_id + 2, :] / bin_count[bin_id]))
            Q_angle.append(get_Q_angle(Q_bin[2 * bin_id:2 * bin_id + 2, :] / bin_count[bin_id]))

    # Convert to numpy arrays
    return (np.array(trace_total), np.array(trace_active), np.array(trace_passive),
            np.array(trace_theoretical_total), np.array(trace_theoretical_active), np.array(trace_theoretical_passive),
            np.array(antisymmetric_xy_total), np.array(antisymmetric_xy_active), np.array(antisymmetric_xy_passive),
            np.array(shear_xx_total), np.array(shear_xx_active), np.array(shear_xx_passive),
            np.array(shear_xy_total), np.array(shear_xy_active), np.array(shear_xy_passive),
            np.array(shear_det_active), np.array(shear_det_passive),
            np.array(shear_angle_active), np.array(shear_angle_passive),
            np.array(Q_det_theoretical), np.array(Q_angle_theoretical),
            np.array(Q_det), np.array(Q_angle),
            np.array(bin_pos_total), np.array(bin_pos_active), np.array(bin_pos_passive))







