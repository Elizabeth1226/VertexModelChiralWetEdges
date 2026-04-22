"""
Laplacian of velocity averaging function

Computes Laplacian of velocity components in cylindrical coordinates
"""

import numpy as np
import os
from read_file_folder_X import read_file_folder_X
from read_file_folder_cVertices import read_file_folder_cVertices
from read_file_folder_cell_bin_categories import read_file_folder_cell_bin_categories
from read_file_folder import read_file_folder
from numerical_differentiation import numerical_differentiation
from numerical_2_differentiation import numerical_2_differentiation

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

def LapV_ave_func(filepath, time, Lb):
    """
    Compute Laplacian of vertex velocity components binned by radial position

    In cylindrical coordinates:
    ∇²v_r = ∂²v_r/∂r² + (1/r)∂v_r/∂r - v_r/r²
    ∇²v_φ = ∂²v_φ/∂r² + (1/r)∂v_φ/∂r - v_φ/r²

    Parameters:
    -----------
    filepath : str
        Full path to the output folder
    time : int
        Current time step
    Lb : float
        Bin size multiplier

    Returns:
    --------
    cbins : ndarray
        Radial bin positions
    LapVr_list : ndarray
        Laplacian of radial velocity in each bin
    LapVphi_list : ndarray
        Laplacian of azimuthal velocity in each bin
    sigma_LapVr : ndarray
        Total error estimate for LapVr (statistical + truncation in quadrature)
    sigma_LapVphi : ndarray
        Total error estimate for LapVphi (statistical + truncation in quadrature)
    """

    data = read_v_xy_file_full(os.path.join(filepath, f'v_xy_{time}.dat'))
    print(f"Data shape: {data.shape}")
    print(f"Columns: {data.shape[1] if data.ndim > 1 else 1}")
    print(data[:,0])

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
    v_fr2 = np.zeros(nbins)    # sum of f_r^2 for variance estimation
    v_fphi2 = np.zeros(nbins)  # sum of f_phi^2 for variance estimation
    count = np.zeros(nbins)
    for v_id in range(Nv):
        f_r = (v_fx[v_id] * v_x_rel[v_id] + v_fy[v_id] * v_y_rel[v_id]) / np.sqrt(v_x_rel[v_id]**2 + v_y_rel[v_id]**2)
        f_phi = (v_fx[v_id] * (-v_y_rel[v_id]) + v_fy[v_id] * v_x_rel[v_id]) / np.sqrt(v_x_rel[v_id]**2 + v_y_rel[v_id]**2)
        v_fr[vCategories[v_id]] += f_r
        v_fphi[vCategories[v_id]] += f_phi
        v_fr2[vCategories[v_id]] += f_r**2
        v_fphi2[vCategories[v_id]] += f_phi**2
        count[vCategories[v_id]] += 1


    vr_list = np.divide(v_fr, count, out=np.zeros_like(v_fr), where=count != 0)
    vphi_list = np.divide(v_fphi, count, out=np.zeros_like(v_fphi), where=count != 0)

    # Standard error of the mean in each bin
    var_vr = np.maximum(
        np.divide(v_fr2, count, out=np.zeros_like(v_fr2), where=count > 0) - vr_list**2, 0.0)
    var_vphi = np.maximum(
        np.divide(v_fphi2, count, out=np.zeros_like(v_fphi2), where=count > 0) - vphi_list**2, 0.0)
    sigma_vr = np.sqrt(np.divide(var_vr, count, out=np.zeros_like(var_vr), where=count > 1))
    sigma_vphi = np.sqrt(np.divide(var_vphi, count, out=np.zeros_like(var_vphi), where=count > 1))


    # Compute Laplacian of radial velocity
    # ∇²v_r = ∂²v_r/∂r² + (1/r)∂v_r/∂r - v_r/r²
    divR_vr, _ = numerical_differentiation(0, vr_list, cbins)
    div2R_vr, _ = numerical_2_differentiation(0, vr_list, cbins)

    # Note: divR_vr has length Nb-1, div2R_vr has length Nb-2
    # We align them at indices 0:Nb-2 for divR_vr and all of div2R_vr
    LapVr_list = (divR_vr[:-1] / cbins[1:-1] + div2R_vr -
                  vr_list[1:-1] / cbins[1:-1]**2)

    # Compute Laplacian of azimuthal velocity
    # ∇²v_φ = ∂²v_φ/∂r² + (1/r)∂v_φ/∂r - v_φ/r²
    divR_vphi, _ = numerical_differentiation(0, vphi_list, cbins)
    div2R_vphi, _ = numerical_2_differentiation(0, vphi_list, cbins)

    # Average at two adjacent points for better smoothing
    LapVphi_list_1 = (divR_vphi[:-1] / cbins[1:-1] + div2R_vphi -
                      vphi_list[1:-1] / cbins[1:-1]**2)
    LapVphi_list_2 = (divR_vphi[1:] / cbins[2:] + div2R_vphi -
                      vphi_list[2:] / cbins[2:]**2)
    LapVphi_list = 0.5 * (LapVphi_list_1 + LapVphi_list_2)

    # ------------------------------------------------------------------ #
    # Error estimation                                                     #
    # ------------------------------------------------------------------ #
    h = cbins[1] - cbins[0]          # uniform bin step size
    r_out = cbins[1:-1]              # radial positions of Laplacian output
    r1    = cbins[1:-1]              # same as r_out (used for vphi)
    r2    = cbins[2:]                # r positions used by LapVphi_list_2

    # --- Statistical error for LapVr ---
    # LapVr[i] = A*vr[i] + B*vr[i+1] + C*vr[i+2], with:
    A_r = -1.0 / (h * r_out) + 1.0 / h**2
    B_r =  1.0 / (h * r_out) - 2.0 / h**2 - 1.0 / r_out**2
    C_r =  1.0 / h**2                        # scalar
    sigma_LapVr_stat = np.sqrt(
        A_r**2 * sigma_vr[:-2]**2 +
        B_r**2 * sigma_vr[1:-1]**2 +
        C_r**2 * sigma_vr[2:]**2)

    # --- Truncation error for LapVr ---
    # Forward-difference (O(h)) in divR_vr: error ≈ h/2 · |d²vr/dr²| / r
    eps_fwd_vr = 0.5 * h * np.abs(div2R_vr) / r_out
    # Central-difference (O(h²)) in div2R_vr: error ≈ h²/12 · |d⁴vr/dr⁴|
    # Estimate d⁴vr/dr⁴ via 5-point stencil; central point in vr_list is j = i+1
    N_vr = len(vr_list)
    d4_vr = np.zeros(len(div2R_vr))
    for i in range(len(div2R_vr)):
        j = i + 1
        if 2 <= j <= N_vr - 3:
            d4_vr[i] = (vr_list[j-2] - 4*vr_list[j-1] + 6*vr_list[j]
                        - 4*vr_list[j+1] + vr_list[j+2]) / h**4
    eps_cen_vr = (h**2 / 12.0) * np.abs(d4_vr)
    sigma_LapVr_trunc = np.sqrt(eps_fwd_vr**2 + eps_cen_vr**2)

    sigma_LapVr = np.sqrt(sigma_LapVr_stat**2 + sigma_LapVr_trunc**2)

    # --- Statistical error for LapVphi ---
    # LapVphi[i] = 0.5*(L1 + L2); combined coefficients of vphi[i..i+2]:
    Dp_i  = 0.5 * (-1.0/(h*r1) + 2.0/h**2)
    Dp_i1 = 0.5 * ( 1.0/(h*r1) - 1.0/(h*r2) - 4.0/h**2 - 1.0/r1**2)
    Dp_i2 = 0.5 * ( 2.0/h**2   + 1.0/(h*r2) - 1.0/r2**2)
    sigma_LapVphi_stat = np.sqrt(
        Dp_i**2  * sigma_vphi[:-2]**2 +
        Dp_i1**2 * sigma_vphi[1:-1]**2 +
        Dp_i2**2 * sigma_vphi[2:]**2)

    # --- Truncation error for LapVphi ---
    # L1 uses divR_vphi[i] at r1; L2 uses divR_vphi[i+1] at r2.
    # Approximate d²vphi/dr² at r2 with div2R_vphi[i+1] (fall back to [i] at boundary).
    div2R_vphi_next = np.empty_like(div2R_vphi)
    div2R_vphi_next[:-1] = div2R_vphi[1:]
    div2R_vphi_next[-1]  = div2R_vphi[-1]
    eps_fwd_vphi = 0.5 * (0.5*h * np.abs(div2R_vphi)      / r1 +
                           0.5*h * np.abs(div2R_vphi_next) / r2)
    N_vphi = len(vphi_list)
    d4_vphi = np.zeros(len(div2R_vphi))
    for i in range(len(div2R_vphi)):
        j = i + 1
        if 2 <= j <= N_vphi - 3:
            d4_vphi[i] = (vphi_list[j-2] - 4*vphi_list[j-1] + 6*vphi_list[j]
                          - 4*vphi_list[j+1] + vphi_list[j+2]) / h**4
    eps_cen_vphi = (h**2 / 12.0) * np.abs(d4_vphi)
    sigma_LapVphi_trunc = np.sqrt(eps_fwd_vphi**2 + eps_cen_vphi**2)

    sigma_LapVphi = np.sqrt(sigma_LapVphi_stat**2 + sigma_LapVphi_trunc**2)

    # Trim cbins to match output length
    cbins = cbins[1:-1]

    return cbins, LapVr_list, LapVphi_list, sigma_LapVr, sigma_LapVphi
