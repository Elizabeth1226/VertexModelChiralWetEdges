"""
Laplacian of velocity time-averaged plotting function

Computes and plots time-averaged Laplacian of velocity components
"""

import numpy as np
import matplotlib.pyplot as plt
from LapV_ave_func import LapV_ave_func


def LapV_time_ave_plot_func(P0, C1, dw, start_time, i, j, k, plot_style, Lb, filepath):
    """
    Compute time-averaged Laplacian of velocity and create error bar plot

    Parameters:
    -----------
    P0 : float
        Target perimeter (used in plot title)
    C1 : float
        Chiral coupling strength (used in plot title)
    dw : int
        Write interval
    start_time : int
        Starting time step
    i : int
        Subplot row index
    j : int
        Subplot column index
    k : int
        Subplot position
    plot_style : list of bool
        [plot_LapVr, plot_LapVphi] - controls which components to plot
    Lb : float
        Bin size multiplier
    filepath : str
        Full path to the output folder

    Returns:
    --------
    mean_T_1 : ndarray
        Mean radial positions (bin centers)
    mean_V_1 : ndarray
        Mean Laplacian of radial velocity
    mean_Vphi_1 : ndarray
        Mean Laplacian of azimuthal velocity
    err_V_1 : ndarray
        Combined error for LapVr (temporal std + mean numerical error in quadrature)
    err_Vphi_1 : ndarray
        Combined error for LapVphi (temporal std + mean numerical error in quadrature)
    """
    # Compute Laplacian of velocity at different time points

    T30_1, V30_1, Vphi30_1, sVr30_1, sVphi30_1 = LapV_ave_func(filepath, start_time, Lb)
    T35_1, V35_1, Vphi35_1, sVr35_1, sVphi35_1 = LapV_ave_func(filepath, start_time + 1*dw, Lb)
    T40_1, V40_1, Vphi40_1, sVr40_1, sVphi40_1 = LapV_ave_func(filepath, start_time + 2*dw, Lb)
    T45_1, V45_1, Vphi45_1, sVr45_1, sVphi45_1 = LapV_ave_func(filepath, start_time + 3*dw, Lb)
    T50_1, V50_1, Vphi50_1, sVr50_1, sVphi50_1 = LapV_ave_func(filepath, start_time + 4*dw, Lb)

    # Transpose position arrays (MATLAB returns column vectors, we want row vectors for consistency)
    T30_1 = T30_1.reshape(-1, 1)
    T35_1 = T35_1.reshape(-1, 1)
    T40_1 = T40_1.reshape(-1, 1)
    T45_1 = T45_1.reshape(-1, 1)
    T50_1 = T50_1.reshape(-1, 1)

    # Find minimum size across all time points (subtract 4 for safety margin)
    S = np.array([len(T30_1), len(T35_1), len(T40_1), len(T45_1), len(T50_1)])
    D = np.min(S) - 4

    # Create matrices with truncated data
    V_mat_1 = np.column_stack([
        V30_1[:D],
        V35_1[:D],
        V40_1[:D],
        V45_1[:D],
        V50_1[:D]
    ])

    Vphi_mat_1 = np.column_stack([
        Vphi30_1[:D],
        Vphi35_1[:D],
        Vphi40_1[:D],
        Vphi45_1[:D],
        Vphi50_1[:D]
    ])

    T_mat_1 = np.column_stack([
        T30_1[:D, 0],
        T35_1[:D, 0],
        T40_1[:D, 0],
        T45_1[:D, 0],
        T50_1[:D, 0]
    ])

    # Sigma matrices (numerical error per snapshot)
    sVr_mat_1 = np.column_stack([
        sVr30_1[:D], sVr35_1[:D], sVr40_1[:D], sVr45_1[:D], sVr50_1[:D]])
    sVphi_mat_1 = np.column_stack([
        sVphi30_1[:D], sVphi35_1[:D], sVphi40_1[:D], sVphi45_1[:D], sVphi50_1[:D]])

    # Compute mean and standard deviation across time points
    mean_T_1 = np.mean(T_mat_1, axis=1)
    mean_V_1 = np.mean(V_mat_1, axis=1)
    mean_Vphi_1 = np.mean(Vphi_mat_1, axis=1)
    std_V_1 = np.std(V_mat_1, axis=1)
    std_Vphi_1 = np.std(Vphi_mat_1, axis=1)

    # Combined error: temporal std and mean numerical error added in quadrature
    mean_sVr_1 = np.mean(sVr_mat_1, axis=1)
    mean_sVphi_1 = np.mean(sVphi_mat_1, axis=1)
    err_V_1    = np.sqrt(std_V_1**2    + mean_sVr_1**2)
    err_Vphi_1 = np.sqrt(std_Vphi_1**2 + mean_sVphi_1**2)

    # Create subplot
    plt.subplot(i, j, k)

    # Plot selected components
    if plot_style[0]:  # Plot radial component
        plt.errorbar(mean_T_1, mean_V_1, yerr=err_V_1, fmt='-s',
                     linewidth=1, markersize=5, label=r'$\nabla^2 V_r$')

    if plot_style[1]:  # Plot azimuthal component
        plt.errorbar(mean_T_1, mean_Vphi_1, yerr=err_Vphi_1, fmt='-s',
                     linewidth=1, markersize=5, label=r'$\nabla^2 V_\phi$')

    # Set axis limits
    plt.xlim([0, 30])
    plt.ylim([-0.01, 0.01])
    plt.ylim(auto=True)  # Auto y-axis

    # Set labels and title
    plt.title(f'$P_0={P0:.1f}$, $c_a^1={C1:.2f}$')
    plt.ylabel('cell radial velocity', fontsize=10)
    plt.grid(True)

    # Set font size for tick labels
    ax = plt.gca()
    ax.tick_params(labelsize=10)

    return mean_T_1, mean_V_1, mean_Vphi_1, err_V_1, err_Vphi_1
