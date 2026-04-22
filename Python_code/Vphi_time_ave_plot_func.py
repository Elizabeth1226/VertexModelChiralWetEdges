"""
Azimuthal velocity time-averaged plotting function

Computes and plots time-averaged azimuthal (tangential) velocity profiles
"""

import numpy as np
import matplotlib.pyplot as plt
from Vphi_ave_func import Vphi_ave_func


def Vphi_time_ave_plot_func(P0, C1, dw, start_time, i, j, k, filepath):
    """
    Compute time-averaged azimuthal velocity and create error bar plot

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
    filepath : str
        Full path to the output folder

    Returns:
    --------
    mean_T_1 : ndarray
        Mean radial positions (bin centers)
    mean_V_1 : ndarray
        Mean azimuthal velocities
    """
    # Compute azimuthal velocity at different time points
    T30_1, V30_1 = Vphi_ave_func(filepath, start_time)
    T35_1, V35_1 = Vphi_ave_func(filepath, start_time + 1*dw)
    T40_1, V40_1 = Vphi_ave_func(filepath, start_time + 2*dw)
    T45_1, V45_1 = Vphi_ave_func(filepath, start_time + 3*dw)
    T50_1, V50_1 = Vphi_ave_func(filepath, start_time + 4*dw)

    # Transpose position arrays (MATLAB returns column vectors, we want row vectors for consistency)
    T30_1 = T30_1.reshape(-1, 1)
    T35_1 = T35_1.reshape(-1, 1)
    T40_1 = T40_1.reshape(-1, 1)
    T45_1 = T45_1.reshape(-1, 1)
    T50_1 = T50_1.reshape(-1, 1)

    # Find minimum size across all time points
    S = np.array([len(T30_1), len(T35_1), len(T40_1), len(T45_1), len(T50_1)])
    D = np.min(S)

    # Create matrices with truncated data
    V_mat_1 = np.column_stack([
        V30_1[:D],
        V35_1[:D],
        V40_1[:D],
        V45_1[:D],
        V50_1[:D]
    ])

    T_mat_1 = np.column_stack([
        T30_1[:D, 0],
        T35_1[:D, 0],
        T40_1[:D, 0],
        T45_1[:D, 0],
        T50_1[:D, 0]
    ])

    # Compute mean and standard deviation across time points
    mean_T_1 = np.mean(T_mat_1, axis=1)
    mean_V_1 = np.mean(V_mat_1, axis=1)
    std_V_1 = np.std(V_mat_1, axis=1)

    # Create subplot
    plt.subplot(i, j, k)
    plt.errorbar(mean_T_1, mean_V_1, yerr=std_V_1, fmt='-s',
                 linewidth=1, markersize=5,
                 markeredgecolor='blue', markerfacecolor=[0.65, 0.85, 0.90])

    # Set axis limits
    plt.xlim([0, 30])
    plt.ylim(auto=True)  # Auto y-axis

    # Set labels and title
    plt.title(f'$P_0={P0:.1f}$, $c_a^1={C1:.2f}$')
    plt.ylabel('tangential velocity', fontsize=5)
    plt.grid(True)

    # Set font size for tick labels
    ax = plt.gca()
    ax.tick_params(labelsize=5)

    return mean_T_1, mean_V_1
