import numpy as np
import matplotlib.pyplot as plt
from read_file_folder import read_file_folder


def V_time_ave_plot_func(Nx, P0, C1, dw, di, disorder, tMAX, seed, i, j, k):
    """
    Plot time-averaged velocity

    Parameters:
    -----------
    Nx, P0, C1, dw, di, disorder, tMAX, seed : numeric
        Parameters for folder path construction
    i, j, k : int
        Subplot position parameters

    Returns:
    --------
    None (creates plot)
    """
    # Read data from multiple time steps
    T5_1, V5_1 = read_file_folder(Nx, P0, C1, dw, di, disorder, tMAX, seed, 'azimuthal_ave_v_45500.dat')
    T10_1, V10_1 = read_file_folder(Nx, P0, C1, dw, di, disorder, tMAX, seed, 'azimuthal_ave_v_46000.dat')
    T15_1, V15_1 = read_file_folder(Nx, P0, C1, dw, di, disorder, tMAX, seed, 'azimuthal_ave_v_46500.dat')
    T20_1, V20_1 = read_file_folder(Nx, P0, C1, dw, di, disorder, tMAX, seed, 'azimuthal_ave_v_47000.dat')
    T25_1, V25_1 = read_file_folder(Nx, P0, C1, dw, di, disorder, tMAX, seed, 'azimuthal_ave_v_47500.dat')
    T30_1, V30_1 = read_file_folder(Nx, P0, C1, dw, di, disorder, tMAX, seed, 'azimuthal_ave_v_48000.dat')
    T35_1, V35_1 = read_file_folder(Nx, P0, C1, dw, di, disorder, tMAX, seed, 'azimuthal_ave_v_48500.dat')
    T40_1, V40_1 = read_file_folder(Nx, P0, C1, dw, di, disorder, tMAX, seed, 'azimuthal_ave_v_49000.dat')
    T45_1, V45_1 = read_file_folder(Nx, P0, C1, dw, di, disorder, tMAX, seed, 'azimuthal_ave_v_49500.dat')
    T50_1, V50_1 = read_file_folder(Nx, P0, C1, dw, di, disorder, tMAX, seed, 'azimuthal_ave_v_50000.dat')

    # Find minimum size
    S = np.array([len(T5_1), len(T10_1), len(T15_1), len(T20_1), len(T25_1),
                  len(T30_1), len(T35_1), len(T40_1), len(T45_1), len(T50_1)])
    D = np.min(S)

    # Create matrices
    V_mat_1 = np.column_stack([V5_1[:D], V10_1[:D], V15_1[:D], V20_1[:D], V25_1[:D],
                                V30_1[:D], V35_1[:D], V40_1[:D], V45_1[:D], V50_1[:D]])
    std_V_1 = np.std(V_mat_1, axis=1)

    T_mat_1 = np.column_stack([T5_1[:D], T10_1[:D], T15_1[:D], T20_1[:D], T25_1[:D],
                                T30_1[:D], T35_1[:D], T40_1[:D], T45_1[:D], T50_1[:D]])
    mean_T_1 = np.mean(T_mat_1, axis=1)
    mean_V_1 = np.mean(V_mat_1, axis=1)

    # Create subplot
    plt.subplot(i, j, k)
    plt.errorbar(mean_T_1, mean_V_1, yerr=std_V_1, fmt="-s", linewidth=1, markersize=5)
    plt.xlim([0, 30])
    plt.ylim([0, 0.01])
    plt.autoscale(enable=True, axis='y')
    plt.title(f'P0={P0:.1f}, ca1={C1:.2f}')
    plt.xlabel('cell position', fontsize=10)
    plt.ylabel('cell speed', fontsize=10)
    plt.gca().tick_params(labelsize=10)
    plt.grid(True)
