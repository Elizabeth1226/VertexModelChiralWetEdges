import numpy as np
import matplotlib.pyplot as plt
from Kirkwood_stress_func import Kirkwood_stress_func


def Kirkwood_stress_ave_time_plot_func(kP, P0, C1, dw, start_time, i, j, k, cell_type, force_type, rho_0, plot_style, color1, color2, filepath):
    """
    Plot time-averaged Kirkwood stress components

    Parameters:
    -----------
    kP : float
        Perimeter stiffness
    P0 : float
        Target perimeter
    C1 : float
        Chiral coupling strength (used in plot title)
    dw : int
        Write interval
    start_time : int
        Starting time step
    i, j, k : int
        Subplot position parameters
    cell_type : list of bool
        [plot_active_cells, plot_passive_cells]
    force_type : int
        Type of force (1-5)
    rho_0 : float
        Density parameter
    plot_style : list of bool
        [total_trace, antisym_xy, shear_xx_xy, shear_det, shear_angle,
         fit_trace, Q_angle_theoretical, Q_det_theoretical]
    color1, color2 : array-like
        Colors for active and passive cells
    filepath : str
        Full path to the output folder

    Returns:
    --------
    yData, xData : numpy arrays
        Y and X data for the plot
    """
    # Call Kirkwood_stress_func at multiple time steps
    results_5 = Kirkwood_stress_func(kP, P0, start_time, force_type, rho_0, filepath)
    results_10 = Kirkwood_stress_func(kP, P0, start_time + 1*dw, force_type, rho_0, filepath)
    results_15 = Kirkwood_stress_func(kP, P0, start_time + 2*dw, force_type, rho_0, filepath)
    results_20 = Kirkwood_stress_func(kP, P0, start_time + 3*dw, force_type, rho_0, filepath)
    results_25 = Kirkwood_stress_func(kP, P0, start_time + 4*dw, force_type, rho_0, filepath)

    # Unpack results
    (trace_total_5, _, _, trace_theoretical_total_5, _, _, _, antisym_xy_active_5, antisym_xy_passive_5,
     _, shear_xx_active_5, shear_xx_passive_5, _, shear_xy_active_5, shear_xy_passive_5,
     shear_det_active_5, shear_det_passive_5, shear_angle_active_5, shear_angle_passive_5,
     Q_det_theoretical_5, Q_angle_theoretical_5, _, _, bin_pos_total_5, bin_pos_active_5, bin_pos_passive_5) = results_5

    (trace_total_10, _, _, trace_theoretical_total_10, _, _, _, antisym_xy_active_10, antisym_xy_passive_10,
     _, shear_xx_active_10, shear_xx_passive_10, _, shear_xy_active_10, shear_xy_passive_10,
     shear_det_active_10, shear_det_passive_10, shear_angle_active_10, shear_angle_passive_10,
     Q_det_theoretical_10, Q_angle_theoretical_10, _, _, bin_pos_total_10, bin_pos_active_10, bin_pos_passive_10) = results_10

    (trace_total_15, _, _, trace_theoretical_total_15, _, _, _, antisym_xy_active_15, antisym_xy_passive_15,
     _, shear_xx_active_15, shear_xx_passive_15, _, shear_xy_active_15, shear_xy_passive_15,
     shear_det_active_15, shear_det_passive_15, shear_angle_active_15, shear_angle_passive_15,
     Q_det_theoretical_15, Q_angle_theoretical_15, _, _, bin_pos_total_15, bin_pos_active_15, bin_pos_passive_15) = results_15

    (trace_total_20, _, _, trace_theoretical_total_20, _, _, _, antisym_xy_active_20, antisym_xy_passive_20,
     _, shear_xx_active_20, shear_xx_passive_20, _, shear_xy_active_20, shear_xy_passive_20,
     shear_det_active_20, shear_det_passive_20, shear_angle_active_20, shear_angle_passive_20,
     Q_det_theoretical_20, Q_angle_theoretical_20, _, _, bin_pos_total_20, bin_pos_active_20, bin_pos_passive_20) = results_20

    (trace_total_25, _, _, trace_theoretical_total_25, _, _, _, antisym_xy_active_25, antisym_xy_passive_25,
     _, shear_xx_active_25, shear_xx_passive_25, _, shear_xy_active_25, shear_xy_passive_25,
     shear_det_active_25, shear_det_passive_25, shear_angle_active_25, shear_angle_passive_25,
     Q_det_theoretical_25, Q_angle_theoretical_25, _, _, bin_pos_total_25, bin_pos_active_25, bin_pos_passive_25) = results_25

    # Find minimum sizes for each category
    S_active = np.array([len(bin_pos_active_5), len(bin_pos_active_10), len(bin_pos_active_15),
                         len(bin_pos_active_20), len(bin_pos_active_25)])
    D1 = np.min(S_active)

    S_passive = np.array([len(bin_pos_passive_5), len(bin_pos_passive_10), len(bin_pos_passive_15),
                          len(bin_pos_passive_20), len(bin_pos_passive_25)])
    D2 = np.min(S_passive)

    S_total = np.array([len(bin_pos_total_5), len(bin_pos_total_10), len(bin_pos_total_15),
                        len(bin_pos_total_20), len(bin_pos_total_25)])
    D3 = np.min(S_total)

    # Create matrices for active cells
    antisym_xy_mat_active = np.column_stack([antisym_xy_active_5[:D1], antisym_xy_active_10[:D1],
                                              antisym_xy_active_15[:D1], antisym_xy_active_20[:D1],
                                              antisym_xy_active_25[:D1]])
    shear_xx_mat_active = np.column_stack([shear_xx_active_5[:D1], shear_xx_active_10[:D1],
                                            shear_xx_active_15[:D1], shear_xx_active_20[:D1],
                                            shear_xx_active_25[:D1]])
    shear_xy_mat_active = np.column_stack([shear_xy_active_5[:D1], shear_xy_active_10[:D1],
                                            shear_xy_active_15[:D1], shear_xy_active_20[:D1],
                                            shear_xy_active_25[:D1]])
    shear_det_mat_active = np.column_stack([shear_det_active_5[:D1], shear_det_active_10[:D1],
                                             shear_det_active_15[:D1], shear_det_active_20[:D1],
                                             shear_det_active_25[:D1]])
    shear_angle_mat_active = np.column_stack([shear_angle_active_5[:D1], shear_angle_active_10[:D1],
                                               shear_angle_active_15[:D1], shear_angle_active_20[:D1],
                                               shear_angle_active_25[:D1]])
    bin_pos_mat_active = np.column_stack([bin_pos_active_5[:D1], bin_pos_active_10[:D1],
                                           bin_pos_active_15[:D1], bin_pos_active_20[:D1],
                                           bin_pos_active_25[:D1]])

    # Create matrices for passive cells
    antisym_xy_mat_passive = np.column_stack([antisym_xy_passive_5[1:D2], antisym_xy_passive_10[1:D2],
                                               antisym_xy_passive_15[1:D2], antisym_xy_passive_20[1:D2],
                                               antisym_xy_passive_25[1:D2]])
    shear_xx_mat_passive = np.column_stack([shear_xx_passive_5[1:D2], shear_xx_passive_10[1:D2],
                                             shear_xx_passive_15[1:D2], shear_xx_passive_20[1:D2],
                                             shear_xx_passive_25[1:D2]])
    shear_xy_mat_passive = np.column_stack([shear_xy_passive_5[1:D2], shear_xy_passive_10[1:D2],
                                             shear_xy_passive_15[1:D2], shear_xy_passive_20[1:D2],
                                             shear_xy_passive_25[1:D2]])
    shear_det_mat_passive = np.column_stack([shear_det_passive_5[1:D2], shear_det_passive_10[1:D2],
                                              shear_det_passive_15[1:D2], shear_det_passive_20[1:D2],
                                              shear_det_passive_25[1:D2]])
    shear_angle_mat_passive = np.column_stack([shear_angle_passive_5[1:D2], shear_angle_passive_10[1:D2],
                                                shear_angle_passive_15[1:D2], shear_angle_passive_20[1:D2],
                                                shear_angle_passive_25[1:D2]])
    bin_pos_mat_passive = np.column_stack([bin_pos_passive_5[1:D2], bin_pos_passive_10[1:D2],
                                            bin_pos_passive_15[1:D2], bin_pos_passive_20[1:D2],
                                            bin_pos_passive_25[1:D2]])

    # Create matrices for total
    Q_angle_theoretical_mat = np.column_stack([Q_angle_theoretical_5[:D3], Q_angle_theoretical_10[:D3],
                                                Q_angle_theoretical_15[:D3], Q_angle_theoretical_20[:D3],
                                                Q_angle_theoretical_25[:D3]])
    Q_det_theoretical_mat = np.column_stack([Q_det_theoretical_5[:D3], Q_det_theoretical_10[:D3],
                                              Q_det_theoretical_15[:D3], Q_det_theoretical_20[:D3],
                                              Q_det_theoretical_25[:D3]])
    trace_total_mat = np.column_stack([trace_total_5[:D3], trace_total_10[:D3],
                                        trace_total_15[:D3], trace_total_20[:D3],
                                        trace_total_25[:D3]])
    trace_theoretical_total_mat = np.column_stack([trace_theoretical_total_5[:D3], trace_theoretical_total_10[:D3],
                                                    trace_theoretical_total_15[:D3], trace_theoretical_total_20[:D3],
                                                    trace_theoretical_total_25[:D3]])
    bin_pos_total_mat = np.column_stack([bin_pos_total_5[:D3], bin_pos_total_10[:D3],
                                          bin_pos_total_15[:D3], bin_pos_total_20[:D3],
                                          bin_pos_total_25[:D3]])

    # Calculate means and standard deviations
    std_trace_total = np.std(trace_total_mat, axis=1)
    std_antisym_xy_active = np.std(antisym_xy_mat_active, axis=1)
    std_shear_xx_active = np.std(shear_xx_mat_active, axis=1)
    std_shear_xy_active = np.std(shear_xy_mat_active, axis=1)
    std_shear_det_active = np.std(shear_det_mat_active, axis=1)
    std_shear_angle_active = np.std(shear_angle_mat_active, axis=1)
    std_trace_theoretical_total = np.std(trace_theoretical_total_mat, axis=1)
    std_antisym_xy_passive = np.std(antisym_xy_mat_passive, axis=1)
    std_shear_xx_passive = np.std(shear_xx_mat_passive, axis=1)
    std_shear_xy_passive = np.std(shear_xy_mat_passive, axis=1)
    std_shear_det_passive = np.std(shear_det_mat_passive, axis=1)
    std_shear_angle_passive = np.std(shear_angle_mat_passive, axis=1)
    std_Q_det_theoretical = np.std(Q_det_theoretical_mat, axis=1)
    std_Q_angle_theoretical = np.std(Q_angle_theoretical_mat, axis=1)

    mean_trace_theoretical_total = np.mean(trace_theoretical_total_mat, axis=1)
    mean_antisym_xy_active = np.mean(antisym_xy_mat_active, axis=1)
    mean_shear_xx_active = np.mean(shear_xx_mat_active, axis=1)
    mean_shear_xy_active = np.mean(shear_xy_mat_active, axis=1)
    mean_shear_det_active = np.mean(shear_det_mat_active, axis=1)
    mean_shear_angle_active = np.mean(shear_angle_mat_active, axis=1)
    mean_bin_pos_active = np.mean(bin_pos_mat_active, axis=1)
    mean_trace_total = np.mean(trace_total_mat, axis=1)
    mean_antisym_xy_passive = np.mean(antisym_xy_mat_passive, axis=1)
    mean_shear_xx_passive = np.mean(shear_xx_mat_passive, axis=1)
    mean_shear_xy_passive = np.mean(shear_xy_mat_passive, axis=1)
    mean_shear_det_passive = np.mean(shear_det_mat_passive, axis=1)
    mean_shear_angle_passive = np.mean(shear_angle_mat_passive, axis=1)
    mean_bin_pos_passive = np.mean(bin_pos_mat_passive, axis=1)
    mean_Q_det_theoretical = np.mean(Q_det_theoretical_mat, axis=1)
    mean_Q_angle_theoretical = np.mean(Q_angle_theoretical_mat, axis=1)
    mean_bin_pos_total = np.mean(bin_pos_total_mat, axis=1)

    # Initialize return values
    xData = None
    yData = None

    # Plotting
    if plot_style[0]:  # total trace
        plt.subplot(i, j, k)
        xData = mean_bin_pos_total
        yData = mean_trace_total
        plt.errorbar(mean_bin_pos_total, mean_trace_total, yerr=std_trace_total,
                     fmt="-s", linewidth=1.5, markersize=5, label='trace tot')
        plt.xlim([0, 30])
        plt.ylim([-0.04, 0])
        plt.autoscale(enable=True, axis='y')
        plt.title(f'P0={P0:.1f}, ca1={C1:.2f}')
        plt.xlabel('cell position', fontsize=10)
        plt.ylabel(r'tr($\sigma$)', fontsize=10)
        plt.gca().tick_params(labelsize=10)
        plt.grid(True)

    if plot_style[6]:  # Q angle theoretical
        plt.subplot(i, j, k)
        xData = mean_bin_pos_total
        yData = mean_Q_angle_theoretical
        plt.errorbar(mean_bin_pos_total, mean_Q_angle_theoretical, yerr=std_Q_angle_theoretical,
                     fmt="-s", linewidth=1.5, markersize=5, label='Q angle theo')
        plt.xlim([0, 30])
        plt.ylim([-0.04, 0])
        plt.autoscale(enable=True, axis='y')
        plt.title(f'P0={P0:.1f}, ca1={C1:.2f}')
        plt.xlabel('cell position', fontsize=10)
        plt.ylabel(r'$Q_{\phi}$ theoretical (radians)', fontsize=10)
        plt.gca().tick_params(labelsize=10)
        plt.grid(True)

    if plot_style[7]:  # Q det theoretical
        plt.subplot(i, j, k)
        xData = mean_bin_pos_total
        yData = np.sqrt(-mean_Q_det_theoretical)
        plt.errorbar(mean_bin_pos_total, np.sqrt(-mean_Q_det_theoretical), yerr=0.5 * std_Q_det_theoretical,
                     fmt="-s", linewidth=1.5, markersize=5, label='det(Q) theo')
        plt.xlim([0, 30])
        plt.ylim([-0.04, 0])
        plt.autoscale(enable=True, axis='y')
        plt.title(f'P0={P0:.1f}, ca1={C1:.2f}')
        plt.xlabel('cell position', fontsize=10)
        plt.ylabel(r'$|Q|$ theoretical (radians)', fontsize=10)
        plt.gca().tick_params(labelsize=10)
        plt.grid(True)

    if cell_type[0]:  # active cells
        plt.subplot(i, j, k)
        xData = mean_bin_pos_active
        yData = mean_antisym_xy_active

        if plot_style[1]:  # antisym xy
            plt.errorbar(mean_bin_pos_active, mean_antisym_xy_active, yerr=std_antisym_xy_active,
                         fmt="-s", linewidth=1.5, markersize=5, label='antisym_xy active')

        if plot_style[2]:  # shear xx xy
            plt.errorbar(mean_bin_pos_active, mean_shear_xx_active, yerr=std_shear_xx_active,
                         fmt="-s", linewidth=1.5, markersize=5, label='shear_xx active')
            plt.errorbar(mean_bin_pos_active, mean_shear_xy_active, yerr=std_shear_xy_active,
                         fmt="-s", linewidth=1.5, markersize=5, label='shear_xy active')

        if plot_style[3]:  # shear det
            plt.errorbar(mean_bin_pos_active, mean_shear_det_active, yerr=std_shear_det_active,
                         fmt="-s", linewidth=1.5, markersize=5, label='det(shear) active')

        if plot_style[4]:  # shear angle
            plt.errorbar(mean_bin_pos_active, mean_shear_angle_active, yerr=std_shear_angle_active,
                         fmt="-s", linewidth=1.5, markersize=5, label='angle shear active')

        plt.xlim([0, 30])
        plt.ylim([-0.04, 0])
        plt.autoscale(enable=True, axis='y')
        plt.title(f'P0={P0:.1f}, ca1={C1:.2f}')
        plt.xlabel('cell position', fontsize=10)
        plt.ylabel(r'$\sigma$', fontsize=10)
        if plot_style[4]:
            plt.ylabel('shear part angle (radian)', fontsize=10)
        plt.gca().tick_params(labelsize=10)
        plt.grid(True)

    if cell_type[1]:  # passive cells
        plt.subplot(i, j, k)
        xData = mean_bin_pos_passive
        yData = mean_antisym_xy_passive

        if plot_style[1]:  # antisym xy
            plt.errorbar(mean_bin_pos_passive, mean_antisym_xy_passive, yerr=std_antisym_xy_passive,
                         fmt="-s", linewidth=1.5, markersize=5, label='antisym_xy passive')

        if plot_style[2]:  # shear xx xy
            plt.errorbar(mean_bin_pos_passive, mean_shear_xx_passive, yerr=std_shear_xx_passive,
                         fmt="-s", linewidth=1.5, markersize=5, label='shear_xx passive')
            plt.errorbar(mean_bin_pos_passive, mean_shear_xy_passive, yerr=std_shear_xy_passive,
                         fmt="-s", linewidth=1.5, markersize=5, label='shear_xy passive')

        if plot_style[3]:  # shear det
            plt.errorbar(mean_bin_pos_passive, mean_shear_det_passive, yerr=std_shear_det_passive,
                         fmt="-s", linewidth=1.5, markersize=5, label='shear_det passive')

        if plot_style[4]:  # shear angle
            plt.errorbar(mean_bin_pos_passive, mean_shear_angle_passive, yerr=std_shear_angle_passive,
                         fmt="-s", linewidth=1.5, markersize=5, label='shear angle passive')

        plt.xlim([0, 30])
        plt.ylim([-0.0005, 0.0005])
        plt.autoscale(enable=True, axis='y')
        plt.title(f'P0={P0:.1f}, ca1={C1:.2f}')
        plt.xlabel('cell position', fontsize=10)
        plt.ylabel(r'$\sigma$', fontsize=10)
        if plot_style[4]:
            plt.ylabel('shear part angle (radian)', fontsize=10)
        plt.gca().tick_params(labelsize=10)
        plt.grid(True)

    return yData, xData
