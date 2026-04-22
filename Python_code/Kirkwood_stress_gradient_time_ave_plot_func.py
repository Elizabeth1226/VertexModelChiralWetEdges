import numpy as np
import matplotlib.pyplot as plt
from Kirkwood_stress_func import Kirkwood_stress_func


def numerical_differentiation(D1, data, positions):
    """Helper function for numerical differentiation with protection against division by zero"""
    gradient = np.zeros(D1 - 1)
    for i in range(D1 - 1):
        denom = positions[i + 1] - positions[i]
        if abs(denom) > 1e-10:  # Avoid division by zero
            gradient[i] = (data[i + 1] - data[i]) / denom
        else:
            gradient[i] = 0.0  # Set gradient to zero when positions are too close
    return gradient


def Kirkwood_stress_gradient_time_ave_plot_func(kP, P0, C1, dw, start_time, i, j, k, rho_0, plot_style, filepath):
    """
    Plot time-averaged gradient of Kirkwood stress

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
    rho_0 : float
        Density parameter
    plot_style : list of bool
        [if_plot_r_component, if_plot_phi_component]
    filepath : str
        Full path to the output folder

    Returns:
    --------
    xData, yData : numpy arrays
        X and Y data for the plot

    NOTE: Only sensible to check force balance of Kirkwood stress calculated from
    total cell force, hence force_type=1
    """
    force_type = 1

    # Call Kirkwood_stress_func at multiple time steps
    results_5 = Kirkwood_stress_func(kP, P0, start_time, force_type, rho_0, filepath)
    results_10 = Kirkwood_stress_func(kP, P0, start_time + 1*dw, force_type, rho_0, filepath)
    results_15 = Kirkwood_stress_func(kP, P0, start_time + 2*dw, force_type, rho_0, filepath)
    results_20 = Kirkwood_stress_func(kP, P0, start_time + 3*dw, force_type, rho_0, filepath)
    results_25 = Kirkwood_stress_func(kP, P0, start_time + 4*dw, force_type, rho_0, filepath)

    # Unpack results (indices 0, 6, 9, 12 correspond to trace_total, antisym_xy_total, shear_xx_total, shear_xy_total)
    trace_total_5 = results_5[0]
    antisym_xy_total_5 = results_5[6]
    shear_xx_total_5 = results_5[9]
    shear_xy_total_5 = results_5[12]
    bin_pos_total_5 = results_5[23]

    trace_total_10 = results_10[0]
    antisym_xy_total_10 = results_10[6]
    shear_xx_total_10 = results_10[9]
    shear_xy_total_10 = results_10[12]
    bin_pos_total_10 = results_10[23]

    trace_total_15 = results_15[0]
    antisym_xy_total_15 = results_15[6]
    shear_xx_total_15 = results_15[9]
    shear_xy_total_15 = results_15[12]
    bin_pos_total_15 = results_15[23]

    trace_total_20 = results_20[0]
    antisym_xy_total_20 = results_20[6]
    shear_xx_total_20 = results_20[9]
    shear_xy_total_20 = results_20[12]
    bin_pos_total_20 = results_20[23]

    trace_total_25 = results_25[0]
    antisym_xy_total_25 = results_25[6]
    shear_xx_total_25 = results_25[9]
    shear_xy_total_25 = results_25[12]
    bin_pos_total_25 = results_25[23]

    # Find minimum size
    S = np.array([len(trace_total_5), len(trace_total_10), len(trace_total_15),
                  len(trace_total_20), len(trace_total_25)])
    D1 = np.min(S)

    # Create matrices
    trace_mat = np.column_stack([trace_total_5[:D1], trace_total_10[:D1],
                                  trace_total_15[:D1], trace_total_20[:D1],
                                  trace_total_25[:D1]])
    antisym_xy_mat = np.column_stack([antisym_xy_total_5[:D1], antisym_xy_total_10[:D1],
                                       antisym_xy_total_15[:D1], antisym_xy_total_20[:D1],
                                       antisym_xy_total_25[:D1]])
    shear_xx_mat = np.column_stack([shear_xx_total_5[:D1], shear_xx_total_10[:D1],
                                     shear_xx_total_15[:D1], shear_xx_total_20[:D1],
                                     shear_xx_total_25[:D1]])
    shear_xy_mat = np.column_stack([shear_xy_total_5[:D1], shear_xy_total_10[:D1],
                                     shear_xy_total_15[:D1], shear_xy_total_20[:D1],
                                     shear_xy_total_25[:D1]])
    bin_pos_mat = np.column_stack([bin_pos_total_5[:D1], bin_pos_total_10[:D1],
                                    bin_pos_total_15[:D1], bin_pos_total_20[:D1],
                                    bin_pos_total_25[:D1]])

    # Calculate means and standard deviations
    std_trace = np.std(trace_mat, axis=1)
    std_antisym_xy = np.std(antisym_xy_mat, axis=1)
    std_shear_xx = np.std(shear_xx_mat, axis=1)
    std_shear_xy = np.std(shear_xy_mat, axis=1)
    mean_trace = np.mean(trace_mat, axis=1)
    mean_antisym_xy = np.mean(antisym_xy_mat, axis=1)
    mean_shear_xx = np.mean(shear_xx_mat, axis=1)
    mean_shear_xy = np.mean(shear_xy_mat, axis=1)
    mean_bin_pos = np.mean(bin_pos_mat, axis=1)

    # Calculate stress tensor components
    T_rr = 0.5 * mean_trace + mean_shear_xx
    T_rphi = mean_antisym_xy + mean_shear_xy
    T_phir = -mean_antisym_xy + mean_shear_xy
    T_phiphi = 0.5 * mean_trace - mean_shear_xx

    # Calculate divergence (with protection against division by zero)
    # Forward difference: dV/dr[i] = (V[i+1] - V[i]) / (R[i+1] - R[i]),
    # evaluated at position R[i].  Output length is D1-1.
    dr = np.diff(mean_bin_pos)                        # R[i+1] - R[i], shape (D1-1,)
    safe_dr = np.where(np.abs(dr) > 1e-10, dr, 1.0)  # avoid division by zero
    nonzero_dr = np.abs(dr) > 1e-10
    dT_rr_dr   = np.where(nonzero_dr, np.diff(T_rr)   / safe_dr, 0.0)
    dT_rphi_dr = np.where(nonzero_dr, np.diff(T_rphi) / safe_dr, 0.0)
    r_grad = mean_bin_pos[:-1]   # derivative i is evaluated at R[i], i = 0 ... D1-2

    # Stress tensor components at the same positions where derivatives are evaluated.
    # Forward difference places dV/dr[i] at R[i], so we use T[i] for the geometric term.
    T_rr_at_r = T_rr[:D1-1]
    T_phiphi_at_r = T_phiphi[:D1-1]
    T_rphi_at_r = T_rphi[:D1-1]

    # Check for zero positions in radial component calculation
    zero_pos_mask_r = r_grad == 0
    if np.any(zero_pos_mask_r):
        zero_indices = np.where(zero_pos_mask_r)[0]
        print(f"Warning: Division by zero avoided in div_T_r calculation at {len(zero_indices)} position(s): indices {zero_indices}. Setting terms to 0.")

    # Calculate divergence in cylindrical coordinates:
    # div_T_r = dT_rr/dr + (T_rr - T_phiphi)/r
    div_T_r = dT_rr_dr + np.divide(T_rr_at_r - T_phiphi_at_r, r_grad,
                                    where=r_grad!=0, out=np.zeros(D1-1))

    # Check for zero positions in azimuthal component calculation
    zero_pos_mask_phi = r_grad == 0
    if np.any(zero_pos_mask_phi):
        zero_indices = np.where(zero_pos_mask_phi)[0]
        print(f"Warning: Division by zero avoided in div_T_phi calculation at {len(zero_indices)} position(s): indices {zero_indices}. Setting terms to 0.")

    # Calculate divergence in cylindrical coordinates:
    # div_T_phi = dT_rphi/dr + (T_rphi + T_phir)/r
    # Use tensor components at the same positions where derivatives are evaluated
    T_phir_at_r = T_phir[:D1-1]
    div_T_phi = dT_rphi_dr + np.divide(T_rphi_at_r + T_phir_at_r, r_grad,
                                        where=r_grad!=0, out=np.zeros(D1-1))


    std_div_T_r = 0.5 * std_trace + std_shear_xx
    std_div_T_phi = std_antisym_xy + std_shear_xy

    # Plotting
    plt.subplot(i, j, k)



    # Shift r_grad to bin centres: forward difference straddles [R[i], R[i+1]],
    # so the representative radial position is the midpoint R[i] + dr[i]/2.
    r_grad = r_grad + 0.5 * dr

    # Use positions where gradients are evaluated
    xData = r_grad

    if plot_style[0]:  # radial component
        yData = div_T_r
        plt.errorbar(r_grad, div_T_r, yerr=std_div_T_r[:D1-1],
                     fmt="-s", linewidth=1, markersize=5, label=r'$(\nabla \cdot \mathbf{T})_r$')
        plt.xlim([0, 30])
        plt.ylim([0, 0.01])
        plt.autoscale(enable=True, axis='y')
        plt.title(f'P0={P0:.1f}, ca1={C1:.2f}')
        plt.xlabel('cell position', fontsize=10)
        plt.ylabel(r'$(\nabla \sigma)_r$', fontsize=10)
        plt.gca().tick_params(labelsize=10)
        plt.grid(True)

    if plot_style[1]:  # azimuthal component
        yData = div_T_phi
        plt.errorbar(r_grad, div_T_phi, yerr=std_div_T_phi[:D1-1],
                     fmt="-s", linewidth=1, markersize=5, label=r'$(\nabla \cdot \mathbf{T})_\phi$')
        plt.xlim([0, 30])
        plt.ylim([0, 0.01])
        plt.autoscale(enable=True, axis='y')
        plt.title(f'P0={P0:.1f}, ca1={C1:.2f}')
        plt.xlabel('cell position', fontsize=10)
        plt.ylabel(r'$(\nabla \sigma)_\phi$', fontsize=10)
        plt.gca().tick_params(labelsize=10)
        plt.grid(True)

    return xData, yData
