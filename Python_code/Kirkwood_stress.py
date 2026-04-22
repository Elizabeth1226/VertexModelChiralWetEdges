"""
Kirkwood stress calculations main script

This script performs force balance analysis by computing:
- Kirkwood stress components (trace, antisymmetric, shear)
- Velocity fields (radial and azimuthal)
- Laplacian of velocity fields
- Divergence of stress tensor
- Comparison with theoretical predictions
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')  # Use non-interactive backend for saving figures
import matplotlib.pyplot as plt
from scipy import stats
import pickle
import os

from Kirkwood_stress_ave_time_plot_func import Kirkwood_stress_ave_time_plot_func
from Kirkwood_stress_gradient_time_ave_plot_func import Kirkwood_stress_gradient_time_ave_plot_func
from Vr_time_ave_plot_func import Vr_time_ave_plot_func
from Vphi_time_ave_plot_func import Vphi_time_ave_plot_func
from LapV_time_ave_plot_func import LapV_time_ave_plot_func


# Global parameters
Nx = 64
gammaW = 0.1
zettaW = 1
dynamics = 1
disorder = 0.02
kP = 0.01
P0_list = np.array([3.6, 3.8, 4.1])
C1_list = np.array([0.05, 0.05, 0.05])
disorder_list = np.array([0.02, 0.02])
Color_list_type1 = np.array([[0, 1, 1], [0, 0.4471, 0.7412], [0.3020, 0.7451, 0.9333]])
Color_list_type2 = np.array([[0.6353, 0.0784, 0.1843], [1.0000, 0, 0], [1.0000, 0.4784, 0.5490]])
color1 = np.array([0, 1, 1])
color2 = np.array([1.0000, 0.4784, 0.5490])
P0 = 3.6
C1 = 0.05
dw = 400
di = 400
tMAX = 20000
seed = 1
count = 1
start_time = 18400

# force_type: 1, cell force; 2, cell elastic force; 3, cell perimeter force; 4, cell area force; 5, cell chiral force
force_type = 1
rho_0 = 1

# plot_style = [total trace, antisym_xy, shear_xx_xy, shear_det, shear_angle, fit_trace, Q_angle_theoretical, Q_det_theoretical]
plot_style_tot_trace = [True, False, False, False, False, False, False, False]
plot_style_sym_antisym_components = [False, True, True, False, False, False, False, False]

r = 4
c = 4

# ************************************************************************
# Plotting Vphi and Vr (commented out in original MATLAB)
# P0_list = np.arange(3.6, 4.3, 0.1)
# C1_list = np.arange(0.01, 0.06, 0.01)
# r = 7
# c = 5

# **********************************************************************************************
# Plotting azimuthal av. components of kirkwood stress from various force contributions (commented out in original)
# active_cell_plot = [True, False]
# passive_cell_plot = [False, True]
# Kirkwood_stress_ave_time_plot_func(Nx, kP, P0, C1, dw, di, disorder, tMAX, seed, start_time, 1, 1, 1, passive_cell_plot, force_type, rho_0, plot_style_sym_antisym_components, color1, color2, dynamics, gammaW, zettaW)
# Kirkwood_stress_ave_time_plot_func(Nx, kP, P0, C1, dw, di, disorder, tMAX, seed, start_time, 1, 1, 1, active_cell_plot, force_type, rho_0, plot_style_tot_trace, color1, color2, dynamics, gammaW, zettaW)
# Kirkwood_stress_ave_time_plot_func(Nx, kP, P0, C1, dw, di, disorder, tMAX, seed, start_time, 1, 1, 1, active_cell_plot, force_type, rho_0, plot_style_sym_antisym_components, color1, color2, dynamics, gammaW, zettaW)
# plt.legend()

# **********************************************************************************************
# Force balance check: div(Kirkwood stress) = gammaW * velocity - zeta * lap(velocity)
r = 3
c = len(P0_list)

plot_LapVr = [True, False]
plot_LapVphi = [False, True]

r_pos_list = []
div_T_r_list = []
div_T_phi_list = []
Vr_list = []
Vphi_list = []
LapVr_list = []
LapVphi_list = []

b0_list = []
gamma_list = []
zeta_list = []

b0r_list = []
gammar_list = []
zetar_list = []

r_pos_list_cell = {}
div_T_r_list_cell = {}
div_T_phi_list_cell = {}
Vr_list_cell = {}
Vphi_list_cell = {}
LapVr_list_cell = {}
LapVphi_list_cell = {}

# Create figure for force balance comparison
plt.figure(figsize=(15, 10))

for count in range(len(P0_list)):
    P0 = P0_list[count]
    C1 = C1_list[count]
    disorder = 0.02

    # Get stress gradients
    r_pos, div_T_r = Kirkwood_stress_gradient_time_ave_plot_func(Nx, kP, P0, C1, dw, di, disorder, tMAX, seed, start_time, r, c, count + 1, rho_0, [True, False], dynamics, gammaW, zettaW)
    _, div_T_phi = Kirkwood_stress_gradient_time_ave_plot_func(Nx, kP, P0, C1, dw, di, disorder, tMAX, seed, start_time, r, c, c + count + 1, rho_0, [False, True], dynamics, gammaW, zettaW)

    # Get velocities
    Vr_bin_pos, Vr = Vr_time_ave_plot_func(Nx, kP, P0, C1, dw, di, disorder, tMAX, seed, r, c, count + 1, dynamics, gammaW, zettaW, start_time)
    plt.legend(['radial stress gradient', 'radial velocity'])

    Vphi_bin_pos, Vphi = Vphi_time_ave_plot_func(Nx, kP, P0, C1, dw, di, disorder, tMAX, seed, r, c, c + count + 1, dynamics, gammaW, zettaW, start_time)
    plt.legend(['azimuthal stress gradient', 'azimuthal velocity'])

    # Get Laplacians of velocity
    LapV_bin_pos, LapVr, LapVphi = LapV_time_ave_plot_func(Nx, kP, P0, C1, dw, di, disorder, tMAX, seed, r, c, c + count + 1, plot_LapVphi, dynamics, gammaW, zettaW, start_time)
    plt.legend()

    _, _, _ = LapV_time_ave_plot_func(Nx, kP, P0, C1, dw, di, disorder, tMAX, seed, r, c, count + 1, plot_LapVr, dynamics, gammaW, zettaW, start_time)
    plt.legend()

    # Store data for later analysis
    r_pos_list_cell[count] = r_pos
    div_T_r_list_cell[count] = div_T_r
    div_T_phi_list_cell[count] = div_T_phi
    Vr_list_cell[count] = Vr
    Vphi_list_cell[count] = Vphi
    LapVr_list_cell[count] = LapVr
    LapVphi_list_cell[count] = LapVphi

    # Fit azimuthal components
    x = Vphi[:len(div_T_phi)]
    y = LapVphi[:len(div_T_phi)]
    z = div_T_phi

    # Create design matrix for linear regression
    X = np.column_stack([np.ones(len(x)), x, y])

    # Perform linear regression
    result = np.linalg.lstsq(X, z, rcond=None)
    coeffs = result[0]

    b0 = coeffs[0]
    b1 = coeffs[1]
    b2 = coeffs[2]

    print(f'-------------------------------')
    print(f'azimuthal P0 {P0}')
    print(f'Intercept (b0) = {b0:.4f}')
    print(f'Coefficient for x (b1) = {b1:.4f}')
    print(f'Coefficient for y (b2) = {b2:.4f}')

    # Fit radial components
    xr = Vr[:len(div_T_r)]
    yr = LapVr[:len(div_T_r)]
    zr = div_T_r

    Xr = np.column_stack([np.ones(len(xr)), xr, yr])
    result_r = np.linalg.lstsq(Xr, zr, rcond=None)
    coeffs_r = result_r[0]

    b0r = coeffs_r[0]
    b1r = coeffs_r[1]
    b2r = coeffs_r[2]

    print(f'-------------------------------')
    print(f'radial P0 {P0}')
    print(f'Intercept (b0) = {b0r:.4f}')
    print(f'Coefficient for x (b1) = {b1r:.4f}')
    print(f'Coefficient for y (b2) = {b2r:.4f}')

    # Plot comparison
    plt.subplot(r, c, 2*c + count + 1)
    plt.plot(r_pos, div_T_r[1:], label='div(Stress)_r')
    plt.plot(r_pos, Vr[1:len(div_T_phi)] * b1r + LapVr[1:len(div_T_phi)] * b2r, label='gamma*v_r-zeta*Lap(v)_r')
    plt.plot(r_pos, z[1:], label='div(Stress)_phi')
    plt.plot(r_pos, x[1:] * b1 + y[1:] * b2, label='gamma*v_phi-zeta*Lap(v)_phi')
    plt.title(f'gamma_r,phi={b1r:.2f}, {b1:.2f} zeta_r,phi={-b2r:.2f}, {-b2:.2f}')
    plt.xlabel('cell position')
    plt.legend()
    plt.xlim([0, 20])

    b0_list.append(b0)
    gamma_list.append(b1)
    zeta_list.append(-b2)

    b0r_list.append(b0r)
    gammar_list.append(b1r)
    zetar_list.append(-b2r)

try:
    plt.savefig('force_balance_comp.png', dpi=300, bbox_inches='tight')
    print(f"Saved force_balance_comp.png to {os.getcwd()}")
except Exception as e:
    print(f"Error saving force_balance_comp.png: {e}")

# Plot gamma and zeta vs P0
plt.figure()
plt.plot(P0_list, gamma_list, label='gamma_phi')
plt.plot(P0_list, zeta_list, label='zeta_phi')
plt.title(f'azimuthal, C1{C1}')
plt.legend()
try:
    plt.savefig('force_balance_r.png', dpi=300, bbox_inches='tight')
    print(f"Saved force_balance_r.png")
except Exception as e:
    print(f"Error saving force_balance_r.png: {e}")

plt.figure()
plt.plot(P0_list, gammar_list, label='gamma_r')
plt.plot(P0_list, zetar_list, label='zeta_r')
plt.title(f'radial, C1{C1}')
plt.legend()
try:
    plt.savefig('force_balance_phi.png', dpi=300, bbox_inches='tight')
    print(f"Saved force_balance_phi.png")
except Exception as e:
    print(f"Error saving force_balance_phi.png: {e}")

# Save data
data_dict = {
    'P0_list': P0_list,
    'C1_list': C1_list,
    'r_pos_list_cell': r_pos_list_cell,
    'div_T_r_list_cell': div_T_r_list_cell,
    'div_T_phi_list_cell': div_T_phi_list_cell,
    'Vr_list_cell': Vr_list_cell,
    'Vphi_list_cell': Vphi_list_cell,
    'LapVr_list_cell': LapVr_list_cell,
    'LapVphi_list_cell': LapVphi_list_cell,
    'gamma_list': gamma_list,
    'zeta_list': zeta_list,
    'gammar_list': gammar_list,
    'zetar_list': zetar_list
}

try:
    with open('data.pkl', 'wb') as f:
        pickle.dump(data_dict, f)
    print(f"Data saved to data.pkl in {os.getcwd()}")
except Exception as e:
    print(f"Error saving data.pkl: {e}")

# ******************************************************************
# Debugging Kirkwood stress function (commented out in original MATLAB)
# This section contains extensive debugging code that was commented out in the original

print("\n=== Analysis Complete ===")
print("Force balance analysis has been completed.")
print("Results saved to:")
print("  - force_balance_comp.png")
print("  - force_balance_r.png")
print("  - force_balance_phi.png")
print("  - data.pkl")
