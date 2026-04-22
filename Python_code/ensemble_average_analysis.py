"""
Ensemble Average Analysis  (equivalent to Section 4 of Kirkwood_stress.ipynb)

Computes and plots:
  1. Ensemble-averaged vertex forces, vertex density, and Lap(v)
  2. Ensemble-averaged Kirkwood stress divergences (total, chiral, perimeter, area forces)
  3. Dissipative Kirkwood stress divergence and combined (total + dissipative) divergence
  4. Divergence of the ensemble-mean stress tensor for each force type
  5. Hypothesis test:  div(stress) - rho_v * vF = -k * zetaW * Lap(vF)
  6. Hypothesis test:  div(T_total + T_diss) = k * vF

Usage:
    python ensemble_average_analysis.py
    or import and call run_ensemble_average_analysis(params) directly.
"""

import sys
import os
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

# ── path setup (adjust if running from a different directory) ──────────────
_SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
if _SCRIPT_DIR not in sys.path:
    sys.path.insert(0, _SCRIPT_DIR)

from Kirkwood_stress_func             import Kirkwood_stress_func
from Kirkwood_stress_dissipation_func import Kirkwood_stress_dissipation_func
from vF_ave_func                      import vF_ave_func
from LapV_ave_func                    import LapV_ave_func
from vertex_density_ave_func          import vertex_density_ave_func


# ═══════════════════════════════════════════════════════════════════════════
# Default simulation parameters  (edit here or pass as kwargs to the runner)
# ═══════════════════════════════════════════════════════════════════════════
DEFAULT_PARAMS = dict(
    iT_list         = list(range(1, 11)),
    Nx              = 60,
    gammaW          = 1,
    zettaW          = 1e-1,
    disorder        = 0.02,
    kP              = 0.01,
    dynamics        = 1,           # 1=WVM vertex friction
    transition_type = 0,           # 0=sharp
    P0              = 3.8,
    C1              = 0.03,
    dw              = 50,
    di              = 50,
    tMAX            = 1000,
    seed            = 1,
    Lb              = 1,
    rho_0           = 1,
    output_root     = '/mnt/users/wangg/VertexModelChiralWetEdges/output',
)


# ═══════════════════════════════════════════════════════════════════════════
# Internal helpers
# ═══════════════════════════════════════════════════════════════════════════

def _folder_path(p, iT):
    return (
        f"{p['output_root']}/"
        f"out_Nx_{p['Nx']}_kA_0.5_kP_{p['kP']}_P0_{p['P0']}_"
        f"disorder_{p['disorder']}_tMAX_{p['tMAX']}_dw_{p['dw']}_"
        f"di_{p['di']}_ca1_{p['C1']}_ca2_0_iC_0_iP_0.25_"
        f"iS_{p['seed']}_iT_{iT}_dynamics_{p['dynamics']}_"
        f"gammaW_{p['gammaW']}_zetaW_{p['zettaW']}_tT_{p['transition_type']}/"
    )


def _trim_and_stack(lst):
    """Trim each array to the shortest length, then stack into a 2-D array."""
    min_len = min(arr.shape[0] for arr in lst)
    return np.array([arr[:min_len] for arr in lst])


def _stress_tensor(trace_total, antisym_xy_total, shear_xx_total, shear_xy_total):
    """Unpack Kirkwood stress-function output into (T_rr, T_rphi, T_phir, T_phiphi)."""
    T_rr     =  0.5 * trace_total + shear_xx_total
    T_rphi   =  antisym_xy_total  + shear_xy_total
    T_phir   = -antisym_xy_total  + shear_xy_total
    T_phiphi =  0.5 * trace_total - shear_xx_total
    return T_rr, T_rphi, T_phir, T_phiphi


def _forward_div(T_rr, T_rphi, T_phir, T_phiphi, bin_pos):
    """Forward-difference divergence of a cylindrical stress tensor."""
    n      = len(bin_pos)
    dr     = np.diff(bin_pos)
    safe   = np.where(np.abs(dr) > 1e-10, dr, 1.0)
    nz     = np.abs(dr) > 1e-10
    dTrr   = np.where(nz, np.diff(T_rr)   / safe, 0.0)
    dTrphi = np.where(nz, np.diff(T_rphi) / safe, 0.0)
    r      = bin_pos[:-1]
    div_r   = dTrr   + np.divide(T_rr[:n-1] - T_phiphi[:n-1], r,
                                  where=r != 0, out=np.zeros(n-1))
    div_phi = dTrphi + np.divide(T_rphi[:n-1] + T_phir[:n-1],  r,
                                  where=r != 0, out=np.zeros(n-1))
    return div_r, div_phi, r          # r = divergence radial grid


def _numerical_diff_err(T_mean, pos):
    """Forward-difference truncation error |T''(r)| × dr at each div-grid point.

    Returns array of length len(pos)-1, aligned with the divergence grid.
    """
    n   = len(T_mean)
    dr  = np.diff(pos)
    d2T = np.zeros(n)

    drl, drr = dr[:-1], dr[1:]
    d2T[1:-1] = 2.0 * (
        T_mean[2:]   / (drr * (drl + drr)) -
        T_mean[1:-1] / (drl * drr) +
        T_mean[:-2]  / (drl * (drl + drr))
    )
    h1, h2 = dr[0], dr[1]
    d2T[0] = 2.0 * (
        T_mean[0] / (h1*(h1+h2)) - T_mean[1] / (h1*h2) + T_mean[2] / (h2*(h1+h2))
    )
    h1, h2 = dr[-2], dr[-1]
    d2T[-1] = 2.0 * (
        T_mean[-3] / (h1*(h1+h2)) - T_mean[-2] / (h1*h2) + T_mean[-1] / (h2*(h1+h2))
    )
    return np.abs(d2T[:-1]) * dr


def _num_diff_err_full(T_mean, pos):
    """Same as _numerical_diff_err but returns length len(pos) (full grid)."""
    n   = len(T_mean)
    dr  = np.diff(pos)

    dr_eff = np.empty(n)
    dr_eff[:-1] = dr
    dr_eff[-1]  = dr[-1]

    d2T = np.zeros(n)
    drl, drr = dr[:-1], dr[1:]
    d2T[1:-1] = 2.0 * (
        T_mean[2:]   / (drr * (drl + drr)) -
        T_mean[1:-1] / (drl * drr) +
        T_mean[:-2]  / (drl * (drl + drr))
    )
    h1, h2 = dr[0], dr[1]
    d2T[0]  = 2.0 * (T_mean[0] / (h1*(h1+h2)) - T_mean[1]  / (h1*h2) + T_mean[2]  / (h2*(h1+h2)))
    h1, h2  = dr[-2], dr[-1]
    d2T[-1] = 2.0 * (T_mean[-3]/ (h1*(h1+h2)) - T_mean[-2] / (h1*h2) + T_mean[-1] / (h2*(h1+h2)))
    return np.abs(d2T) * dr_eff


def _first_deriv(T, pos):
    """First derivative: forward difference (i=0..n-2), backward (i=n-1)."""
    dr = np.diff(pos)
    dT = np.empty(len(T))
    dT[:-1] = np.diff(T) / dr
    dT[-1]  = (T[-1] - T[-2]) / dr[-1]
    return dT


def _first_deriv_std(T_std, pos):
    """Propagated std through _first_deriv."""
    dr    = np.diff(pos)
    sigma = np.empty(len(T_std))
    sigma[:-1] = np.sqrt(T_std[1:]**2 + T_std[:-1]**2) / dr
    sigma[-1]  = np.sqrt(T_std[-1]**2 + T_std[-2]**2)  / dr[-1]
    return sigma


def _div_from_mean(Trr_m, Tpp_m, Trp_m, Tpr_m,
                   Trr_s, Tpp_s, Trp_s, Tpr_s, pos):
    """Divergence of the ensemble-mean stress tensor with propagated uncertainty.

    Returns: div_r, div_phi, sigma_r, sigma_phi  — all at positions pos.
    """
    safe_pos = np.where(pos > 0, pos, np.inf)

    dTrr_dr = _first_deriv(Trr_m, pos)
    dTrp_dr = _first_deriv(Trp_m, pos)
    div_r   = dTrr_dr + (Trr_m - Tpp_m) / safe_pos
    div_phi = dTrp_dr + (Trp_m + Tpr_m) / safe_pos

    num_r   = _num_diff_err_full(Trr_m, pos)
    num_phi = _num_diff_err_full(Trp_m, pos)

    sig_dTrr    = _first_deriv_std(Trr_s, pos)
    sig_dTrp    = _first_deriv_std(Trp_s, pos)
    sig_alg_r   = np.sqrt(Trr_s**2 + Tpp_s**2) / safe_pos
    sig_alg_phi = np.sqrt(Trp_s**2 + Tpr_s**2) / safe_pos

    sig_r   = np.sqrt(sig_dTrr**2 + sig_alg_r**2   + num_r**2)
    sig_phi = np.sqrt(sig_dTrp**2 + sig_alg_phi**2 + num_phi**2)
    return div_r, div_phi, sig_r, sig_phi


def _shaded(ax, x, mean, std, fmt, label, color):
    """Plot line + shaded ±1σ band."""
    ax.plot(x, mean, fmt, label=label, ms=4)
    ax.fill_between(x, mean - std, mean + std, alpha=0.2, color=color)


# ═══════════════════════════════════════════════════════════════════════════
# Hypothesis fit functions
# ═══════════════════════════════════════════════════════════════════════════

def _hypothesis_fit(div_mean, sigma_div,
                    vd_mean,  sigma_vd,
                    vF_mean,  sigma_vF,
                    LapV_mean, sigma_LapV,
                    zetaW, r_pos_div,
                    cbins_LapV_ref, cbins_vd_ref, cbins_ref, R_fit,
                    label=''):
    """
    Test: div(stress) - rho_v * vF = -k * zetaW * Lap(vF)

    Arrays interpolated to cbins_LapV_ref before fitting.
    Fit restricted to r < R_fit; error-weighted least squares (no intercept).

    Returns: r_c, LHS, sigma_LHS, X, sigma_X, k, sigma_k
    """
    r_c     = cbins_LapV_ref[:LapV_mean.shape[0]]
    n_div   = div_mean.shape[0]
    r_div   = r_pos_div[:n_div]

    div_i   = np.interp(r_c, r_div,                       div_mean)
    sig_di  = np.interp(r_c, r_div,                       sigma_div[:n_div])
    vd_i    = np.interp(r_c, cbins_vd_ref[:vd_mean.shape[0]],   vd_mean)
    sig_vdi = np.interp(r_c, cbins_vd_ref[:sigma_vd.shape[0]],  sigma_vd)
    vF_i    = np.interp(r_c, cbins_ref[:vF_mean.shape[0]],      vF_mean)
    sig_vFi = np.interp(r_c, cbins_ref[:sigma_vF.shape[0]],     sigma_vF)
    lap_i   = LapV_mean[:len(r_c)]
    sig_li  = sigma_LapV[:len(r_c)]

    LHS       = div_i - vd_i * vF_i
    sigma_LHS = np.sqrt(sig_di**2 + (sig_vdi * vF_i)**2 + (vd_i * sig_vFi)**2)
    sigma_LHS = np.where(sigma_LHS > 0, sigma_LHS, 1e-30)

    X       = -zetaW * lap_i
    sigma_X =  zetaW * sig_li

    mask   = r_c < R_fit
    Lm, sLm = LHS[mask], sigma_LHS[mask]
    Xm, sXm = X[mask],   sigma_X[mask]

    W          = 1.0 / sLm**2
    denom      = np.sum(W * Xm**2)
    k          = np.sum(W * Lm * Xm) / denom
    sigma_k_WLS = 1.0 / np.sqrt(denom)
    sigma_k_X   = np.sqrt(np.sum((W * sXm * (Lm - 2*k*Xm))**2)) / denom
    sigma_k     = np.sqrt(sigma_k_WLS**2 + sigma_k_X**2)

    print(f'{label:32s}  k = {k:.4f} +/- {sigma_k:.4f}'
          f'  (fit: {mask.sum()}/{len(r_c)} pts, r < {R_fit:.2f})')
    return r_c, LHS, sigma_LHS, X, sigma_X, k, sigma_k


def _combined_stress_fit(div_comb_mean,
                         sigma_div_total, r_pos_div_total,
                         sigma_div_diss,  r_pos_div_diss,
                         vF_mean, sigma_vF,
                         r_pos_div, cbins_ref, R_fit,
                         label=''):
    """
    Test: div(T_total + T_diss) = k * vF

    LHS uncertainty is the quadrature sum of the total-force and dissipative
    stress divergence uncertainties, propagated explicitly:
        sigma_LHS = sqrt(sigma_div_total^2 + sigma_div_diss^2)

    All arrays interpolated to cbins_ref (vF grid) before fitting.
    Fit restricted to r < R_fit; error-weighted least squares (no intercept).

    Returns: r_c, div_interp, sigma_LHS, vF_interp, sigma_vF_interp, k, sigma_k
    """
    r_c   = cbins_ref[:min(vF_mean.shape[0], div_comb_mean.shape[0])]
    n_div = div_comb_mean.shape[0]
    r_div = r_pos_div[:n_div]

    div_i    = np.interp(r_c, r_div, div_comb_mean)

    # Propagate total-force and dissipative stress divergence stds in quadrature
    sig_tot  = np.interp(r_c, r_pos_div_total[:len(sigma_div_total)], sigma_div_total)
    sig_diss = np.interp(r_c, r_pos_div_diss[:len(sigma_div_diss)],   sigma_div_diss)
    sig_di   = np.sqrt(sig_tot**2 + sig_diss**2)
    sig_di   = np.where(sig_di > 0, sig_di, 1e-30)

    vF_i    = vF_mean[:len(r_c)]
    sig_vFi = sigma_vF[:len(r_c)]

    mask   = r_c < R_fit
    Dm, sDm = div_i[mask],  sig_di[mask]
    Vm, sVm = vF_i[mask],   sig_vFi[mask]

    W          = 1.0 / sDm**2
    denom      = np.sum(W * Vm**2)
    k          = np.sum(W * Dm * Vm) / denom
    sigma_k_WLS = 1.0 / np.sqrt(denom)
    sigma_k_vF  = np.sqrt(np.sum((W * sVm * (Dm - 2*k*Vm))**2)) / denom
    sigma_k     = np.sqrt(sigma_k_WLS**2 + sigma_k_vF**2)

    print(f'{label:32s}  k = {k:.4f} +/- {sigma_k:.4f}'
          f'  (fit: {mask.sum()}/{len(r_c)} pts, r < {R_fit:.2f})')
    return r_c, div_i, sig_di, vF_i, sig_vFi, k, sigma_k


# ═══════════════════════════════════════════════════════════════════════════
# Main analysis function
# ═══════════════════════════════════════════════════════════════════════════

def run_ensemble_average_analysis(params=None):
    """Run the full ensemble-average analysis and produce all plots.

    Parameters
    ----------
    params : dict, optional
        Simulation parameters.  Any key not provided falls back to DEFAULT_PARAMS.

    Returns
    -------
    results : dict
        All ensemble-averaged arrays and fit results, keyed by name.
    """
    p = {**DEFAULT_PARAMS, **(params or {})}

    iT_list = p['iT_list']
    N       = len(iT_list)
    Nx      = p['Nx']
    R       = 0.5 * 0.25 * Nx      # active droplet radius
    r_a     = 2.0 * R               # plotting x-limit
    R_fit   = 4 * R                 # fit cutoff

    # ── per-iT accumulation lists ─────────────────────────────────────────
    v_fr_all         = []
    v_fphi_all       = []
    vertex_density_all = []
    uncertainty_vd_all = []
    LapVr_all         = []
    LapVphi_all       = []
    sigma_LapVr_all   = []
    sigma_LapVphi_all = []

    # stress tensor components per iT
    T_rr_total_all     = []; T_rphi_total_all   = []
    T_phir_total_all   = []; T_phiphi_total_all = []
    T_rr_chiral_all    = []; T_rphi_chiral_all  = []
    T_phir_chiral_all  = []; T_phiphi_chiral_all= []
    T_rr_perim_all     = []; T_rphi_perim_all   = []
    T_phir_perim_all   = []; T_phiphi_perim_all = []
    T_rr_area_all      = []; T_rphi_area_all    = []
    T_phir_area_all    = []; T_phiphi_area_all  = []
    T_rr_diss_all      = []; T_rphi_diss_all    = []
    T_phir_diss_all    = []; T_phiphi_diss_all  = []

    # divergence per iT
    div_r_total_all    = []; div_phi_total_all  = []
    div_r_chiral_all   = []; div_phi_chiral_all = []
    div_r_perim_all    = []; div_phi_perim_all  = []
    div_r_area_all     = []; div_phi_area_all   = []
    div_r_diss_all     = []; div_phi_diss_all   = []
    div_r_comb_all     = []; div_phi_comb_all   = []

    # reference grids (captured from first iT)
    cbins_ref     = None
    cbins_vd_ref  = None
    cbins_LapV_ref = None
    bin_pos_ref   = None
    r_grad_ref    = None

    # ── loop over seeds ───────────────────────────────────────────────────
    for iT in iT_list:
        print(f"Processing iT {iT}...")
        fp = _folder_path(p, iT)

        # vertex forces
        cbins_s, v_fr_s, v_fphi_s = vF_ave_func(fp, p['tMAX'], p['Lb'])
        v_fr_all.append(v_fr_s);  v_fphi_all.append(v_fphi_s)
        if cbins_ref is None:
            cbins_ref = cbins_s

        # vertex density
        vd_s, bc_vd_s, unc_vd_s = vertex_density_ave_func(fp, p['tMAX'], p['Lb'])
        vertex_density_all.append(vd_s); uncertainty_vd_all.append(unc_vd_s)
        if cbins_vd_ref is None:
            cbins_vd_ref = bc_vd_s

        # Laplacian of velocity
        cbins_lap_s, LapVr_s, LapVphi_s, sig_LapVr_s, sig_LapVphi_s = \
            LapV_ave_func(fp, p['tMAX'], p['Lb'])
        LapVr_all.append(LapVr_s); LapVphi_all.append(LapVphi_s)
        sigma_LapVr_all.append(sig_LapVr_s); sigma_LapVphi_all.append(sig_LapVphi_s)
        if cbins_LapV_ref is None:
            cbins_LapV_ref = cbins_lap_s

        # helper: call Kirkwood_stress_func for one force type
        def _kstress(force_type):
            (tr, tr_a, tr_p,
             tr_th, tr_th_a, tr_th_p,
             as_xy, as_xy_a, as_xy_p,
             sh_xx, sh_xx_a, sh_xx_p,
             sh_xy, sh_xy_a, sh_xy_p,
             sh_det_a, sh_det_p, sh_ang_a, sh_ang_p,
             Q_det_th, Q_ang_th, Q_det, Q_ang,
             bp_tot, bp_act, bp_pas) = Kirkwood_stress_func(
                p['kP'], p['P0'], p['tMAX'], force_type, p['rho_0'], fp, p['Lb'])
            return _stress_tensor(tr, as_xy, sh_xx, sh_xy), bp_tot

        # total force (force_type=1)
        (T_rr, T_rphi, T_phir, T_phiphi), bp_tot = _kstress(1)
        T_rr_total_all.append(T_rr); T_rphi_total_all.append(T_rphi)
        T_phir_total_all.append(T_phir); T_phiphi_total_all.append(T_phiphi)
        div_r_tot, div_phi_tot, r_grad = _forward_div(T_rr, T_rphi, T_phir, T_phiphi, bp_tot)
        div_r_total_all.append(div_r_tot); div_phi_total_all.append(div_phi_tot)
        if r_grad_ref is None: r_grad_ref = r_grad
        if bin_pos_ref is None: bin_pos_ref = bp_tot

        # chiral force (force_type=5)
        (T_rr, T_rphi, T_phir, T_phiphi), bp = _kstress(5)
        T_rr_chiral_all.append(T_rr); T_rphi_chiral_all.append(T_rphi)
        T_phir_chiral_all.append(T_phir); T_phiphi_chiral_all.append(T_phiphi)
        dr, dp, _ = _forward_div(T_rr, T_rphi, T_phir, T_phiphi, bp)
        div_r_chiral_all.append(dr); div_phi_chiral_all.append(dp)

        # perimeter force (force_type=3)
        (T_rr, T_rphi, T_phir, T_phiphi), bp = _kstress(3)
        T_rr_perim_all.append(T_rr); T_rphi_perim_all.append(T_rphi)
        T_phir_perim_all.append(T_phir); T_phiphi_perim_all.append(T_phiphi)
        dr, dp, _ = _forward_div(T_rr, T_rphi, T_phir, T_phiphi, bp)
        div_r_perim_all.append(dr); div_phi_perim_all.append(dp)

        # area force (force_type=4)
        (T_rr, T_rphi, T_phir, T_phiphi), bp = _kstress(4)
        T_rr_area_all.append(T_rr); T_rphi_area_all.append(T_rphi)
        T_phir_area_all.append(T_phir); T_phiphi_area_all.append(T_phiphi)
        dr, dp, _ = _forward_div(T_rr, T_rphi, T_phir, T_phiphi, bp)
        div_r_area_all.append(dr); div_phi_area_all.append(dp)

        # dissipative stress
        (tr_d, tr_d_a, tr_d_p,
         tr_th_d, tr_th_d_a, tr_th_d_p,
         as_d, as_d_a, as_d_p,
         sh_xx_d, sh_xx_d_a, sh_xx_d_p,
         sh_xy_d, sh_xy_d_a, sh_xy_d_p,
         sh_det_d_a, sh_det_d_p, sh_ang_d_a, sh_ang_d_p,
         Q_det_th_d, Q_ang_th_d, Q_det_d, Q_ang_d,
         bp_d, bp_d_a, bp_d_p) = Kirkwood_stress_dissipation_func(
            p['kP'], p['P0'], p['tMAX'], p['rho_0'], fp, p['Lb'], p['gammaW'])
        T_rr_d, T_rphi_d, T_phir_d, T_phiphi_d = _stress_tensor(tr_d, as_d, sh_xx_d, sh_xy_d)
        T_rr_diss_all.append(T_rr_d); T_rphi_diss_all.append(T_rphi_d)
        T_phir_diss_all.append(T_phir_d); T_phiphi_diss_all.append(T_phiphi_d)
        dr_d, dp_d, _ = _forward_div(T_rr_d, T_rphi_d, T_phir_d, T_phiphi_d, bp_d)
        div_r_diss_all.append(dr_d); div_phi_diss_all.append(dp_d)

        # combined (total + dissipative) per iT — aligned to shorter grid
        n_c = min(len(div_r_tot), len(dr_d))
        div_r_comb_all.append(div_r_tot[:n_c] + dr_d[:n_c])
        div_phi_comb_all.append(div_phi_tot[:n_c] + dp_d[:n_c])

    print("Done processing all iTs.")

    # ── trim & stack ──────────────────────────────────────────────────────
    v_fr_all        = _trim_and_stack(v_fr_all)
    v_fphi_all      = _trim_and_stack(v_fphi_all)
    vertex_density_all = _trim_and_stack(vertex_density_all)
    uncertainty_vd_all = _trim_and_stack(uncertainty_vd_all)
    LapVr_all       = _trim_and_stack(LapVr_all)
    LapVphi_all     = _trim_and_stack(LapVphi_all)
    sigma_LapVr_all = _trim_and_stack(sigma_LapVr_all)
    sigma_LapVphi_all = _trim_and_stack(sigma_LapVphi_all)

    def _stack_all(*pairs):
        return [_trim_and_stack(lst) for lst in pairs]

    (T_rr_total_all, T_rphi_total_all, T_phir_total_all, T_phiphi_total_all,
     T_rr_chiral_all, T_rphi_chiral_all, T_phir_chiral_all, T_phiphi_chiral_all,
     T_rr_perim_all,  T_rphi_perim_all,  T_phir_perim_all,  T_phiphi_perim_all,
     T_rr_area_all,   T_rphi_area_all,   T_phir_area_all,   T_phiphi_area_all,
     T_rr_diss_all,   T_rphi_diss_all,   T_phir_diss_all,   T_phiphi_diss_all,
     div_r_total_all, div_phi_total_all,
     div_r_chiral_all, div_phi_chiral_all,
     div_r_perim_all,  div_phi_perim_all,
     div_r_area_all,   div_phi_area_all,
     div_r_diss_all,   div_phi_diss_all,
     div_r_comb_all,   div_phi_comb_all) = _stack_all(
        T_rr_total_all, T_rphi_total_all, T_phir_total_all, T_phiphi_total_all,
        T_rr_chiral_all, T_rphi_chiral_all, T_phir_chiral_all, T_phiphi_chiral_all,
        T_rr_perim_all,  T_rphi_perim_all,  T_phir_perim_all,  T_phiphi_perim_all,
        T_rr_area_all,   T_rphi_area_all,   T_phir_area_all,   T_phiphi_area_all,
        T_rr_diss_all,   T_rphi_diss_all,   T_phir_diss_all,   T_phiphi_diss_all,
        div_r_total_all, div_phi_total_all,
        div_r_chiral_all, div_phi_chiral_all,
        div_r_perim_all,  div_phi_perim_all,
        div_r_area_all,   div_phi_area_all,
        div_r_diss_all,   div_phi_diss_all,
        div_r_comb_all,   div_phi_comb_all,
    )

    # align reference grids to trimmed arrays
    n_trim      = T_rr_total_all.shape[1]
    bin_pos_ref = bin_pos_ref[:n_trim]
    r_grad_ref  = r_grad_ref[:n_trim - 1]

    # ── ensemble means & stds ─────────────────────────────────────────────
    def _ms(a): return np.mean(a, axis=0), np.std(a, axis=0)

    v_fr_mean,   v_fr_std   = _ms(v_fr_all)
    v_fphi_mean, v_fphi_std = _ms(v_fphi_all)

    vd_mean, vd_std = _ms(vertex_density_all)
    vd_sem          = vd_std / np.sqrt(N)
    unc_vd_combined = np.sqrt(vd_sem**2 + np.mean(uncertainty_vd_all**2, axis=0))

    LapVr_mean,   LapVr_std   = _ms(LapVr_all)
    LapVphi_mean, LapVphi_std = _ms(LapVphi_all)
    sig_LapVr_combined  = np.sqrt(LapVr_std**2  / N + np.mean(sigma_LapVr_all**2,  axis=0))
    sig_LapVphi_combined= np.sqrt(LapVphi_std**2/ N + np.mean(sigma_LapVphi_all**2,axis=0))

    T_rr_total_mean,     T_rr_total_std     = _ms(T_rr_total_all)
    T_rphi_total_mean,   T_rphi_total_std   = _ms(T_rphi_total_all)
    T_phir_total_mean,   T_phir_total_std   = _ms(T_phir_total_all)
    T_phiphi_total_mean, T_phiphi_total_std = _ms(T_phiphi_total_all)

    T_rr_chiral_mean,     T_rr_chiral_std     = _ms(T_rr_chiral_all)
    T_rphi_chiral_mean,   T_rphi_chiral_std   = _ms(T_rphi_chiral_all)
    T_phir_chiral_mean,   T_phir_chiral_std   = _ms(T_phir_chiral_all)
    T_phiphi_chiral_mean, T_phiphi_chiral_std = _ms(T_phiphi_chiral_all)

    T_rr_perim_mean,     T_rr_perim_std     = _ms(T_rr_perim_all)
    T_rphi_perim_mean,   T_rphi_perim_std   = _ms(T_rphi_perim_all)
    T_phir_perim_mean,   T_phir_perim_std   = _ms(T_phir_perim_all)
    T_phiphi_perim_mean, T_phiphi_perim_std = _ms(T_phiphi_perim_all)

    T_rr_area_mean,     T_rr_area_std     = _ms(T_rr_area_all)
    T_rphi_area_mean,   T_rphi_area_std   = _ms(T_rphi_area_all)
    T_phir_area_mean,   T_phir_area_std   = _ms(T_phir_area_all)
    T_phiphi_area_mean, T_phiphi_area_std = _ms(T_phiphi_area_all)

    T_rr_diss_mean,     T_rr_diss_std     = _ms(T_rr_diss_all)
    T_rphi_diss_mean,   T_rphi_diss_std   = _ms(T_rphi_diss_all)
    T_phir_diss_mean,   T_phir_diss_std   = _ms(T_phir_diss_all)
    T_phiphi_diss_mean, T_phiphi_diss_std = _ms(T_phiphi_diss_all)

    div_r_total_mean,  div_r_total_std  = _ms(div_r_total_all)
    div_phi_total_mean,div_phi_total_std= _ms(div_phi_total_all)
    div_r_chiral_mean, div_r_chiral_std = _ms(div_r_chiral_all)
    div_phi_chiral_mean,div_phi_chiral_std = _ms(div_phi_chiral_all)
    div_r_perim_mean,  div_r_perim_std  = _ms(div_r_perim_all)
    div_phi_perim_mean,div_phi_perim_std= _ms(div_phi_perim_all)
    div_r_area_mean,   div_r_area_std   = _ms(div_r_area_all)
    div_phi_area_mean, div_phi_area_std = _ms(div_phi_area_all)
    div_r_diss_mean,   div_r_diss_std   = _ms(div_r_diss_all)
    div_phi_diss_mean, div_phi_diss_std = _ms(div_phi_diss_all)
    div_r_comb_mean,   div_r_comb_std   = _ms(div_r_comb_all)
    div_phi_comb_mean, div_phi_comb_std = _ms(div_phi_comb_all)

    print("Averaging complete.")

    # ── numerical-differentiation uncertainty for divergences ─────────────
    pos_T = bin_pos_ref

    def _div_combined_std(div_std, T_rr_mean, T_rphi_mean):
        num_r   = _numerical_diff_err(T_rr_mean,   pos_T)
        num_phi = _numerical_diff_err(T_rphi_mean, pos_T)
        return (np.sqrt(div_std[0]**2 + num_r**2),
                np.sqrt(div_std[1]**2 + num_phi**2))

    div_r_total_cstd, div_phi_total_cstd = _div_combined_std(
        (div_r_total_std,  div_phi_total_std),  T_rr_total_mean,  T_rphi_total_mean)
    div_r_chiral_cstd, div_phi_chiral_cstd = _div_combined_std(
        (div_r_chiral_std, div_phi_chiral_std), T_rr_chiral_mean, T_rphi_chiral_mean)
    div_r_perim_cstd, div_phi_perim_cstd = _div_combined_std(
        (div_r_perim_std,  div_phi_perim_std),  T_rr_perim_mean,  T_rphi_perim_mean)
    div_r_area_cstd, div_phi_area_cstd = _div_combined_std(
        (div_r_area_std,   div_phi_area_std),   T_rr_area_mean,   T_rphi_area_mean)

    # dissipative stress uses its own grid (same bin_pos_ref by construction)
    pos_T_diss = bin_pos_ref[:T_rr_diss_mean.shape[0]]
    _num_r_d   = _numerical_diff_err(T_rr_diss_mean,   pos_T_diss)
    _num_phi_d = _numerical_diff_err(T_rphi_diss_mean, pos_T_diss)
    div_r_diss_cstd   = np.sqrt(div_r_diss_std**2   + _num_r_d**2)
    div_phi_diss_cstd = np.sqrt(div_phi_diss_std**2 + _num_phi_d**2)

    # combined stress numerical-diff uncertainty from T_rr_total + T_rr_diss
    n_comb = min(T_rr_total_mean.shape[0], T_rr_diss_mean.shape[0])
    T_rr_comb_mean   = T_rr_total_mean[:n_comb]   + T_rr_diss_mean[:n_comb]
    T_rphi_comb_mean = T_rphi_total_mean[:n_comb] + T_rphi_diss_mean[:n_comb]
    pos_T_comb = bin_pos_ref[:n_comb]
    _num_r_c   = _numerical_diff_err(T_rr_comb_mean,   pos_T_comb)
    _num_phi_c = _numerical_diff_err(T_rphi_comb_mean, pos_T_comb)
    div_r_comb_cstd   = np.sqrt(div_r_comb_std**2   + _num_r_c**2)
    div_phi_comb_cstd = np.sqrt(div_phi_comb_std**2 + _num_phi_c**2)

    # ── divergence of the ensemble-mean stress tensor ─────────────────────
    (div_r_total_fm,   div_phi_total_fm,
     sig_r_total_fm,   sig_phi_total_fm)   = _div_from_mean(
        T_rr_total_mean,   T_phiphi_total_mean,
        T_rphi_total_mean, T_phir_total_mean,
        T_rr_total_std,    T_phiphi_total_std,
        T_rphi_total_std,  T_phir_total_std, pos_T)

    (div_r_chiral_fm,  div_phi_chiral_fm,
     sig_r_chiral_fm,  sig_phi_chiral_fm)  = _div_from_mean(
        T_rr_chiral_mean,   T_phiphi_chiral_mean,
        T_rphi_chiral_mean, T_phir_chiral_mean,
        T_rr_chiral_std,    T_phiphi_chiral_std,
        T_rphi_chiral_std,  T_phir_chiral_std, pos_T)

    (div_r_perim_fm,   div_phi_perim_fm,
     sig_r_perim_fm,   sig_phi_perim_fm)   = _div_from_mean(
        T_rr_perim_mean,   T_phiphi_perim_mean,
        T_rphi_perim_mean, T_phir_perim_mean,
        T_rr_perim_std,    T_phiphi_perim_std,
        T_rphi_perim_std,  T_phir_perim_std, pos_T)

    (div_r_area_fm,    div_phi_area_fm,
     sig_r_area_fm,    sig_phi_area_fm)    = _div_from_mean(
        T_rr_area_mean,   T_phiphi_area_mean,
        T_rphi_area_mean, T_phir_area_mean,
        T_rr_area_std,    T_phiphi_area_std,
        T_rphi_area_std,  T_phir_area_std, pos_T)

    print("div(mean stress) computed for all force types.")

    # ── radial grids for plotting ─────────────────────────────────────────
    r_div = r_grad_ref[:div_r_total_mean.shape[0]]

    # ── hypothesis test 1: div - rho_v*vF = -k*zetaW*Lap(vF) ─────────────
    v_fr_sem   = v_fr_std   / np.sqrt(N)
    v_fphi_sem = v_fphi_std / np.sqrt(N)

    print(f'\n=== Hypothesis 1: div(stress) - rho_v*vF = -k*zetaW*Lap(vF) ===')
    print(f'   zetaW = {p["zettaW"]},  R_fit = {R_fit:.3f}')

    kw = dict(cbins_LapV_ref=cbins_LapV_ref,
              cbins_vd_ref=cbins_vd_ref,
              cbins_ref=cbins_ref,
              R_fit=R_fit)

    (r_c_r, LHS_r, sig_LHS_r, X_r, sig_X_r, k_r, sig_k_r) = _hypothesis_fit(
        div_r_total_mean, div_r_total_cstd,
        vd_mean,          unc_vd_combined,
        v_fr_mean,        v_fr_sem,
        LapVr_mean,       sig_LapVr_combined,
        p['zettaW'], r_div, label='Radial (total force)  ', **kw)

    (r_c_phi, LHS_phi, sig_LHS_phi, X_phi, sig_X_phi, k_phi, sig_k_phi) = _hypothesis_fit(
        div_phi_total_mean, div_phi_total_cstd,
        vd_mean,            unc_vd_combined,
        v_fphi_mean,        v_fphi_sem,
        LapVphi_mean,       sig_LapVphi_combined,
        p['zettaW'], r_div, label='Azimuthal (total force)', **kw)

    print(f'\nk (radial)    = {k_r:.4f} +/- {sig_k_r:.4f}')
    print(f'k (azimuthal) = {k_phi:.4f} +/- {sig_k_phi:.4f}')

    # ── hypothesis test 2: div(T_total + T_diss) = k * vF ────────────────
    r_div_comb  = r_grad_ref[:div_r_comb_mean.shape[0]]
    r_div_total = r_grad_ref[:div_r_total_mean.shape[0]]
    r_div_diss  = r_grad_ref[:div_r_diss_mean.shape[0]]

    print(f'\n=== Hypothesis 2: div(T_total + T_diss) = k * vF ===')
    print(f'   R_fit = {R_fit:.3f}')

    (r_c2_r, div2_r, sig_div2_r, vF2_r, sig_vF2_r, k2_r, sig_k2_r) = \
        _combined_stress_fit(
            div_r_comb_mean,
            div_r_total_cstd,  r_div_total,
            div_r_diss_cstd,   r_div_diss,
            v_fr_mean,         v_fr_sem,
            r_div_comb, cbins_ref, R_fit,
            label='Radial (total+diss)   ')

    (r_c2_phi, div2_phi, sig_div2_phi, vF2_phi, sig_vF2_phi, k2_phi, sig_k2_phi) = \
        _combined_stress_fit(
            div_phi_comb_mean,
            div_phi_total_cstd, r_div_total,
            div_phi_diss_cstd,  r_div_diss,
            v_fphi_mean,        v_fphi_sem,
            r_div_comb, cbins_ref, R_fit,
            label='Azimuthal (total+diss)')

    print(f'\nk radial    = {k2_r:.4f} +/- {sig_k2_r:.4f}')
    print(f'k azimuthal = {k2_phi:.4f} +/- {sig_k2_phi:.4f}')

    # ═══════════════════════════════════════════════════════════════════════
    # Plots
    # ═══════════════════════════════════════════════════════════════════════
    out_dir = os.path.join(_SCRIPT_DIR, 'wet_force_balance_out')
    os.makedirs(out_dir, exist_ok=True)

    def _save(name):
        path = os.path.join(out_dir, name)
        plt.savefig(path, dpi=150, bbox_inches='tight')
        print(f'Saved {path}')

    # ── Plot 1: ensemble-averaged div(stress) vs vertex forces ────────────
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    for ax, div_r, std_r, div_phi, std_phi, vF_mean, vF_std, vF_lbl, comp in [
        (axes[0], div_r_total_mean,  div_r_total_cstd,
                  div_r_chiral_mean, div_r_chiral_cstd,
                  v_fr_mean, v_fr_std, r'$\langle v_{F,r} \rangle$', 'Radial'),
        (axes[1], div_phi_total_mean,  div_phi_total_cstd,
                  div_phi_chiral_mean, div_phi_chiral_cstd,
                  v_fphi_mean, v_fphi_std, r'$\langle v_{F,\phi} \rangle$', 'Azimuthal'),
    ]:
        _shaded(ax, r_div, div_r,   std_r,   'b-o',  'Total Force',  'b')
        _shaded(ax, r_div, div_phi, std_phi, 'm-.d', 'Chiral Force', 'm')
        _shaded(ax, cbins_ref, vF_mean, vF_std, 'k-s', vF_lbl, 'k')
        ax.set_xlim(0, r_a); ax.set_xlabel(r'$r$'); ax.legend(); ax.grid(True)
        ax.set_title(f'{comp}: div(σ) vs vertex force')

    plt.suptitle(
        f'Ensemble-averaged Kirkwood stress divergence vs vertex forces\n'
        f'(N={N} seeds, P0={p["P0"]}, C1={p["C1"]}; band = ensemble std ⊕ num-diff unc.)',
        fontsize=11)
    plt.tight_layout()
    _save('01_div_stress_vs_vF.pdf')
    plt.show()

    # ── Plot 2: ⟨∇·σ⟩ vs ∇·⟨σ⟩ for each force type (4 × 2) ──────────────
    force_specs = [
        ('Total Force',     'b',
         div_r_total_mean,  div_r_total_cstd,  div_phi_total_mean,  div_phi_total_cstd,
         div_r_total_fm,    sig_r_total_fm,    div_phi_total_fm,    sig_phi_total_fm),
        ('Chiral Force',    'm',
         div_r_chiral_mean, div_r_chiral_cstd, div_phi_chiral_mean, div_phi_chiral_cstd,
         div_r_chiral_fm,   sig_r_chiral_fm,   div_phi_chiral_fm,   sig_phi_chiral_fm),
        ('Perimeter Force', 'g',
         div_r_perim_mean,  div_r_perim_cstd,  div_phi_perim_mean,  div_phi_perim_cstd,
         div_r_perim_fm,    sig_r_perim_fm,    div_phi_perim_fm,    sig_phi_perim_fm),
        ('Area Force',      'c',
         div_r_area_mean,   div_r_area_cstd,   div_phi_area_mean,   div_phi_area_cstd,
         div_r_area_fm,     sig_r_area_fm,     div_phi_area_fm,     sig_phi_area_fm),
    ]

    fig, axes = plt.subplots(4, 2, figsize=(13, 16))
    for row, (lbl, col,
              dr_ea, std_r_ea, dp_ea, std_p_ea,
              dr_fm, std_r_fm, dp_fm, std_p_fm) in enumerate(force_specs):
        for c, (div_ea, s_ea, div_fm, s_fm, vF_m, vF_s, vF_lbl) in enumerate([
            (dr_ea, std_r_ea, dr_fm, std_r_fm, v_fr_mean,   v_fr_std,   r'$\langle v_{F,r}\rangle$'),
            (dp_ea, std_p_ea, dp_fm, std_p_fm, v_fphi_mean, v_fphi_std, r'$\langle v_{F,\phi}\rangle$'),
        ]):
            ax = axes[row, c]
            ax.plot(r_div, div_ea, '-o', ms=3, color=col,
                    label=r'$\langle\nabla\cdot\sigma\rangle$')
            ax.fill_between(r_div, div_ea - s_ea, div_ea + s_ea, alpha=0.25, color=col)
            ax.plot(pos_T, div_fm, '--s', ms=3, color=col,
                    label=r'$\nabla\cdot\langle\sigma\rangle$')
            ax.fill_between(pos_T, div_fm - s_fm, div_fm + s_fm, alpha=0.10, color=col)
            _shaded(ax, cbins_ref, vF_m, vF_s, 'k-^', vF_lbl, 'k')
            ax.set_xlim(0, r_a); ax.set_xlabel(r'$r$'); ax.legend(fontsize=9); ax.grid(True)
            ax.set_title(f'{lbl} — {"r" if c==0 else "φ"} component')

    plt.suptitle(
        r'$\langle\nabla\cdot\sigma\rangle$ vs $\nabla\cdot\langle\sigma\rangle$'
        f'\n(P0={p["P0"]}, C1={p["C1"]}, N={N} seeds)', fontsize=11)
    plt.tight_layout()
    _save('02_ensemble_avg_vs_div_of_mean.pdf')
    plt.show()

    # ── Plot 3: vertex density and Lap(v) ─────────────────────────────────
    fig, axes = plt.subplots(1, 3, figsize=(18, 5))
    r_vd  = cbins_vd_ref[:vd_mean.shape[0]]
    r_lap = cbins_LapV_ref[:LapVr_mean.shape[0]]
    _shaded(axes[0], r_vd,  vd_mean,       unc_vd_combined,     'b-o', r'$\langle\rho_v\rangle$',     'b')
    _shaded(axes[1], r_lap, LapVr_mean,   sig_LapVr_combined,   'r-s', r'$\langle\nabla^2 v_r\rangle$','r')
    _shaded(axes[2], r_lap, LapVphi_mean, sig_LapVphi_combined, 'g-^', r'$\langle\nabla^2 v_\phi\rangle$','g')
    for ax, ylabel, title in [
        (axes[0], r'Vertex density $\rho_v$', 'Ensemble-averaged vertex density'),
        (axes[1], r'$\nabla^2 v_r$',          r'Ensemble-averaged $\nabla^2 v_r$'),
        (axes[2], r'$\nabla^2 v_\phi$',        r'Ensemble-averaged $\nabla^2 v_\phi$'),
    ]:
        ax.set_xlim(0, r_a); ax.set_xlabel(r'$r$'); ax.set_ylabel(ylabel)
        ax.set_title(title); ax.legend(); ax.grid(True)
    plt.suptitle(f'Ensemble-averaged vertex density and ∇²v  (N={N}, P0={p["P0"]}, C1={p["C1"]})')
    plt.tight_layout()
    _save('03_vertex_density_and_LapV.pdf')
    plt.show()

    # ── Plot 4: hypothesis 1 — div - rho_v*vF vs -k*zetaW*Lap(vF) ────────
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    for ax, r_c, LHS, sLHS, X, sX, k, sk, comp in [
        (axes[0], r_c_r,   LHS_r,   sig_LHS_r,   X_r,   sig_X_r,   k_r,   sig_k_r,   'Radial'),
        (axes[1], r_c_phi, LHS_phi, sig_LHS_phi, X_phi, sig_X_phi, k_phi, sig_k_phi, 'Azimuthal'),
    ]:
        ax.errorbar(r_c, LHS,  yerr=sLHS, fmt='b-o', ms=4, capsize=3,
                    label=r'$\nabla\cdot\sigma - \rho_v\,v_F$  (LHS)')
        ax.errorbar(r_c, k*X, yerr=np.abs(k)*sX, fmt='r--s', ms=4, capsize=3,
                    label=rf'$k(-\zeta_W\nabla^2 v_F)$,  $k={k:.3f}\pm{sk:.3f}$')
        ax.fill_between(r_c, (k-sk)*X, (k+sk)*X, alpha=0.15, color='r')
        ax.axvline(R_fit, color='gray', lw=1.2, ls='--',
                   label=f'fit cutoff $R_{{\\rm fit}}={R_fit:.1f}$')
        ax.axhline(0, color='k', lw=0.7, ls=':')
        ax.set_xlim(0, r_a); ax.set_xlabel(r'$r$')
        ax.set_ylabel(f'{comp} component')
        ax.set_title(f'{comp}: fit $r<R_{{\\rm fit}}$ only')
        ax.legend(fontsize=9); ax.grid(True)
    plt.suptitle(
        r'Ensemble-avg: $\nabla\cdot\sigma - \rho_v\,v_F = -k\,\zeta_W\,\nabla^2 v_F$'
        f'\n(N={N} seeds, P0={p["P0"]}, C1={p["C1"]},'
        f' $\\zeta_W$={p["zettaW"]}, fit $r<R_{{\\rm fit}}={R_fit:.1f}$)',
        fontsize=12)
    plt.tight_layout()
    _save('04_hypothesis1_div_stress_LapV.pdf')
    plt.show()

    # ── Plot 5: hypothesis 2 — div(T_total + T_diss) vs k*vF ─────────────
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    for ax, r_c, div_c, sd_c, vF_c, svF_c, k, sk, comp in [
        (axes[0], r_c2_r,   div2_r,   sig_div2_r,   vF2_r,   sig_vF2_r,   k2_r,   sig_k2_r,   'Radial'),
        (axes[1], r_c2_phi, div2_phi, sig_div2_phi, vF2_phi, sig_vF2_phi, k2_phi, sig_k2_phi, 'Azimuthal'),
    ]:
        ax.errorbar(r_c, div_c,  yerr=sd_c,          fmt='b-o',  ms=4, capsize=3,
                    label=r'$\nabla\cdot(\sigma^{\rm tot}+\sigma^{\rm diss})$')
        ax.errorbar(r_c, k*vF_c, yerr=np.abs(k)*svF_c, fmt='r--s', ms=4, capsize=3,
                    label=rf'$k\,v_F$,  $k={k:.3f}\pm{sk:.3f}$')
        ax.fill_between(r_c, (k-sk)*vF_c, (k+sk)*vF_c, alpha=0.15, color='r')
        ax.axvline(R_fit, color='gray', lw=1.2, ls='--',
                   label=f'fit cutoff $R_{{\\rm fit}}={R_fit:.1f}$')
        ax.axhline(0, color='k', lw=0.7, ls=':')
        ax.set_xlim(0, r_a); ax.set_xlabel(r'$r$')
        ax.set_ylabel(f'{comp} component')
        ax.set_title(f'{comp}: $\\nabla\\cdot(\\sigma^{{\\rm tot}}+\\sigma^{{\\rm diss}}) = k\\,v_F$')
        ax.legend(fontsize=9); ax.grid(True)
    plt.suptitle(
        r'Ensemble-avg: $\nabla\cdot(\sigma^{\rm tot}+\sigma^{\rm diss}) = k\,v_F$'
        f'\n(N={N} seeds, P0={p["P0"]}, C1={p["C1"]}, fit $r<R_{{\\rm fit}}={R_fit:.1f}$)',
        fontsize=12)
    plt.tight_layout()
    _save('05_hypothesis2_div_total_diss_vs_vF.pdf')
    plt.show()

    # ── return all key results ─────────────────────────────────────────────
    return dict(
        # grids
        r_div=r_div, r_div_comb=r_div_comb, bin_pos_ref=bin_pos_ref,
        cbins_ref=cbins_ref, cbins_vd_ref=cbins_vd_ref, cbins_LapV_ref=cbins_LapV_ref,
        R=R, r_a=r_a, R_fit=R_fit,
        # vertex forces
        v_fr_mean=v_fr_mean, v_fr_std=v_fr_std,
        v_fphi_mean=v_fphi_mean, v_fphi_std=v_fphi_std,
        # vertex density
        vd_mean=vd_mean, unc_vd_combined=unc_vd_combined,
        # Lap(v)
        LapVr_mean=LapVr_mean, sig_LapVr_combined=sig_LapVr_combined,
        LapVphi_mean=LapVphi_mean, sig_LapVphi_combined=sig_LapVphi_combined,
        # divergences (ensemble avg of per-iT div)
        div_r_total_mean=div_r_total_mean,   div_r_total_cstd=div_r_total_cstd,
        div_phi_total_mean=div_phi_total_mean, div_phi_total_cstd=div_phi_total_cstd,
        div_r_diss_mean=div_r_diss_mean,     div_r_diss_cstd=div_r_diss_cstd,
        div_phi_diss_mean=div_phi_diss_mean, div_phi_diss_cstd=div_phi_diss_cstd,
        div_r_comb_mean=div_r_comb_mean,     div_r_comb_cstd=div_r_comb_cstd,
        div_phi_comb_mean=div_phi_comb_mean, div_phi_comb_cstd=div_phi_comb_cstd,
        # hypothesis 1 fit
        k_r=k_r, sig_k_r=sig_k_r, k_phi=k_phi, sig_k_phi=sig_k_phi,
        # hypothesis 2 fit
        k2_r=k2_r, sig_k2_r=sig_k2_r, k2_phi=k2_phi, sig_k2_phi=sig_k2_phi,
    )


# ═══════════════════════════════════════════════════════════════════════════
if __name__ == '__main__':
    results = run_ensemble_average_analysis()
