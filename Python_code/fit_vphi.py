"""
fit_vphi.py
-----------
Fit the analytic v_phi(r; gamma, eta, r0) profile to ensemble-averaged
simulation data for a list of zettaW values, then plot the fitted parameters
(r0, gamma, eta) as a function of zettaW.

Analytic model (from Kirkwood_stress.ipynb):
    k  = sqrt(gamma / eta)
    v_phi(r) = (epsilon0 * r0 / eta) * K1(k*r0) * I1(k*r)   for r <= r0
               (epsilon0 * r0 / eta) * I1(k*r0) * K1(k*r)   for r  > r0

Because the data are bin-averaged, the model is compared against the
area-weighted radial average over each bin:
    <v_phi>_bin = int_{r1}^{r2} v_phi(r) r dr / int_{r1}^{r2} r dr

Uncertainties are estimated as the standard error of the mean (SEM) of the
individual vertex azimuthal velocities within each radial bin.

Usage
-----
Run from the Python_code directory:
    python fit_vphi.py
"""

import os, sys
import numpy as np
import matplotlib.pyplot as plt
from scipy.special import iv, kv
from scipy.optimize import curve_fit
from scipy.integrate import quad

# ── path so local modules are importable ───────────────────────────────────────
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from read_file_folder_X import read_file_folder_X
from read_file_folder_cVertices import read_file_folder_cVertices
from read_file_folder_cell_bin_categories import read_file_folder_cell_bin_categories
from read_file_folder import read_file_folder

# ── simulation parameters ──────────────────────────────────────────────────────
Nx           = 60
gammaW       = 1
disorder     = 0.02
kP           = 0.01
dynamics     = 1
transition_type = 0
P0_list      = [3, 3.1, 3.2, 3.3, 3.4, 3.5, 3.6, 3.7, 3.8, 3.9, 4, 4.1, 4.2]
C1           = 0.03          # epsilon_0 in the analytic solution
dw           = 50
di           = 50
tMAX         = 1000
seed         = 1
iT_list      = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
Lb           = 1             # bin-width multiplier used in vF_ave_func

epsilon0 = float(C1)         # epsilon_0 in analytic solution

# zettaW_list = [1, 2, 4, 8, 16, 32, 64, 128]
zettaW_list = [1e-1,1, 2, 4, 8, 16, 32, 64, 128]

# ══════════════════════════════════════════════════════════════════════════════
# 1.  Data loading with per-bin statistics
# ══════════════════════════════════════════════════════════════════════════════

def load_vphi_with_uncertainty(filepath, time, Lb):
    """
    Load azimuthal velocity data and return per-bin mean, SEM, and bin edges.

    Returns
    -------
    cbins      : (nbins,)   bin edges
    v_fphi_mean: (nbins-1,) mean azimuthal velocity per bin
    v_fphi_sem : (nbins-1,) standard error of the mean per bin
    counts     : (nbins-1,) number of vertices per bin
    perioY     : float      system box height
    """
    # ── raw velocity file ──────────────────────────────────────────────────
    data  = np.loadtxt(os.path.join(filepath, f'v_xy_{time}.dat'), comments='#')
    v_fx  = data[:, 4]
    v_fy  = data[:, 5]
    v_x   = data[:, 2]
    v_y   = data[:, 3]
    Nv    = len(v_fx)

    # ── cell centres & chiral-cluster centre ──────────────────────────────
    c_v_xpos_list = read_file_folder_X(filepath, f'X_{time}.dat')
    c_v_ypos_list = read_file_folder_X(filepath, f'Y_{time}.dat')
    cVertices     = read_file_folder_cVertices(filepath, f'c_vertices_{time}.dat')
    cActivity     = read_file_folder_cell_bin_categories(filepath, f'c_activity_{time}.dat')
    _, perioY     = read_file_folder(filepath, f'perio_{time}.dat')
    perioY        = float(np.asarray(perioY).flat[0])
    Nc            = len(cVertices)

    xpos_center = np.zeros(Nc)
    ypos_center = np.zeros(Nc)
    xpos_active, ypos_active = [], []

    for cid in range(Nc):
        row    = cVertices[cid]
        xlist  = c_v_xpos_list[cid]
        ylist  = c_v_ypos_list[cid]
        Nv_c   = int(row[1])
        xpos_center[cid] = np.mean(xlist[:Nv_c])
        ypos_center[cid] = np.mean(ylist[:Nv_c])
        if cActivity[cid] == 1:
            xpos_active.append(xpos_center[cid])
            ypos_active.append(ypos_center[cid])

    cx = np.mean(xpos_active)
    cy = np.mean(ypos_active)

    dx = v_x - cx
    dy = v_y - cy
    r  = np.sqrt(dx**2 + dy**2)

    # azimuthal component of each vertex's force-velocity
    f_phi = (-v_fx * dy + v_fy * dx) / r

    # ── radial bins ───────────────────────────────────────────────────────
    bin_step = Lb * perioY / np.sqrt(Nc)
    cbins    = np.arange(0, perioY + bin_step, bin_step)
    nbins    = len(cbins)

    # assign vertices to bins
    bin_ids = np.digitize(r, cbins) - 1          # 0-based bin index
    bin_ids = np.clip(bin_ids, 0, nbins - 2)     # stay within [0, nbins-2]

    # per-bin statistics
    v_fphi_mean = np.zeros(nbins - 1)
    v_fphi_sem  = np.zeros(nbins - 1)
    counts      = np.zeros(nbins - 1, dtype=int)

    for b in range(nbins - 1):
        mask = bin_ids == b
        n    = mask.sum()
        counts[b] = n
        if n > 0:
            vals = f_phi[mask]
            v_fphi_mean[b] = vals.mean()
            v_fphi_sem[b]  = vals.std(ddof=1) / np.sqrt(n) if n > 1 else np.nan

    return cbins, v_fphi_mean, v_fphi_sem, counts, perioY


def load_vphi_ensemble(P0, zettaW):
    """
    Load v_fphi for all available iT in iT_list, ensemble-average across trials.

    For each iT the per-bin mean is computed; the ensemble mean and SEM are
    then taken across the N_iT realisations:
        ensemble_mean[b] = mean over iT of v_fphi_mean_iT[b]
        ensemble_sem[b]  = std over iT of v_fphi_mean_iT[b] / sqrt(N_iT)

    Returns
    -------
    cbins          : (nbins,)   bin edges (from first valid iT)
    ensemble_mean  : (nbins-1,) ensemble-averaged bin means
    ensemble_sem   : (nbins-1,) SEM across iT realisations
    mean_counts    : (nbins-1,) mean vertex count per bin across iT
    perioY         : float      system box height
    n_iT_used      : int        number of iT realisations that were loaded
    """
    means_list  = []
    counts_list = []
    cbins_ref   = None
    perioY_ref  = None

    for iT in iT_list:
        fp = (f'/mnt/users/wangg/VertexModelChiralWetEdges/output/'
              f'out_Nx_{Nx}_kA_0.5_kP_{kP}_P0_{P0}_disorder_{disorder}_'
              f'tMAX_{tMAX}_dw_{dw}_di_{di}_ca1_{C1}_ca2_0_iC_0_iP_0.25_'
              f'iS_{seed}_iT_{iT}_dynamics_{dynamics}_gammaW_{gammaW}_'
              f'zetaW_{zettaW}_tT_{transition_type}/')
        if not os.path.isdir(fp):
            continue
        if not os.path.isfile(os.path.join(fp, f'v_xy_{tMAX}.dat')):
            continue

        cbins, v_fphi_mean, _, counts, perioY = \
            load_vphi_with_uncertainty(fp, tMAX, Lb)

        # align to common bin length using the reference from the first iT
        if cbins_ref is None:
            cbins_ref  = cbins
            perioY_ref = perioY
        n = min(len(v_fphi_mean), len(cbins_ref) - 1)
        means_list.append(v_fphi_mean[:n])
        counts_list.append(counts[:n])

    if not means_list:
        return None

    n_iT   = len(means_list)
    n_bins = min(len(m) for m in means_list)
    arr    = np.array([m[:n_bins] for m in means_list])   # (n_iT, n_bins)
    cnt    = np.array([c[:n_bins] for c in counts_list])  # (n_iT, n_bins)

    ensemble_mean = arr.mean(axis=0)
    ensemble_sem  = arr.std(axis=0, ddof=1) / np.sqrt(n_iT) if n_iT > 1 \
                    else np.full(n_bins, np.nan)
    mean_counts   = cnt.mean(axis=0)

    return cbins_ref[:n_bins + 1], ensemble_mean, ensemble_sem, \
           mean_counts.astype(int), perioY_ref, n_iT


# ══════════════════════════════════════════════════════════════════════════════
# 2.  Analytic model & area-weighted bin average
# ══════════════════════════════════════════════════════════════════════════════

def v_phi_analytic(r, gamma, eta, r0, epsilon0):
    """Point-wise analytic azimuthal velocity (scalar or 1-D array)."""
    k   = np.sqrt(gamma / eta)
    kr0 = k * r0
    A   = epsilon0 * r0 / eta * kv(1, kr0)
    B   = epsilon0 * r0 / eta * iv(1, kr0)
    r   = np.asarray(r, dtype=float)
    out = np.empty_like(r)
    inn = r <= r0
    out[inn]  = A * iv(1, k * r[inn])
    out[~inn] = B * kv(1, k * r[~inn])
    return out


def v_phi_bin_avg(r1, r2, gamma, eta, r0, epsilon0):
    """
    Area-weighted average of v_phi over [r1, r2]:
        <v_phi> = int_{r1}^{r2} v_phi(r)*r dr / int_{r1}^{r2} r dr
    Split integral at r0 if it falls inside the bin.
    """
    denom = 0.5 * (r2**2 - r1**2)
    if denom == 0:
        return 0.0

    def integrand(r):
        return v_phi_analytic(np.atleast_1d(r), gamma, eta, r0, epsilon0)[0] * r

    if r1 < r0 < r2:
        num1, _ = quad(integrand, r1, r0)
        num2, _ = quad(integrand, r0, r2)
        num = num1 + num2
    else:
        num, _ = quad(integrand, r1, r2)

    return num / denom


# ══════════════════════════════════════════════════════════════════════════════
# 3.  Fit one dataset
# ══════════════════════════════════════════════════════════════════════════════

def fit_one(P0, zettaW):
    """
    Ensemble-average v_fphi over iT_list, then fit gamma/eta/r0/epsilon0.
    """
    loaded = load_vphi_ensemble(P0, zettaW)
    if loaded is None:
        print(f"  [skip] no data found for zettaW={zettaW}, P0={P0}")
        return None

    cbins, v_fphi_mean, v_fphi_sem, counts, perioY, n_iT = loaded

    bin_centres = 0.5 * (cbins[:-1] + cbins[1:])
    r_fit_max   = 0.5 * Nx

    valid = ((counts >= 2) & np.isfinite(v_fphi_sem) & (v_fphi_sem > 0)
             & (bin_centres <= r_fit_max))

    r_data    = bin_centres[valid]
    y_data    = v_fphi_mean[valid]
    yerr_data = v_fphi_sem[valid]
    valid_indices = np.where(valid)[0]

    # initial guess for r0 from system geometry
    r0_init = 0.5 * 0.25 * perioY

    def model_for_fit(x_dummy, gamma_p, eta, r0, epsilon0_p):
        result = np.empty(len(valid_indices))
        for j, idx in enumerate(valid_indices):
            r1, r2 = cbins[idx], cbins[idx + 1]
            result[j] = v_phi_bin_avg(r1, r2, gamma_p, eta, r0, epsilon0_p)
        return result

    p0     = [float(gammaW), float(zettaW), r0_init, epsilon0]
    bounds = ([1e-3, 1e-3, 0.5,  1e-6],
              [1e4,  1e4,  r_fit_max, 1e2])

    try:
        popt, pcov = curve_fit(
            model_for_fit,
            np.zeros(valid.sum()),
            y_data,
            p0=p0,
            sigma=yerr_data,
            absolute_sigma=True,
            bounds=bounds,
            maxfev=10000,
        )
        perr = np.sqrt(np.diag(pcov))
        gamma_fit, eta_fit, r0_fit, epsilon0_fit = popt
        gamma_err, eta_err, r0_err, epsilon0_err = perr
        success = True
    except RuntimeError as e:
        print(f"  [zettaW={zettaW}] curve_fit did not converge: {e}")
        gamma_fit, eta_fit, r0_fit, epsilon0_fit = p0
        gamma_err = eta_err = r0_err = epsilon0_err = np.nan
        success = False

    print(f"  zettaW={zettaW:9g}  (N_iT={n_iT}):  "
          f"γ={gamma_fit:.3g}±{gamma_err:.1g}  "
          f"η={eta_fit:.3g}±{eta_err:.1g}  "
          f"r₀={r0_fit:.3g}±{r0_err:.1g}  "
          f"ε₀={epsilon0_fit:.3g}±{epsilon0_err:.1g}")

    # bin-averaged model on valid bins (for profile plot)
    v_model_bins = model_for_fit(None, gamma_fit, eta_fit, r0_fit, epsilon0_fit)

    return dict(
        zettaW=zettaW, n_iT=n_iT,
        cbins=cbins, bin_centres=bin_centres,
        r_data=r_data, y_data=y_data, yerr_data=yerr_data,
        valid=valid, valid_indices=valid_indices,
        perioY=perioY, r0_init=r0_init,
        gamma_fit=gamma_fit, eta_fit=eta_fit, r0_fit=r0_fit,
        epsilon0_fit=epsilon0_fit,
        gamma_err=gamma_err, eta_err=eta_err, r0_err=r0_err,
        epsilon0_err=epsilon0_err,
        v_model_bins=v_model_bins, success=success,
    )



# ══════════════════════════════════════════════════════════════════════════════
# 4.  Run fits over all P0 × zettaW
# ══════════════════════════════════════════════════════════════════════════════

import matplotlib.colors as mcolors

# unified profile-plot colours (same meaning in every subplot)
C_DATA = 'steelblue'
C_FIT  = 'crimson'
C_BIN  = 'darkorange'

# one colour per P0 for the comparison plot
P0_colors = plt.cm.tab20(np.linspace(0, 1, len(P0_list)))

# all_results[P0] = list of result dicts for that P0
all_results = {}

for P0 in P0_list:
    print(f"\n── P0 = {P0} ──────────────────────────────────────")
    results = []
    for zw in zettaW_list:
        res = fit_one(P0, zw)
        if res is not None:
            results.append(res)
    all_results[P0] = results


# ══════════════════════════════════════════════════════════════════════════════
# 5.  Per-P0 profile plots (all zettaW overlaid on one axis)
# ══════════════════════════════════════════════════════════════════════════════

zw_global = np.array(zettaW_list, dtype=float)
cmap_zw   = plt.cm.viridis
norm_zw   = mcolors.LogNorm(vmin=zw_global.min(), vmax=zw_global.max())

for P0, results in all_results.items():
    if not results:
        continue
    fig_prof, ax = plt.subplots(figsize=(8, 5))

    for res in results:
        color = cmap_zw(norm_zw(res['zettaW']))
        cbins = res['cbins']
        r_theory = np.linspace(cbins[1] * 0.5, cbins[-1], 500)
        v_theory = v_phi_analytic(r_theory, res['gamma_fit'],
                                  res['eta_fit'], res['r0_fit'], res['epsilon0_fit'])
        ax.errorbar(res['r_data'], res['y_data'], yerr=res['yerr_data'],
                    fmt='o', ms=3, color=color, ecolor=color, alpha=0.5,
                    elinewidth=0.8, capsize=0)
        ax.plot(r_theory, v_theory, '-', color=color, lw=1.8,
                label=rf'$\zeta_W={res["zettaW"]:g}$')

    ax.axhline(0, color='k', ls=':', lw=0.6)
    ax.set_xlabel(r'$r$', fontsize=12)
    ax.set_ylabel(r'$v_\phi(r)$', fontsize=12)
    ax.set_title(rf'$v_\phi$ fits  ($P_0={P0}$, $\varepsilon_0={epsilon0}$, $\gamma_W={gammaW}$)',
                 fontsize=11)
    ax.legend(fontsize=8, ncol=2)
    ax.grid(alpha=0.3)
    sm = plt.cm.ScalarMappable(cmap=cmap_zw, norm=norm_zw)
    sm.set_array([])
    fig_prof.colorbar(sm, ax=ax, label=r'$\zeta_W$')
    fig_prof.tight_layout()
    outfile = (f'/mnt/users/wangg/VertexModelChiralWetEdges/output/'
               f'vphi_fit_profiles_P0_{P0}.pdf')
    fig_prof.savefig(outfile, dpi=150, bbox_inches='tight')
    print(f"Profile (overlaid) saved → {outfile}")
    plt.close(fig_prof)


# ══════════════════════════════════════════════════════════════════════════════
# 5b.  Per-P0 subplot grid (one panel per zettaW)
# ══════════════════════════════════════════════════════════════════════════════

for P0, results in all_results.items():
    if not results:
        continue
    n  = len(results)
    nc = min(4, n)
    nr = int(np.ceil(n / nc))

    fig_sub, axes_sub = plt.subplots(nr, nc, figsize=(4.5 * nc, 4 * nr))
    axes_sub = np.array(axes_sub).reshape(-1)

    for ax_s, res in zip(axes_sub, results):
        cbins    = res['cbins']
        r_theory = np.linspace(cbins[1] * 0.5, cbins[-1], 500)
        v_theory = v_phi_analytic(r_theory, res['gamma_fit'],
                                  res['eta_fit'], res['r0_fit'], res['epsilon0_fit'])

        ax_s.errorbar(res['r_data'], res['y_data'], yerr=res['yerr_data'],
                      fmt='o', ms=4, color=C_DATA, ecolor=C_DATA, alpha=0.7,
                      elinewidth=1.0, capsize=2, label='data ± SEM', zorder=3)
        ax_s.plot(r_theory, v_theory, '-', color=C_FIT, lw=2.0,
                  label='analytic fit', zorder=4)
        ax_s.plot(res['r_data'], res['v_model_bins'], 's', ms=5,
                  markerfacecolor='white', markeredgecolor=C_BIN,
                  markeredgewidth=1.5, label='bin avg (analytic)', zorder=5)
        ax_s.axvline(res['r0_fit'], color=C_FIT, ls='--', lw=0.9, alpha=0.6)
        ax_s.axhline(0, color='k', ls=':', lw=0.5)
        ax_s.set_facecolor('#f9f9f9')
        ax_s.set_title(rf'$\zeta_W={res["zettaW"]:g}$  ($N_{{iT}}={res["n_iT"]}$)  '
                       rf'$\gamma={res["gamma_fit"]:.2g}$, '
                       rf'$\eta={res["eta_fit"]:.2g}$, '
                       rf'$r_0={res["r0_fit"]:.2g}$, '
                       rf'$\varepsilon_0={res["epsilon0_fit"]:.2g}$',
                       fontsize=8.5)
        ax_s.set_xlabel(r'$r$', fontsize=9)
        ax_s.set_ylabel(r'$v_\phi$', fontsize=9)
        ax_s.tick_params(labelsize=7)
        ax_s.grid(alpha=0.3)

    axes_sub[0].legend(fontsize=7, loc='upper right')
    for ax_s in axes_sub[n:]:
        ax_s.set_visible(False)

    fig_sub.suptitle(rf'$v_\phi$ fits  ($P_0={P0}$, $\varepsilon_0={epsilon0}$, $\gamma_W={gammaW}$)',
                     fontsize=12)
    fig_sub.tight_layout()
    outfile = (f'/mnt/users/wangg/VertexModelChiralWetEdges/output/'
               f'vphi_fit_profiles_subplots_P0_{P0}.pdf')
    fig_sub.savefig(outfile, dpi=150, bbox_inches='tight')
    print(f"Profile (subplots)  saved → {outfile}")
    plt.close(fig_sub)


# ══════════════════════════════════════════════════════════════════════════════
# 6.  Per-P0 parameter vs zettaW plots
# ══════════════════════════════════════════════════════════════════════════════

for P0, results in all_results.items():
    if not results:
        continue
    zw_arr        = np.array([r['zettaW']       for r in results])
    r0_arr        = np.array([r['r0_fit']       for r in results])
    r0_earr       = np.array([r['r0_err']       for r in results])
    gamma_arr     = np.array([r['gamma_fit']    for r in results])
    gamma_earr    = np.array([r['gamma_err']    for r in results])
    eta_arr       = np.array([r['eta_fit']      for r in results])
    eta_earr      = np.array([r['eta_err']      for r in results])
    epsilon0_arr  = np.array([r['epsilon0_fit'] for r in results])
    epsilon0_earr = np.array([r['epsilon0_err'] for r in results])

    fig_par, axes_par = plt.subplots(4, 1, figsize=(7, 12), sharex=True)
    for ax_p, (vals, errs, label, color) in zip(axes_par, [
            (r0_arr,       r0_earr,       r'$r_0$',           'teal'),
            (gamma_arr,    gamma_earr,    r'$\gamma$',        'steelblue'),
            (eta_arr,      eta_earr,      r'$\eta$',          'crimson'),
            (epsilon0_arr, epsilon0_earr, r'$\varepsilon_0$', 'darkorchid'),
    ]):
        ax_p.errorbar(zw_arr, vals, yerr=errs, fmt='o-', ms=6,
                      color=color, ecolor=color, capsize=4, lw=1.5)
        ax_p.set_ylabel(label, fontsize=13)
        ax_p.set_xscale('log')
        ax_p.grid(alpha=0.3, which='both')
        ax_p.tick_params(labelsize=10)
    axes_par[-1].set_xlabel(r'$\zeta_W$', fontsize=13)
    fig_par.suptitle(rf'Fitted parameters vs $\zeta_W$  ($P_0={P0}$, $\varepsilon_0={epsilon0}$, $\gamma_W={gammaW}$)',
                     fontsize=11)
    fig_par.tight_layout()
    outfile = (f'/mnt/users/wangg/VertexModelChiralWetEdges/output/'
               f'vphi_fit_params_P0_{P0}.pdf')
    fig_par.savefig(outfile, dpi=150, bbox_inches='tight')
    print(f"Params figure       saved → {outfile}")
    plt.close(fig_par)


# ══════════════════════════════════════════════════════════════════════════════
# 7.  Comparison: fitted parameters vs zettaW for all P0 on the same axes
# ══════════════════════════════════════════════════════════════════════════════

fig_cmp, axes_cmp = plt.subplots(4, 1, figsize=(8, 13), sharex=True)
param_keys = [
    ('r0_fit',       'r0_err',       r'$r_0$'),
    ('gamma_fit',    'gamma_err',    r'$\gamma$'),
    ('eta_fit',      'eta_err',      r'$\eta$'),
    ('epsilon0_fit', 'epsilon0_err', r'$\varepsilon_0$'),
]

for (P0, results), color in zip(all_results.items(), P0_colors):
    if not results:
        continue
    zw_arr = np.array([r['zettaW'] for r in results])
    for ax_c, (fkey, ekey, ylabel) in zip(axes_cmp, param_keys):
        vals = np.array([r[fkey] for r in results])
        errs = np.array([r[ekey] for r in results])
        ax_c.errorbar(zw_arr, vals, yerr=errs, fmt='o-', ms=6,
                      color=color, ecolor=color, capsize=4, lw=1.5,
                      label=rf'$P_0={P0}$')

for ax_c, (_, _, ylabel) in zip(axes_cmp, param_keys):
    ax_c.set_ylabel(ylabel, fontsize=13)
    ax_c.grid(alpha=0.3, which='both')
    ax_c.tick_params(labelsize=10)

axes_cmp[0].legend(fontsize=10, loc='best')
axes_cmp[-1].set_xlabel(r'$\zeta_W$', fontsize=13)
fig_cmp.suptitle(rf'Fitted parameters vs $\zeta_W$ for different $P_0$'
                 rf'  ($\varepsilon_0={epsilon0}$, $\gamma_W={gammaW}$)',
                 fontsize=11)
fig_cmp.tight_layout()
outfile_cmp = ('/mnt/users/wangg/VertexModelChiralWetEdges/output/'
               'vphi_fit_params_comparison.pdf')
fig_cmp.savefig(outfile_cmp, dpi=150, bbox_inches='tight')
print(f"\nComparison figure   saved → {outfile_cmp}")

plt.show()
