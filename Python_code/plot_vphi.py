import numpy as np
import matplotlib.pyplot as plt
from scipy.special import i0, i1, k0, k1

# ── Parameters ────────────────────────────────────────────────────────────────
eps0 = 0.03   # ε₀
r0   = 7.5    # r₀

# List of (γ, η) pairs to plot
pairs = [
    (1.0, 1.0),
    (1.0, 2.0),
    (2.0, 1.0),
]

# ── Coefficients ──────────────────────────────────────────────────────────────
def compute_AB(gamma, eta, eps0=eps0, r0=r0):
    alpha = np.sqrt(gamma / eta)   # √(γ/η)

    x1 = alpha * r0                # √(γ/η) · r₀  (argument for I1/K1 Wronskian)
    x2 = r0 / alpha                # √(η/γ) · r₀  (argument in A and B numerators)

    I1  = i1(x1)
    K1  = k1(x1)
    # I1'(x) = I0(x) - I1(x)/x,   K1'(x) = -K0(x) - K1(x)/x
    I1p = i0(x1) - I1 / x1
    K1p = -k0(x1) - K1 / x1

    D = I1 * K1p - I1p * K1       # Wronskian denominator

    A = -eps0 * k1(x2) / (np.sqrt(gamma * eta) * D)
    B = -eps0 * i1(x2) / (np.sqrt(gamma * eta) * D)
    return A, B, alpha

# ── v_φ(r) ────────────────────────────────────────────────────────────────────
def v_phi(r, A, B, alpha, r0=r0):
    r = np.asarray(r, dtype=float)
    out = np.full_like(r, np.nan)
    inner = r < r0
    outer = r > r0
    out[inner] = A * i1(alpha * r[inner])
    out[outer] = B * k1(alpha * r[outer])
    return out

# ── Plot ──────────────────────────────────────────────────────────────────────
r = np.linspace(0.05 * r0, 4.0 * r0, 2000)

fig, ax = plt.subplots(figsize=(7, 5))
for gamma, eta in pairs:
    A, B, alpha = compute_AB(gamma, eta)
    ax.plot(r, v_phi(r, A, B, alpha), label=rf'$\gamma={gamma:g},\,\eta={eta:g}$')

ax.axvline(r0, color='gray', ls='--', lw=0.8, label=rf'$r_0={r0}$')
ax.axhline(0,  color='k',    ls=':',  lw=0.6)
ax.set_xlabel(r'$r$')
ax.set_ylabel(r'$v_\phi(r)$')
ax.set_title(r'Azimuthal velocity profile $v_\phi(r)$')
ax.legend()
ax.grid(alpha=0.3)
plt.tight_layout()
plt.show()
