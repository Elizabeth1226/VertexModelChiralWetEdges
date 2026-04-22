"""
Numerical differentiation function

Computes first-order numerical derivative using finite differences
"""

import numpy as np


def numerical_differentiation(unused_param, V, R):
    """
    Compute numerical first derivative

    Parameters:
    -----------
    unused_param : any
        Unused parameter (kept for MATLAB compatibility)
    V : array-like
        Values to differentiate
    R : array-like
        Positions/coordinates

    Returns:
    --------
    result : ndarray
        Numerical derivative dV/dR
    r : ndarray
        Corresponding R positions
    """
    # Initialize result arrays
    D1 = len(V)
    result = np.zeros(D1 - 1)
    r = np.zeros(D1 - 1)

    # First point: forward difference
    if D1 > 1:
        delta = R[1] - R[0]
        if abs(delta) > 1e-10:
            result[0] = (V[1] - V[0]) / delta
        r[0] = R[0]

    # Interior points: central difference
    for ind in range(1, D1 - 1):
        # Compute derivative using central difference
        delta = R[ind + 1] - R[ind - 1]
        if abs(delta) > 1e-10:
            result[ind] = (V[ind + 1] - V[ind - 1]) / delta
        r[ind] = R[ind]
        
    return result, r
