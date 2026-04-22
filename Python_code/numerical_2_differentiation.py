"""
Second-order numerical differentiation function

Computes second-order numerical derivative using central finite differences
"""

import numpy as np


def numerical_2_differentiation(unused_param, V, R):
    """
    Compute numerical second derivative

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
        Numerical second derivative d²V/dR²
    r : ndarray
        Corresponding R positions
    """
    # Initialize result arrays
    D1 = len(V)
    result = np.zeros(D1 - 2)
    r = np.zeros(D1 - 2)

    # Loop over interior points (need neighbors on both sides)
    for ind in range(1, D1 - 1):
        # Compute second derivative using central difference
        delta = R[ind + 1] - R[ind]
        result[ind - 1] = (V[ind - 1] - 2 * V[ind] + V[ind + 1]) / (delta ** 2)
        r[ind - 1] = R[ind + 1]

    return result, r
