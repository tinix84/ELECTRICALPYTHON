###################################################################
#   SYSTEMSOLUTION.PY
#
#   A library of functions, constants and more that are
#   related to solving power flows in Electrical Engineering.
#
#   February 19, 2019
#
#   Written by Joe Stanley
#
#   Newton Raphson Method Solution derived using references from:
#   - SCIPY's Online Guides
#   - Documentation available at:
#     http://hplgit.github.io/Programming-for-Computations/pub/p4c/._p4c-solarized-Python031.html
#
#   Included Functions
#   - Newton Raphson System Solver:      NewtonRaphson
#
####################################################################

# Import Libraries
import numpy as np
from scipy.optimize import newton

# Define 
tint = "<class 'int'>"
tfloat = "<class 'float'>"

# Define Newton-Raphson Method
def NewtonRaphson(F, J, X0, eps=1e-4, mxiter=100):
    """
    Solve nonlinear system F=0 by Newton's method.
    J is the Jacobian of F. Both F and J must be functions of x.
    At input, x holds the start value. The iteration continues
    until ||F|| < eps.
    
    Required Arguments:
    -------------------
    F:     The Non-Linear System
    J:     The Jacobian of F
    X0:    The Initial Value (or initial guess)
    
    Optional Arguments:
    -------------------
    eps:     Epsilon - The error value, default=0.0001
    mxiter:  Maximum Iterations - The highest number of iterations allowed,
             default=100
    
    Returns:
    --------
    X0:                 The computed result
    
    Optional Returns:
    -----------------
    iteration_counter:  The number of iterations completed before returning
                        either due to solution being found, or max iterations
                        being surpassed.
    """
    # Test for one-variable inputs
    ftype = str(type(F(X0)))
    jtype = str(type(J(X0)))
    if(ftype == tint or ftype == tfloat): # System is size-1
        if((jtype != tint) and (jtype != tfloat)): # Jacobian isn't size-1
            raise ValueError("ERROR: The Jacobian isn't size-1.")
        return( newton( F, X0, J ) )
    
    # Test for valid argument sizes
    f0sz = len(F(X0))
    j0sz = len(J(X0))
    if(f0sz!=j0sz): # Size mismatch
        raise ValueError("ERROR: The arguments return arrays or lists"+
                        " of different sizes: f0="+str(f0sz)+"; j0="+str(j0sz))
    
    F_value = F(X0)
    F_norm = np.linalg.norm(F_value, ord=2)  # L2 norm of vector
    iteration_counter = 0
    while abs(F_norm) > eps and iteration_counter < mxiter:
        delta = np.linalg.solve(J(X0), -F_value)
        X0 = X0 + delta
        F_value = F(X0)
        F_norm = np.linalg.norm(F_value, ord=2)
        iteration_counter += 1

    # Here, either a solution is found, or too many iterations
    if abs(F_norm) > eps:
        iteration_counter = -1
    return(X0, iteration_counter)