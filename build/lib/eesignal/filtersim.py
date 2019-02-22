#################################################################################
#   FILTERSIM.PY
#
#   This file contains a variety of functions and constants related to signals.
#   These items will commonly  be used in Electrical Engineering Applications.
#
#   February 21, 2019
#
#   Written by Joe Stanley
#   Special Thanks To and Code Support From:
#   Dr. Dennis Sullivan
#
#   Included Functions:
#   - Z-Domain Filter Simulator       zfiltersim
#################################################################################

import numpy as np
import matplotlib.pyplot as plt

def zfiltersim( fin, filter, NN=1000, title="",plotinput=True):
    """
    ZFILTERSIM Function
    
    Purpose:
    --------
    Given an input function and filter parameters (specified in
    the z-domain) this function will plot the input function over
    NN time-steps of an unspecified size (the step-size must be
    specified outside the scope of this function). This function
    will also plot the resultant (output) function over the NN
    time-steps after having been filtered by that filter which
    is specified.
    
    Required Arguments:
    -------------------
    fin:    The input function, must be callable with specified
            step-size.
    filter: The filter parameter set as described at the end of
            this help/comment section.
    
    Optional Arguments:
    NN:         The number of time-steps to be plotted; default=1000
    title:      The title presented on each plot; default=""
    plotinput:  An argument to control whether the input is plotted
                separately, default=True.
    
    Returns:
    --------
    NONE, this function generates plots.
    
    # NOTICE: ------------------------------------
    # The *filter* argument should be provided as:
    # [[ a11, a12, b10, b11, b12],
    #  [ a21, a22, b20, b21, b22],
    #  [           ...          ],
    #  [ an1, an2, bn0, bn1, bn2]]
    # Where each row corresponds to a 1- or 2-pole
    # filter with the structure as follows:
    #          b0 + b1z^-1 + b2z^-2
    # H(z) = -----------------------------
    #           1 - a1z^-1 - a2z^-2
    # --------------------------------------------
    """
    # Start with arrays set to zero
    x = np.zeros(NN)
    y = np.zeros(NN)
    
    # ----- The input  -----
    for k in range(NN):
        x[k] = fin(k)
    if plotinput:
        plt.figure()
        plt.title(title)
        plt.plot(x)
        plt.show()
    
    # Identify how many rows were provided
    rows, cols = filter.shape
    # Operate with each individual filter set
    x_tmp = np.copy( x )
    nsteps = NN - 4
    for row_n in range(rows):
        row = filter[row_n] # Capture individual row
        A1 = row[0]
        A2 = row[1]
        B0 = row[2]
        B1 = row[3]
        B2 = row[4]
        T = 3
        for _ in range(nsteps):
            T = T + 1
            # Apply Filtering Specified by Individual Row
            y[T] = A1*y[T-1] + A2*y[T-2] + B0*x_tmp[T] + B1*x_tmp[T-1] +  B2*x_tmp[T-2]
        # Copy New output into temporary input
        x_tmp = np.copy( y )
    # Copy finalized output into *ytime* for plotting
    ytime = np.copy( x_tmp )
    # Plot Filtered Output
    plt.figure(1)
    plt.plot(x,'k--')
    plt.plot(ytime,'k')
    plt.title(title)
    plt.grid(which='both')
    plt.show()
    