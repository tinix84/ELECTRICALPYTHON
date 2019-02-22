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
#   - System Response Plotter:        sysresponse
#################################################################################

import numpy as np
import matplotlib.pyplot as plt

def zfiltersim( fin, filter, NN=1000, title="",plotinput=True,legend=True):
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
    legend:     An argument to control whether the legend is shown,
                default=True.
    
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
    plt.plot(x,'k--',label="Input")
    plt.plot(ytime,'k',label="Output")
    plt.title(title)
    plt.grid(which='both')
    if legend: plt.legend()
    plt.show()
    
# Define System Response Plotter function
def sys_response(system,npts=1000,dt=0.01,combine=True,gtitle="",xlim=False,
				stepResponse=True,rampResponse=False,parabolicResponse=False,sv=False):
	""" System Response Plotter Function
	
	Given a transfer function, plots the response against step, ramp, and
	parabolic inputs and plots the error for each of these functions.
	
	Required Arguments
	------------------
	system:		The Transfer Function; can be provided as the following:
				- 1 (instance of lti)
				- 2 (num, den)
				- 3 (zeros, poles, gain)
				- 4 (A, B, C, D)
	
	Optional Arguments
	------------------
	npts:				Number of steps to calculate over; default is 1000.
	dt:					Difference between each data point, default is 0.01.
	combine:			If combination of numerator and denominator is needed.
						This value should be set to "True" if the parts should be
						combined to show the complete system with feedback.
						Default is True.
	gtitle:				Additional string to be added to plot titles;
						default is "".
	stepResponse:		Plot the step-response and corresponding error;
						default is True.
	rampResponse:		Plot the ramp-response and corresponding error;
						default is False.
	parabolicResponse:	Plot the parabolic-response and corresponding error;
						default is False.
	xlim:				Limit in x-axis for graph plot. Accepts tuple of: (xmin, xmax).
						Default is False.
	sv:					Save the figures plotted. Default is False.
	
	Returns
	-------
	NONE:		Generates the plot of the desired responses,
				does not return numerical values.
	"""
	# Define Time Axis
	TT = np.arange(0,npts*dt,dt)
	
	# Condition system input to ensure proper execution
	system = sys_condition(system,combine)	
	
	# Allocate space for all outputs
	step = np.zeros(npts)
	ramp = np.zeros(npts)
	parabola = np.zeros(npts)
	errS = np.zeros(npts)
	errR = np.zeros(npts)
	errP = np.zeros(npts)
	
	# Generate Inputs
	for i in range(npts):
		step[i] = 1.0
		ramp[i] = (dt*i)
		parabola[i] = (dt*i)**(2)
	
	# Simulate Response for each input (step, ramp, parabola)
	# All 'x' values are variables that are considered don't-care
	x, y1, x = sig.lsim((system),step,TT)
	x, y2, x = sig.lsim((system),ramp,TT)
	x, y3, x = sig.lsim((system),parabola,TT)
	
	# Calculate error over all points
	for k in range(npts):
		errS[k] = step[k] - y1[k]
		errR[k] = ramp[k] - y2[k]
		errP[k] = parabola[k] - y3[k]
	
	# Plot responses if allowed
	if (stepResponse):
		plt.figure()
		plt.subplot(121)
		plt.title("Step Response "+gtitle)
		plt.plot(TT,y1,'k--', label="Step Response")
		plt.plot(TT,step,'k', label="Step Function")
		plt.grid()
		plt.legend()
		plt.xlabel("Time (seconds)")
		if xlim != False:
			plt.xlim(xlim)
		plt.subplot(122)
		plt.title("Step Response Error "+gtitle)
		plt.plot(TT,errS,'k', label="Error")
		plt.grid()
		plt.legend()
		plt.xlabel("Time (seconds)")
		if xlim != False:
			plt.xlim(xlim)
		plt.subplots_adjust(wspace=0.3)
		if sv:
			plt.savefig("Step Response ("+gtitle+").png")
		plt.show()
	if (rampResponse):
		plt.figure()
		plt.subplot(121)
		plt.title("Ramp Response "+gtitle)
		plt.plot(TT,y2,'k--', label="Ramp Response")
		plt.plot(TT,ramp,'k', label="Ramp Function")
		plt.grid()
		plt.legend()
		plt.xlabel("Time (seconds)")
		if xlim != False:
			plt.xlim(xlim)
		plt.subplot(122)
		plt.title("Ramp Response Error "+gtitle)
		plt.plot(TT,errR,'k', label="Error")
		plt.grid()
		plt.legend()
		plt.xlabel("Time (seconds)")
		if xlim != False:
			plt.xlim(xlim)
		plt.subplots_adjust(wspace=0.3)
		if sv:
			plt.savefig("Ramp Response ("+gtitle+").png")
		plt.show()
	if (parabolicResponse):
		plt.figure()
		plt.subplot(121)
		plt.title("Parabolic Response "+gtitle)
		plt.plot(TT,y3,'k--', label="Parabolic Response")
		plt.plot(TT,parabola,'k', label="Parabolic Function")
		plt.grid()
		plt.legend()
		plt.xlabel("Time (seconds)")
		if xlim != False:
			plt.xlim(xlim)
		plt.subplot(122)
		plt.title("Parabolic Response Error "+gtitle)
		plt.plot(TT,errP,'k', label="Error")
		plt.grid()
		plt.legend()
		plt.xlabel("Time (seconds)")
		if xlim != False:
			plt.xlim(xlim)
		plt.subplots_adjust(wspace=0.3)
		if sv:
			plt.savefig("Parabolic Response ("+gtitle+").png")
		plt.show()
    
# End of FILTERSIM.PY