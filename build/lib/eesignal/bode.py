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
#   - Transfer Function Bode Plot Generator         bode
#   - S-Domain Bode Plot Generator                  sbode
#   - Z-Domain Bode Plot Generator                  zbode
#################################################################################

import matplotlib.pyplot as plt
import numpy as np
from math import pi, exp, cos, sin, log, sqrt

# Define System Bode Plotting Function
def bode(system,mn=-2,mx=3,npts=100,gtitle="",xlim=False,ylim=False,sv=False):
	""" System Bode Plotting Function
	
	A simple function to generate the Bode Plot for magnitude
	and frequency given a transfer function system.
	
	Required Arguments
	------------------
	system:		The Transfer Function; can be provided as the following:
				- 1 (instance of lti)
				- 2 (num, den)
				- 3 (zeros, poles, gain)
				- 4 (A, B, C, D)
				
	Optional Arguments
	------------------
	mn:			The minimum frequency (as an exponent to 10, e.g. 10^mn)
				to be calculated for. Default is -2.
	mx:			The maximum frequency (as an exponent to 10, e.g. 10^mx)
				to be calculated for. Default is 3.
	npts:		The number of points over which to calculate the system.
				Default is 100.
	gtitle:		Additional string to be added to plot titles;
				default is "".
	xlim:		Limit in x-axis for graph plot. Accepts tuple of: (xmin, xmax).
				Default is False.
	ylim:		Limit in y-axis for graph plot. Accepts tuple of: (xmin, xmax).
				Default is False.
	sv:			Save the plots as PNG files. Default is False.
	
	Returns
	-------
	NONE:	Generates plot of magnitude and phase, does not return
			any numerical values.
	"""
	# Condition system input to ensure proper execution
	system = sys_condition(system,False)
	
	# Generate the frequency range to calculate over
	wover = np.logspace(mn,mx,npts)
	
	# Calculate the bode system
	w, mag, ang = sig.bode(system, wover)
	
	# Plot Magnitude
	magTitle = "Magnitude "+gtitle
	plt.title(magTitle)
	plt.plot(w, mag)
	plt.xscale("log")
	plt.grid(which="both")
	plt.ylabel("Magnitude (dB)")
	plt.xlabel("Frequency (rad/sec)")
	if xlim!=False:
		plt.xlim(xlim)
	if ylim!=False:
		plt.ylim(ylim)
	if sv:
		plt.savefig(magTitle+".png")
	plt.show()

	# Plot Angle
	angTitle = "Angle "+gtitle
	plt.title(angTitle)
	plt.plot(w, ang)
	plt.xscale("log")
	plt.grid(which="both")
	plt.ylabel("Angle (degrees)")
	plt.xlabel("Frequency (rad/sec)")
	if xlim!=False:
		plt.xlim(xlim)
	if ylim!=False:
		plt.ylim(ylim)
	if sv:
		plt.savefig(angTitle+".png")
	plt.show()

def sbode(f,NN=1000,title=""):
	"""
	SBODE Function
	
	Required Arguments:
	-------------------
	f:	The Input Function, must be callable
	
	Optional Arguments:
	-------------------
	NN:	The Interval over which to be generated, default=1000
	title:	The Title to be used for all plots, default=""
	
	Returns:
	--------
	NONE - Plots Generated
	"""
	W = np.linspace(0,1000,NN)
	H = np.zeros(NN, dtype = np.complex)

	for n in range(0,NN):
		s = 1j*W[n]
		H[n] = f(s)
	
	plt.figure(1)
	plt.semilogx(W,20*np.log10(abs(H)),'k')
	plt.ylabel('|H| dB')
	plt.title(title+" Magnitude")
	plt.grid(which='both')

	aaa = np.angle(H)
	for n in range(NN):
		if aaa[n] > pi:
			aaa[n] = aaa[n] - 2*pi

	plt.figure(2)
	plt.title(title+" Phase")
	plt.semilogx(W,(180/pi)*aaa,'k')
	plt.ylabel('H phase (degrees)')
	plt.xlabel('$\omega$ (rad/s)')
	plt.grid(which='both')
	plt.show()


def zbode(f,dt=0.01,NN=1000,title=""):
	"""
	FBODE Function
	
	Required Arguments:
	-------------------
	f:	The Input Function, must be callable
	
	Optional Arguments:
	-------------------
	dt:	The time-step used, default=0.01
	NN:	The Interval over which to be generated, default=1000
	title:	The Title to be used for all plots, default=""
	
	Returns:
	--------
	NONE - Plots Generated
	"""
	phi = np.linspace(0,2*pi,NN)

	z = np.zeros(NN, dtype = np.complex)
	H = np.zeros(NN, dtype = np.complex)
	for n in range(0,NN):
		z[n] = exp(1j*phi[n])
		H[n] = dt*f(z[n])
			
	plt.figure(1)
	plt.semilogx((180/pi)*phi,20*np.log10(abs(H)),'k')
	plt.ylabel('|H| dB')
	#plt.text(6,-15,'$\phi$ = {}'.format(round(.57,3)),fontsize=12)
	plt.title(title+" Magnitude")
	plt.grid(which='both')
	plt.show()

	aaa = np.angle(H)
	for n in range(NN):
		if aaa[n] > pi:
			aaa[n] = aaa[n] - 2*pi

	plt.figure(2)
	plt.semilogx((180/pi)*phi,(180/pi)*aaa,'k')
	plt.ylabel('H (degrees)')
	#plt.text(6,-15,'$\phi$ = {}'.format(round(.57,3)),fontsize=12)
	plt.grid(which='both')
	plt.xlabel('$\phi$ (degrees)')
	plt.title(title+" Phase")
	plt.show()
