# Originally Written by Dr. Dennis Sullivan
# Adapted by Joe Stanley

import matplotlib.pyplot as plt
import numpy as np
from math import pi, exp, cos, sin, log, sqrt

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
