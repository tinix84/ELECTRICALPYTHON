#################################################################################
#   EESIGNAL.PY
#
#   This file contains a variety of functions and constants related to signals.
#   These items will commonly  be used in Electrical Engineering Applications.
#
#   September 3, 2018
#   August 31, 2018
#
#   Written by Joe Stanley
#   Special Thanks To and Code Support From:
#   Steven Weeks
#   Dr. Dennis Sullivan
#
#   Special thanks to stackOverflow user: gg349 whose work is well documented
#   and used for the FFT calculations used throughout this library.
#
#   Included Functions
#   - FFT Coefficient Calculator:       fft_coef
#   - FFT Plotting Function:            fft_plot
#   - RMS Calculator:                   rms
#   - State Space Simulator:            statespace
#   - Step Function:                    step
#   - Phase Margin:                     pm
#   - Gain Margin:                      gm
#   - System Response Plotter:          sys_response
#   - Multi-Argument Convolution:       convolve
#   - Phase Lead System:                phase_lead
#   - Butterworth Min Order Solver:     but_minord
#   - Butterworth Filter Plotter:       filter_plt
#   - Butterworth Filter Generator:     butter_gen
#   - Chebyshev Filter I Term Solver:   cheb_I_terms
#   - Chebyshev Filter II Term Solver:  cheb_II_terms
#   - Filter Conversion:                filter_convert
#   - Filter Polynomial Factoring:      filter_factor
#   - Convolution Bar Graph Visualizer: convbar
#
#   Private Functions ( Those not Intended for Use Outside of Library )
#   - TF System Conditioning:           sys_condition
#   - Tupple to Matrix Converter:       tuple_to_matrix
#   - Numpy Array to Matrix Converter:  nparr_to_matrix
#
#   Private Classes ( Those not Intended for Use Outside of Library )
#   - Function Concatinator:			c_func_concat
#
#   Constants
#   - NaN (Not A Number):				NAN
#   - Type: Numpy Matrix:				matrix
#   - Type: Tuple:						tuple
#   - Type: Numpy Array:				ndarr
#   - Type: Integer:					tint
#   - Type: Float:						tfloat
#   - Type: Function Handle:			tfun
#
#   Submodules
#   - Bode Plot Generator               BODE.PY
#################################################################################
name = "eesignal"
ver = "1.2.1"

# Import Submodules
from . import bode
from . import filtersim

# Import necessary libraries
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad as integrate
from scipy import signal as sig
import sympy as sym
from sympy.abc import s as s

# Define constants
NAN = float('nan')
matrix = "<class 'numpy.matrixlib.defmatrix.matrix'>"
tuple = "<class 'tuple'>"
ndarr = "<class 'numpy.ndarray'>"
tint = "<class 'int'>"
tfloat = "<class 'float'>"
tfun = "<class 'function'>"
tnfloat = "<class 'numpy.float64'>"

# Define Convolution Bar-Graph Function:
def convbar(h, x, outline=True):
	"""
	CONVBAR Function:
	
	INPUTS:
	-------
	h: Impulse Response - Given as Array (Prefferably Numpy Array)
	x: Input Function - Given as Array (Prefferably Numpy Array)
	
	RETURNS:
	--------
	None.
	
	PLOTS:
	------
	Impulse Response: The bar-graph plotted version of h.
	Input Function:   The bar-graph plotted version of x.
	Convolved Output: The bar-graph plotted version of the convolution of h and x.
	"""
	
	# The impulse response
	M = len(h)
	t = np.arange(M)
	# Plot
	plt.subplot(121)
	if(outline): plt.plot(t,h,color='red')
	plt.bar(t,h,color='black')
	plt.xticks([0,5,9])
	plt.ylabel('h')
	plt.title('Impulse Response')
	plt.grid()

	# The input function
	N = len(x)
	s = np.arange(N)
	# Plot
	plt.subplot(122)
	if(outline): plt.plot(s,x,color='red')
	plt.bar(s,x,color='black')
	plt.xticks([0,10,19])
	plt.title('Input Function')
	plt.grid()
	plt.ylabel('x')

	# The output
	L = M+N-1
	w = np.arange(L)
	plt.figure(3)
	y = np.convolve(h,x)
	if(outline): plt.plot(w,y,color='red')
	plt.bar(w,y,color='black')
	plt.ylabel('y')
	plt.grid()
	plt.title('Convolved Output')
	plt.show()


# Define convolution function
def convolve(tuple):
	""" Multi-Argument Convolution Function
	
	Given a tuple of terms, convolves all terms in tuple to
	return one tuple as a numpy array.
	
	Arguments
	---------
	tuple:		Tuple of terms to be convolved.
				i.e. ( [1, 2], [3, 4], ..., [n-1, n] )
	
	Returns
	-------
	c:			The convolved set of the individual terms.
				i.e. np.array([ x1, x2, x3, ..., xn ])
	"""
	c = sig.convolve(tuple[0],tuple[1])
	if (len(tuple) > 2):
		# Iterate starting with second element and continuing
		for i in range(2,len(tuple)):
			c = sig.convolve(c,tuple[i])
	return(c)

# Define Butterworth Minimum Order Solver:
def but_minord(mxDev, w, wc=1):
	""" Minimum Order Solving Function
	
	Finds the minimum order allowable to create the butterworth
	filter to match the specified criteria.
	
	Arguments:
	----------
	mxDev:		The maximum allowable deviation.
	w:			The frequency at which mxDev occurs.
	wc:			The cuttoff frequency (default is 1).
	
	Returns:
	--------
	order:		An integer representing the minimum order."""
	
	# Scale to cuttoff frequency
	w = w/wc
	
	# Find the H2 value
	H2 = mxDev**2
	
	# Find Order
	order = int(round( np.log10( (1/H2) - 1 )/(2*np.log10(w)) ))
	
	# Return output
	return( order )
	
# Define Chebyshev Pole Calculator
def cheb_poles(n, a, b,type=1):
	""" Chebyshev Poles Calculation Function
	
	Purpose: Calculate and return a polynomial set (numpy array)
	describing the poles of a Chebyshev Filter.
	
	Required Arguments:
	-------------------
	n:		Filter Order
	a:		Calculated descriptive term "a"
	b:		Calculated descriptive term "b"
	
	Optional Arguments:
	-------------------
	type:	The Filter type, either 1 or 2 denoting
			Chebyshev type I or type II.
			
	Returns:
	--------
	totPole:	The complete system of poles returned as a
				numpy array representing polynomial coefficients.
				[s^m, ... s^2, s^1, s^0]
	"""
	totPole = np.array([1])
	ang = 180
	# Determine if order is odd or even
	if (n%2 == 0): # Even Order
		for i in range( int(n/2) ):
			div = 180/n
			if i == 0:
				ang -= div/2
			else:
				ang -= div
			ang_r = np.radians(ang)
			s = a*np.cos(ang_r)+1j*b*np.sin(ang_r)
			if type == 1: # Type I Cheb Filter
				pole = np.polymul( [1, -s], [1, -np.conj(s)] )
			elif type == 2: # Type 2 Cheb Filter
				pole = np.polymul( 1/np.array([1, -s]), 1/np.array([1, -np.conj(s)]) )
			else:
				print("WARNING: Cheb Filter must be of type 1 or 2.")
			totPole = np.polymul(totPole, pole)
	else: # Odd Order
		for i in range( int((n-1)/2) ):
			div = 180/n
			if i == 0:
				totPole = np.array([1, 1])
				if type == 2:
					totPole = 1/totPole
			ang -= div
			ang_r = np.radians(ang)
			s = a*np.cos(ang_r)+1j*b*np.sin(ang_r)
			if type == 1: # Type I Cheb Filter
				pole = np.polymul( [1, -s], [1, -np.conj(s)] )
			elif type == 2: # Type 2 Cheb Filter
				pole = np.polymul( 1/np.array([1, -s]), 1/np.array([1, -np.conj(s)]) )
			else:
				print("WARNING: Cheb Filter must be of type 1 or 2.")
			totPole = np.polymul(totPole, pole)
	return(totPole)
	
# Define Chebyshev Zero Calculator
def cheb_zeros(n):
	""" Chebyshev Zeros Calculation Function
	
	Purpose: Calculate and return a polynomial set (numpy array)
	describing the zeros of a Chebyshev Filter.
	
	Required Arguments:
	-------------------
	n:		Filter Order
			
	Returns:
	--------
	wk:		Each omega returned as a list
	zeros:	The complete system of zeros returned as a
			numpy array representing polynomial coefficients.
			[s^m, ... s^2, s^1, s^0] """
			
	zeros = np.array([1])
	wk = np.array([])
	for i in range(n):
		k = 2*i +1
		w = 1/np.cos( k*np.pi/(2*n) )
		if abs(w) < 1e10: # Test for value less than infinity
			zeros = np.polymul( zeros, [1, w] )
		wk = np.append(wk, w)
	for i in range(len(zeros)):
		if abs(zeros.item(i))<1e-10:
			zeros[i] = 0
		if zeros.item(i) < 0:
			zeros[i] = abs(zeros.item(i))
	return(wk, zeros)

# Define Chebyshev I Filter Terms Solver
def cheb_I_terms(ws, Hs, Hp, n=False):
	""" Chebyshev Type I Term Function
	
	Purpose: A function to calculate specific terms used
	in the process of design and development of a Chebyshev
	type I filter.
	
	Required Arguments:
	-------------------
	ws:		Stop-Band Frequency
	Hs:		Stop-Band Magnitude
	Hp:		Pass-Band Magnitude
	
	Optional Arguments:
	-------------------
	n:		Filter Order, used to force and bypass calculation of n.
	
	Returns:
	ep:		Epsilon of system
	n:		System Order
	alpha:	Alpha of system
	a:		A constant of system
	b:		B constant of system
	"""
	# Calculate Epsilon
	ep = np.sqrt( 1/Hp**2 - 1 )
	
	# Determine minimum order n if not provided
	if ( n == False ):
		n = (np.arccosh( (1/ep) * np.sqrt( 1/(Hs**2) - 1 ) ) * 
			( 1 / np.arccosh( ws ) ))

	# Cast as int and round as necessary
	nn = int( round(n, 0) )
	if nn < n: # fractional part exists
		n = nn + 1
	else:
		n = nn
	
	# Calculate alpha, a, b
	alpha = 1/ep + np.sqrt( 1 + 1/ep**2 )
	a = 1/2 * ( alpha**(1/n) - alpha**(-1/n) )
	b = 1/2 * ( alpha**(1/n) + alpha**(-1/n) )
	
	# Return Epsilon, n, alpha, a, b as tuple
	return(ep,n,alpha,a,b)

# Define Chebyshev II Filter Terms Solver
def cheb_II_terms(wp, Hs, Hp, n=False):
	""" Chebyshev Type II Term Function
	
	Purpose: A function to calculate specific terms used
	in the process of design and development of a Chebyshev
	type I filter.
	
	Required Arguments:
	-------------------
	wp:		Pass-Band Frequency
	Hs:		Stop-Band Magnitude
	Hp:		Pass-Band Magnitude
	
	Optional Arguments:
	-------------------
	n:		Filter Order, used to force and bypass calculation of n.
	
	Returns:
	ep:		Epsilon of system
	n:		System Order
	alpha:	Alpha of system
	a:		A constant of system
	b:		B constant of system
	"""
	# Calculate Epsilon
	ep = np.sqrt( Hs**2 / (1-Hs**2) )
	
	# Determine minimum order n if not provided
	if ( n == False ):
		n = (np.arccosh( np.sqrt( (1/ep**2) * (1/(1-Hp**2)) ) ) * 
			( 1 / np.arccosh( 1/wp ) ))
	
	# Cast as int and round as necessary
	nn = int( round(n, 0) )
	if nn < n: # fractional part exists
		n = nn + 1
	else:
		n = nn
	
	# Calculate alpha, a, b
	alpha = 1/ep + np.sqrt( 1 + 1/ep**2 )
	a = 1/2 * ( alpha**(1/n) - alpha**(-1/n) )
	b = 1/2 * ( alpha**(1/n) + alpha**(-1/n) )
	
	# Return Epsilon, n as tuple
	return(ep,n,alpha,a,b)
	
# Define System Factor Generator
def filter_factor(sys):
	""" Filter Factorization Function
	
	Purpose: accept a system polynomial and factor it as needed
	to create a set of 1st and 2nd order polynomial factors.
	
	Arguments:
	----------
	sys:	The polynomial passed to the factorization function.
	
	Returns:
	--------
	poles:	The numpy array of polynomial factors.
	
	"""
	# Find roots to evaluate expression
	expr = np.roots(sys)
	
	# Generate Arbitrary Array with zeros to be removed later
	poles = np.array([[0,0,0],[0,0,0]])

	# Iteratively Compute Factors
	for i in range(len(expr)):
		# If value is complex, use value and next value
		if np.iscomplex(expr.item(i)):
			# Perform calculation only if even term in set
			if (i%2 == 0):
				mult = np.array(np.polymul([1,-expr.item(i)],[1,-expr.item(i+1)])).real
				poles = np.append(poles,[mult],axis=0)
		# If value is real (all else)
		else:
			# If value is 0 (or 0 for all intents and purposes)
			# then factor must be s, (or s^2, which isn't accounted for)
			if expr.item(i) == 0 or abs(expr.item(0)) < 1e-12:
				poles = np.append(poles,[[0,1,0]],axis=0)
			# All other values indicate an (s-value) factor
			else:
				poles = np.append(poles,[[0,1,-expr.item(i)]],axis=0)
	# Remove Leading zero factors
	poles = poles[2:]
	
	# Return resulting factor terms
	return(poles)

# Define System Conditioning Function
def sys_condition(system,feedback):
	if ( len(system) == 2 ):		# System found to be num and den
		num = system[0]
		den = system[1]
		# Convolve numerator or denominator as needed
		if (str(type(num)) == tuple):
			num = convolve(num)		# Convolve terms in numerator
		if (str(type(den)) == tuple):
			den = convolve(den)		# Convolve terms in denominator
		if feedback: # If asked to add the numerator to the denominator
			ld = len(den) # Length of denominator
			ln = len(num) # Length of numerator
			if(ld > ln):
				num = np.append(np.zeros(ld-ln),num) # Pad beginning with zeros
			if(ld < ln):
				den = np.append(np.zeros(ln-ld),den) # Pad beginning with zeros
			den = den + num # Add numerator and denominator
		for i in range( len( num ) ):
			if (num[i] != 0):
				num = num[i:]		# Slice zeros off the front of the numerator
				break 				# Break out of for loop
		for i in range( len( den ) ):
			if (den[i] != 0):
				den = den[i:]		# Slice zeros off the front of the denominator
				break 				# Break out of for loop
		system = (num,den)  # Repack system
	return(system) # Return the conditioned system

# Define Filter to band-pass function
def filter_convert( sys, convn, convd=1, debug=False, TFprint=False):
	""" Filter Conversion Function
	
	Purpose: This function is developed to perform the polynomial
	shift and conversion as dictated by the inputs. This function
	is to be used to shift transfer function systems or to convert
	the system as prescribed.
	
	Required Arguments:
	-------------------
	sys:		The tuple of the system, (numerator, denominator)
	convn:		The numerator of the conversion factor
	
	Optional Arguments:
	-------------------
	convd:		The denominator of the conversion factor, default=1
	debug:		Print debugging information, default=False
	TFprint:	Print the resulting transfer function, default=False
	
	Returns:
	--------
	num:		The newly updated numerator polynomial
	den:		The newly updated denominator polynomial
	
	"""
	# Condition Symbolic Conversion terms
	convn = sym.expand( convn * s**0 )
	convd = sym.expand( convd * s**0 )
	
	# Pull numerator and denominator terms
	num = sys[0]
	den = sys[1]
	
	# Convert to symbolic system
	den = sum(co*s**i for i, co in enumerate(reversed(den)))
	num = sum(co*s**i for i, co in enumerate(reversed(num)))

	# Substitute (s^2 + wc^2)/s for s
	den1 = den.subs(s,(convn/convd))
	den2 = sym.expand(den1)
	num1 = num.subs(s,(convn/convd))
	num2 = sym.expand(num1)

	# Find the maximum order of the system
	enable = sym.degree(sym.poly(convd,s))
	order = sym.degree(sym.poly(den,s))
	m = s**(order*enable)

	# Multiply by maximum order to make all exponents positive
	exp_den = sym.expand( den2 * m )
	exp_num = sym.expand( num2 * m )
	
	# Find the leading coefficient of the system
	LC = sym.LC(sym.poly(exp_den))
	
	# Multiply by 1 / leading coefficient to reduce LC to 1
	den3 = sym.expand( exp_den * 1/LC )
	num3 = sym.expand( exp_num * 1/LC )

	# Generate output as numpy array
	final_den = sym.poly(den3)
	final_num = sym.poly(num3)
	den = np.asarray(final_den.all_coeffs()).astype(np.double)
	num = np.asarray(final_num.all_coeffs()).astype(np.double)

	# Print debug information if needed
	if debug:
		print(num1,num2)
		print(den1,den2)
		print(m)
		print(LC)
		print(final_num)
		print(final_den)
		print(den)
		print(num)
	if TFprint:
		print("\nTransfer Function:\n")
		# Calculate number of spaces needed
		nchar_den = len(str(exp_den))
		nchar_num = len(str(exp_num))
		leftover = nchar_den - nchar_num
		# Generate Numerator
		numstr = ""
		for i in range(int(leftover/2)):
			numstr += " "
		print(numstr + str(exp_num))
		# Generate fractional division
		frac = ""
		for i in range(int(nchar_den)):
			frac += "-"
		print(frac)
		# Print Denominator
		print(str(exp_den))
	
	# Return
	return(num, den)
	
# Define Analog Filter Plotting Function
def filter_plt(system,mn=-1,mx=3,npts=1000,yticks=False,forceticks=False,gtitle="",
				xlim=False,ylim=False,ysize=10,xticks=False,xsize=False,sv=False):
	""" Filter Plotting Function:
	
	Purpose: Generates a magnitude only bode plot specifically for filter design.
	
	Required Arguments:
	-------------------
	system:		The tupled system of (numerator, denominator)
	
	Optional Arguments:
	-------------------
	mn:			The minimum value (10^mn) to be calculated; default=-1.
	mx:			The maximum value (10^mx) to be calculated; default=3.
	npts:		The number of points to calculate over; default=1000.
	yticks:		An array of the points to plot as tickmarks; default=False.
	forceticks:	A value that will allow forced elimination of preexisting
				tickmarks and only leave those provided as yticks or xticks;
				default=False.
	gtitle:		A string to be presented as the plot title; default=""
	xlim:		The list representing the minimum and maximum values over
				which to plot in x-axis; default=False
	ylim:		Same as xlim, but for y-axis.
	ysize:		The font size for y-tickmarks; default is 10
	xticks:		Same as yticks, but for x-axis.
	xsize:		Same as ysize, but for x-axis.
	sv:			A value that will allow the saving of the plotted figure as
				a PNG image with filename: [sv] [gtitle].png, default=False.
	
	Returns:
	--------
	N/A
	
	"""
	# Generate omega to be plotted over
	w = np.logspace(mn,mx,npts)
	
	# Condition system input to ensure proper execution
	system = sys_condition(system,False)
	
	# Generate transfer function
	x, H = sig.freqs(system[0],system[1],w) # x is a don't care
	H = np.abs(H)

	# Plot
	plt.plot(w,H) # Generate Plot of |H| over w
	plt.title(gtitle)
	plt.xscale("log") # Plot over log-scale on x-axis
	plt.grid(which="both") # Display grid in both axes
	plt.xlabel("w (rad/sec)")
	plt.ylabel("| H(s) |")
	if(yticks!=False): # If a set of "yticks" are given, apply them
		ytix = plt.yticks()[0] # Gather original "yticks" as the plot generated automatically
		ytix = np.append(ytix,yticks) # Append new "yticks" to existing
		if(forceticks!=False): # If only supplied "yticks" are desired...
			ytix=yticks # ...only display those ticks.
		plt.yticks(ytix,fontsize=ysize) # Set "yticks" as tickmarks on plot and set fontsize
	if(xticks!=False): # If a set of "xticks" are given, apply them
		xtix = plt.xticks()[0] # Gather original "xticks" as the plot generated automatically
		xtix = np.append(xtix,xticks) # Append new "xticks" to existing
		if(forceticks!=False): # If only supplied "xticks" are desired...
			xtix=xticks # ...only display those ticks.
		plt.xticks(xtix,fontsize=xsize) # Set "xticks" as tickmarks on plot and set fontsize
	if(xlim!=False): # If a limit is provided, apply it
		plt.xlim(xlim) 
	if(ylim!=False):
		plt.ylim(ylim) # Set y-limit on plot
	if sv!=False:
		plt.savefig(sv+" "+gtitle+".png")
	plt.show() # Display Plot

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

# Define Gain Margin Calculator Function
def gm(tf,mn=-2,mx=3,npts=1000,err=1e-12,printout=False,ret=True,find=-180):
	""" Gain Margin Calculator
	
	Given a transfer function, calculates the gain margin (gm) and the
	frequency at which the gain margin occurs (wg).
	
	Required Arguments:
	-------------------
	tf:			The Transfer Function; can be provided as the following:
				- 1 (instance of lti)
				- 2 (num, den)
				- 3 (zeros, poles, gain)
				- 4 (A, B, C, D)
	
	Optional Arguments:
	-------------------
	mn:			The minimum frequency (as an exponent to 10, e.g. 10^mn)
				to be calculated for. Default is -2.
	mx:			The maximum frequency (as an exponent to 10, e.g. 10^mx)
				to be calculated for. Default is 3.
	npts: 		The number of points over which to calculate the system.
				Default is 100.
	err:		The maximum allowable error for an aproximation of zero
				(i.e. the difference between found value and zero).
				Default is 1e-12.
	printout:	If set to true, will automatically print the values.
				Default is False.
	ret:		If set to true, will return the gain margin and frequency.
				Default is True.
	find:		The value to be searched for. For gain margin, the default
				is -180 (degrees).
	
	Returns:
	--------
	wg:			First returned argument; the frequency at which gain margin
				occurs (in radians per second e.g. rad/s).
	gm:			Second and final returned argument; Gain Margin (in dB).
	
	"""
	# Initialize while loop control
	valid = True
	
	# Condition system input to ensure proper execution
	tf = sys_condition(tf,False)

	# Initialize values given numerator and denominator
	wover = np.logspace(mn,mx,npts)
	w, mag, ang = sig.bode((tf),wover)

	while(valid):
		# Re-Initialize variables
		isP = False
		isN = False
		gm = NAN
		wg = NAN
		
		# Perform iterative loop recursively set
		# smaller and smaller boundaries
		for i in range(len(mag)-1):
			if (ang[i] >= find):
				isP = True             			# Positive value found
				Pi = w[i]              			# Store frequency for positive value
				if abs(ang[i]-find) < err:		# if less than error, stop
					valid = False      			# stop while loop
					wg = Pi            			# store Phase Margin angle
					gm = mag[i]        			# store Phase Margin
					break              			# break out of for loop
			elif (ang[i] < find):
				isN = True            			# Negative value found
				Ni = w[i]              			# Store frequency for negative value
				if abs(ang[i]-find) < err:		# if less than error, stop
					valid = False      			# stop while loop
					wg = Ni            			# store gain Margin angle
					gm = mag[i]        			# store gain Margin
					break              			# break out of for loop

			if (isP and isN): # If both positive and negative values found
				if(Pi < Ni):          # Positive comes before negative
					x1 = np.log10(Pi) # Convert into value for logspace
					x2 = np.log10(Ni) # Convert into value for logspace
				elif (Ni < Pi):       # Negative comes before positive
					x1 = np.log10(Ni) # Convert into value for logspace
					x2 = np.log10(Pi) # Convert into value for logspace
				valid = True                            # Reset valid value
				wzoom = np.logspace(x1,x2,npts)       # Generate zoomed logspace
				w, mag, ang = sig.bode((tf),wzoom) # Generate Bode values
				break                                   # Break out of for loop
			else: # Not both positive and negative values were found
				valid = False # stop while loop
	
	if printout:
		print("Gain Margin:",gm,"dB at",wg,"rad/sec")
	
	if ret:
		return(wg, gm) # Return both the gain margin frequency (where it occurs) and the gain margin.

# Define Phase Margin Calculator Function
def pm(tf,mn=-2,mx=3,npts=1000,err=1e-12,printout=False,ret=True,find=0):
	""" Phase Margin Calculator
	
	Given a transfer function, calculates the phase margin (pm) and the
	frequency at which the phase margin occurs (wp).
	
	Required Arguments:
	-------------------
	tf:			The Transfer Function; can be provided as the following:
				- 1 (instance of lti)
				- 2 (num, den)
				- 3 (zeros, poles, gain)
				- 4 (A, B, C, D)
	
	Optional Arguments:
	-------------------
	mn:			The minimum frequency (as an exponent to 10, e.g. 10^mn)
				to be calculated for. Default is -2.
	mx:			The maximum frequency (as an exponent to 10, e.g. 10^mx)
				to be calculated for. Default is 3.
	npts: 		The number of points over which to calculate the system.
				Default is 100.
	err:		The maximum allowable error for an aproximation of zero
				(i.e. the difference between found value and zero).
				Default is 1e-12.
	printout:	If set to true, will automatically print the values.
				Default is False.
	ret:		If set to true, will return the phase margin and frequency.
				Default is True.
	find:		The value to be searched for. Default for phase margin is
				0 (dB).
	
	Returns:
	--------
	wp:			First returned argument; the frequency at which phase margin
				occurs (in radians per second e.g. rad/s).
	pm:			Second and final returned argument; Phase Margin (in degrees).
	
	"""
	# Initialize while loop control
	valid = True
	
	# Condition system input to ensure proper execution
	tf = sys_condition(tf,False)

	# Initialize values given numerator and denominator
	wover = np.logspace(mn,mx,npts)
	w, mag, ang = sig.bode((tf),wover)

	while(valid):
		# Re-Initialize variables
		isP = False
		isN = False
		pm = NAN
		wp = NAN
		
		# Perform iterative loop recursively set
		# smaller and smaller boundaries
		for i in range(len(mag)-1):
			if (mag[i] >= find):
				isP = True             			# Positive value found
				Pi = w[i]              			# Store frequency for positive value
				if abs(mag[i]-find) < err:		# if less than error, stop
					valid = False      			# stop while loop
					wp = Pi            			# store Phase Margin angle
					pm = ang[i]        			# store Phase Margin
					break              			# break out of for loop
			elif (mag[i] < find):
				isN = True             			# Negative value found
				Ni = w[i]              			# Store frequency for negative value
				if abs(mag[i]-find) < err:		# if less than error, stop
					valid = False      			# stop while loop
					wp = Ni            			# store Phase Margin angle
					pm = ang[i]        			# store Phase Margin
					break              			# break out of for loop

			if (isP and isN): # If both positive and negative values found
				if(Pi < Ni):          # Positive comes before negative
					x1 = np.log10(Pi) # Convert into value for logspace
					x2 = np.log10(Ni) # Convert into value for logspace
				elif (Ni < Pi):       # Negative comes before positive
					x1 = np.log10(Ni) # Convert into value for logspace
					x2 = np.log10(Pi) # Convert into value for logspace
				valid = True                            # Reset valid value
				wzoom = np.logspace(x1,x2,npts)       # Generate zoomed logspace
				w, mag, ang = sig.bode((tf),wzoom) # Generate Bode values
				break                                   # Break out of for loop
			else: # Not both positive and negative values were found
				valid = False # stop while loop
				
	pm -= 180     # subtract 180 degrees to find true phase margin
	if (pm < 0):  # If phase margin is less than zero, shift by 360 degrees
		pm += 360
	
	if printout:
		print("Phase Margin:",pm,"degrees at",wp,"rad/sec")
	
	if ret:
		return(wp, pm) # Return both the phase margin frequency (where it occurs) and the phase margin.

# Define Phase Lead Compensation Calculator
def phase_lead(system,desired,tolerance=5,printout=False,ret=True,plot=False):
	""" Phase Lead Compensation Calculator
	
	Given a transfer-function system, and a desired phase-margin,
	calculate the wp and wz to be used in a feedback system of type I
	(defined above or in EESignal.py help) necessary to achieve the
	desired phase-margin.
	
	Required Arguments:
	-------------------
	system:				The Transfer Function; can be provided as the following:
						- 1 (instance of lti)
						- 2 (num, den)
						- 3 (zeros, poles, gain)
						- 4 (A, B, C, D)
	desired:			The goal Phase-Margin.
	
	Optional Arguments:
	-------------------
	tolerance:			The additional phase given to make bring
						the output closer to the desired result.
						Default is 5.
	printout:			Print out the all values in order of calculation.
						Default is False.
	ret:				Determines if there are values returned.
						Default is True.
	plot:				Determines if plots are generated, will generate
						two graphs, original and corrected. Default is False.
						
	Return:
	-------
	wp:					The Pole value of Phase Lead circuit { G(s) }
	wz:					The Zero value of Phase Lead circuit { G(s) }
	pm:					The phase margin of the new system.
	
	"""
	# Find original phase-margin
	wpm, phm = pm(system)

	# From original phase margin, find phi-m
	phim = desired-phm+tolerance

	# Calculate alpha:
	alpha = (1+np.sin(np.radians(phim)))/(1-np.sin(np.radians(phim)))

	# Calculate wm
	wmp = -10*np.log10(alpha)
	wm,x = pm(system,find=wmp)

	# Calculate wp and wz
	wp = np.sqrt(alpha)*wm
	wz = wm/np.sqrt(alpha)
	
	# Condition system
	system = sys_condition(system,False)
	
	# Add feedback control to system
	if(len(system) == 2):
		num = system[0]
		den = system[1]
	else:
		num = sig.TransferFunction(system).num
		den = sig.TransferFunction(system).den
	nwz = np.array([1,wz])
	nwp = np.array([1,wp])
	num = convolve((num, nwz)) * wp
	den = convolve((den, nwp)) * wz
	sys = (num, den)
	
	# Calculate new Phase Margin
	nwp, npm = pm(sys)

	if printout:
		print("Original Phase Margin:")
		print("Phase Margin:",phm,"degrees at",wpm,"rad/sec")
		if plot:
			bode(system)
		print("Phi-M:",phim)
		print("Alpha:",alpha)
		print("Magnitude where Wm appears:",wmp,"dB")
		print("Wm:",wm,"rad/sec")
		print("Wp:",wp)
		print("Wz:",wz)
		if plot:
			bode(sys)
		print("Phase Margin:",npm,"degrees at",nwp,"rad/sec")
	
	if ret:
		return(wp,wz,npm)

# Define Function Concatinator Class
class c_func_concat:
	def __init__(self,funcs): # Initialize class with tupple of functions
		self.nfuncs = len(funcs) # Determine how many functions are in tuple
		self.func_reg = {} # Create empty keyed list of function handles
		for key in range(self.nfuncs): # Iterate adding to key
			self.func_reg[key] = funcs[key] # Fill keyed list with functions

	def func_c(self,x): # Concatenated Function
		rets = np.array([]) # Create blank numpy array to store function outputs
		for i in range(self.nfuncs):
			y = self.func_reg[i](x) # Calculate each function at value x
			rets = np.append(rets, y) # Add value to return array
		rets = np.asmatrix(rets).T # Convert array to matrix, then transpose
		return(rets)

# Define Step function
def step(t):
	return( np.heaviside( t, 1) )

# Tuple to Matrix Converter
def tuple_to_matrix(x,yx):
	n = yx(x) # Evaluate function at specified point
	n = np.asmatrix(n) # Convert tuple output to matrix
	n = n.T # Transpose matrix
	return(n)

# Numpy Array to Matrix Converter
def nparr_to_matrix(x,yx):
	n = yx(x) # Evaluate function at specified point
	n = np.asmatrix(n) # Convert np.arr output to matrix
	if n.shape[1] != 1: # If there is more than 1 column
		n = np.matrix.reshape(n,(n.size,1)) # Reshape
	return(n)

def statespace(A,B,x=0,f=0,solution=2,C=False,D=False,npts=9999,NN=10000,dt=0.01,
		xlim=False,ylim=False,gtitle="",ret=False,plot=True,pltfn=False,sv=False):
	""" Plots the state-space simulation of an arbitrary set of matricies.

	Required Arguments:
	--------------------
	A :			Matrix A; if not type=numpy.matrix, converts to numpy.matrix
	B :			Matrix B; if not type=numpy.matrix, converts to numpy.matrix

	Optional Arguments:
	--------------------
	x :			Matrix x; if not type=numpy.matrix, converts to numpy.matrix
	f :			Forcing Function; must be provided as callable function that
				will return any/all forcing function Arguments needed as
				numpy matrix (preferred), numpy array, or tuple.
				Forcing function(s) can be provided as tuple of function
				handles, system will automatically concatenate their output
				to a matrix that can be handled.
	solution:	Determines What Type of Solution Simulation is Required;
				Default of solution is 2
				0=zero-input	( No Forcing Function )
				1=zero-state	( No Initial Conditions )
				2=total			( Both Initial Conditions and Forcing Function )
				3=total, output ( Both ICs and FFs, also plot combined output )
	npts: 		Changes the range of simulation; defualt=9999
	NN:			Number of descrete points; default=10,000
	dt:			Delta-t, step-size; default=0.01
	xlim:		Limit in x-axis for graph plot.
	ylim:		Limit in y-axis for graph plot.
	gtitle:		Additional String for Plot Title
	ret:		If true: returns state space terms
	plot:		If true: Plots individual state space terms
	pltfn:		If true: Plots original Forcing Functions

	Figures:
	--------
	Forcing Functions:		The plot of forcing functions, only provided if pltfn is true.
	State Variables:		The plot of state variables, always provided if plot is true.
	Combined Output:		The plot of the combined terms in the output, provided if C and D are not False.

	"""

	# Test for NN and npts
	if (npts >= NN):
		print("WARNING: NN must be greater than npts; NN="+str(NN)+"npts="+str(npts))
		print("Autocorrecting npts to be NN-1.")
		npts = NN-1

	# Test for C and D matricies
	mC = str(type(C))
	mD = str(type(D))
	if (((mC==matrix) or (mC==ndarr) or (mC==tuple)) and
		((mD==matrix) or (mD==ndarr) or (mD==tuple))):
		if (solution!=3):
			print("WARNING: C and D matricies provided, but solution requested "+
					"does not include combined output.")
	elif (((mC==matrix) or (mC==ndarr) or (mC==tuple)) and
		(mD!=matrix) and (mD!=ndarr) and (mD!=tuple)):
		if (D==False):
			print("WARNING: D matrix not provided; D now assumed to be 0.")
			D = np.matrix('0')
	else:
		C = np.matrix('0')
		D = np.matrix('0')
		if (solution==3):
			print("WARNING: Combined output requested, but no matricies C and D given.")
			print("         Solution being set to: 2 - Complete Simulation")
			solution = 2


	# Create values for input testing
	mA = str(type(A))
	mB = str(type(B))
	mx = str(type(x))
	mC = str(type(C))
	mD = str(type(D))
	if (str(type(f)) == tfun): # if f is a function, test as one
		mF = str(type(f(1))) # f should return: int, float, tuple, np.arr, np.matrix
	elif (str(type(f)) == tuple): # if f is tupple of arguments
		if (str(type(f[0])) == tfun): #if first argument is a function
			c_funcs = c_func_concat(f) # concatinate functions into one
			mF = "MultiFunctions" # label as multiple concatenated functions
		else:
			mF = "NA" # Can't handle function type
	else:
		mF = "NA" # Can't handle function type

	# Test for x input
	if (mx!=matrix) and (mx!=ndarr) and (mx!=tuple):
		if x==0: # No specified initial conditions
			if (mA==matrix): # Use A matrix as reference
				rA = A.shape[0]
				x = np.asmatrix(np.zeros(rA)).T
				mx = str(type(x))
				print("WARNING: No input x (Initial Condition) given.")
				if (solution!=1) and (solution!=3):
					solution = 1
					print("\n         Solution type changed to 1: Zero-State.")
			elif (mB==matrix): # Use B matrix as reference
				rB = B.shape[0]
				x = np.asmatrix(np.zeros(rB)).T
				mx = str(type(x))
				print("WARNING: No input x (Initial Condition) given.")
				if (solution!=1) and (solution!=3):
					solution = 1
					print("\n         Solution type changed to 1: Zero-State.")
			else:
				raise ValueError("ERROR: No x matrix (Initial Conditions) given,"+
								"\nNot enough additional information to infer x matrix.")

	# Test for matrix inputs
	if (mA!=matrix) or (mB!=matrix) or (mx!=matrix) or (mC!=matrix) or (mD!=matrix):
		# One or more arguments are not of type numpy.matrix
		# Convert to type matrix
		print("WARNING: Converting one or more input matricies to type: numpy matrix")
		A = np.asmatrix(A)
		B = np.asmatrix(B)
		x = np.asmatrix(x)

	# Gather dimensions of inputs
	rA, cA = A.shape
	rB, cB = B.shape
	rx, cx = x.shape
	rC, cC = C.shape
	rD, cD = D.shape
	rF, cF = 1, 1 # Defualt for a function returning one value

	if (mF==tuple): # If function returns tuple
		print("WARNING: Converting Forcing Function from type: tuple")
		fn = lambda x: tuple_to_matrix(x, f) # Use conversion function
		rF, cF = fn(1).shape # Prepare for further testing
	elif (mF==ndarr): # If function returns numpy array
		print("WARNING: Converting Forcing Function from type: numpy array")
		fn = lambda x: nparr_to_matrix(x, f) # Use conversion function
		rF, cF = fn(1).shape # Prepare for further testing
	elif (mF==tint) or (mF==tfloat) or (mF==tnfloat): # If function returns int or float or numpy float
		fn = f # Pass function handle
	elif (mF==matrix): # If function returns matrix
		fn = f # Pass function handle
		rF, cF = fn(1).shape # Prepare for further testing
	elif (mF=="MultiFunctions"): # There are multiple functions in one argument
		print("WARNING: Casting multiple forcing functions to output type: numpy matrix")
		fn = c_funcs.func_c # Gather function handle from function concatenation class
		rF, cF = fn(1).shape # Prepare for further testing
	elif (mF=="NA"): # Function doesn't meet requirements
		raise ValueError("ERROR: Forcing function does not meet requirements."+
						"\nFunction doesn't return data type: int, float, numpy.ndarray"+
						"\n or numpy.matrixlib.defmatrix.matrix. Nor does function "+
						"\ncontain tuple of function handles. Please review function.")

	# Test for size correlation between matricies
	if (cA != rA): # A isn't nxn matrix
		raise ValueError("ERROR: Matrix 'A' is not NxN matrix.")
	elif (rA != rB): # A and B matricies don't have same number of rows
		if (B.size % rA) == 0: # Elements in B divisible by rows in A
			print("WARNING: Reshaping 'B' matrix to match 'A' matrix.")
			B = np.matrix.reshape(B,(rA,int(B.size/rA))) # Reshape Matrix
		else:
			raise ValueError("ERROR: 'A' matrix dimensions don't match 'B' matrix dimensions.")
	elif (rA != rx): # A and x matricies don't have same number of rows
		if (x.size % rA) == 0: # Elements in x divisible by rows in A
			print("WARNING: Reshaping 'x' matrix to match 'A' matrix.")
			x = np.matrix.reshape(x,(rA,1)) # Reshape Matrix
		else:
			raise ValueError("ERROR: 'A' matrix dimensions don't match 'B' matrix dimensions.")
	elif (cB != rF) or (cF != 1): # Forcing Function matrix doesn't match B matrix
		raise ValueError("ERROR: 'B' matrix dimensions don't match forcing function dimensions.")
	elif (solution==3) and (cC != cA) or (rC != 1): # Number of elements in C don't meet requirements
		raise ValueError("ERROR: 'C' matrix dimensions don't match state-space variable dimensions.")
	elif (solution==3) and ((cD != rF) or (rD != 1)): # Number of elements in D don't meet requirements
		if (cD == rD) and (cD == 1) and (D[0] == 0): # D matrix is set to [0]
			D = np.asmatrix(np.zeros(rF)) # Re-create D to meet requirements
			print("WARNING: Autogenerating 'D' matrix of zeros to match forcing functions.")
		else:
			raise ValueError("ERROR: 'D' matrix dimensions don't match forcing function dimensions.")

	# Test for forcing function
	if (f==0) and (solution!=0):
		print("WARNING: No forcing function provided.\n         "+
				"Solution type changed to 0: Zero-Input")
		solution = 0 # Change to Zero-Input calculation

	# Start by defining Constants
	T = 0
	TT = np.arange(0,(dt*(NN)),dt)
	yout = 0

	# Define list of strings for plot output
	soltype = ["(Zero-Input)","(Zero-State)","(Complete Simulation)","(Complete Sim., Combined Output)"]

	# Create a keyed list of state-space variables
	xtim = {}
	xtim_len = rA # Number of Rows in A matrix
	for n in range(xtim_len):
		key = n #Each key should be the iterative variable
		xtim_init = np.zeros(NN) #Define the initial array
		xtim[key] = xtim_init #Create each xtim

	# Create a keyed list of function outputs
	if (mF!=tint) and (mF!=tfloat):
		fn_arr = {}
		for n in range(rF):
			key = n #Each key should be the iterative variable
			fn_init = np.zeros(NN) #Define the initial array
			fn_arr[key] = fn_init #Create each fn_arr
			fnc = rF
	else:
		fn_arr = np.zeros(NN) #Create the fn_arr
		fnc = 1

	# When asked to find zero-state, set all ICs to zero
	if solution == 1:
		for n in range(xtim_len):
			x[n] = 0 #Set each value to zero

	# Finite-Difference Simulation
	for i in range(0,npts):
		for n in range(xtim_len):
			xtim[n][i] = x[n] #xtim[state-variable][domain] = x[state-variable]
		# Create Forcing Function output

		if fnc > 1: # More than one forcing function
			for n in range(fnc):
				fn_arr[n][i] = np.asarray(fn(T))[n][0]
		else: # only one forcing function
			fn_arr[i] = fn(T)

		if solution == 0: #Zero-input, no added function input
			x = x + dt*A*x
		else: #Zero-state or Total, add function input
			x = x + dt*A*x + dt*B*fn(T)
			if solution==3:
				yout = yout + dt*D*fn(T)

		T = T+dt #Add discrete increment to T

	# Plot Forcing Functions
	if (pltfn):
		fffig = plt.figure("Forcing Functions")
		if fnc > 1:
			for x in range(fnc):
				plt.plot(TT,fn_arr[x],label="f"+str(x+1))
		else:
			plt.plot(TT,fn_arr,label="f1")
		if xlim!=False:
			plt.xlim(xlim)
		if ylim!=False:
			plt.ylim(ylim)
		plt.title("Forcing Functions "+gtitle)
		plt.xlabel("Time (seconds)")
		plt.legend(title="Forcing Functions")
		plt.grid()
		if sv:
			plt.savefig('Simulation Forcing Functions.png')
		if plot:
			plt.show()

	# Plot each state-variable over time
	stvfig = plt.figure("State Variables")
	for x in range(xtim_len):
		plt.plot(TT,xtim[x],label="x"+str(x+1))
	if xlim!=False:
			plt.xlim(xlim)
	if ylim!=False:
		plt.ylim(ylim)
	plt.title("Simulated Output Terms "+soltype[solution]+gtitle)
	plt.xlabel("Time (seconds)")
	plt.legend(title="State Variable")
	plt.grid()
	if sv:
		plt.savefig('Simulation Terms.png')
	if plot:
		plt.show()

	# Plot combined output
	if (solution==3):
		cofig = plt.figure("Combined Output")
		C = np.asarray(C) # convert back to array for operation
		for i in range(cC):
			yout = yout + xtim[i]*C[0][i] # Sum all st-space var mult. by their coeff
		yout = np.asarray(yout) # convert output to array for plotting purposes
		plt.plot(TT,yout[0])
		if xlim!=False:
			plt.xlim(xlim)
		if ylim!=False:
			plt.ylim(ylim)
		plt.title("Combined Output "+gtitle)
		plt.xlabel("Time (seconds)")
		plt.grid()
		if sv:
			plt.savefig('Simulation Combined Output.png')
		if plot:
			plt.show()

	# Return Variables if asked to
	if ret:
		return(TT, xtim)


# RMS Calculating Function
def rms(f, T):
	""" Calculates the RMS value of the provided function.

	Arguments
	----------
	f : the periodic function, a callable like f(t)
	T : the period of the function f, so that f(0)==f(T)

	Returns
	-------
	RMS : the RMS value of the function (f) over the interval ( 0, T )

	"""
	fn = lambda x: f(x)**2
	integral = integrate(fn,0,T)
	RMS = np.sqrt(1/T*integral)
	return(RMS)

# FFT Coefficient Calculator Function
def fft_coef(f, N, T=1, return_complex=False):
	"""Calculates the first 2*N+1 Fourier series coeff. of a periodic function.

	Given a periodic, function f(t) with period T, this function returns the
	coefficients a0, {a1,a2,...},{b1,b2,...} such that:

	f(t) ~= a0/2+ sum_{k=1}^{N} ( a_k*cos(2*pi*k*t/T) + b_k*sin(2*pi*k*t/T) )

	If return_complex is set to True, it returns instead the coefficients
	{c0,c1,c2,...}
	such that:

	f(t) ~= sum_{k=-N}^{N} c_k * exp(i*2*pi*k*t/T)

	where we define c_{-n} = complex_conjugate(c_{n})

	Refer to wikipedia for the relation between the real-valued and complex
	valued coeffs at http://en.wikipedia.org/wiki/Fourier_series.

	Arguments
	----------
	f : the periodic function, a callable like f(t)
	T : the period of the function f, so that f(0)==f(T)
	N_max : the function will return the first N_max + 1 Fourier coeff.

	Returns
	-------
	if return_complex == False, the function returns:

	a0 : float
	a,b : numpy float arrays describing respectively the cosine and sine coeff.

	if return_complex == True, the function returns:

	c : numpy 1-dimensional complex-valued array of size N+1

	"""
	# From Shanon theoreom we must use a sampling freq. larger than the maximum
	# frequency you want to catch in the signal.
	f_sample = 2 * N
	# we also need to use an integer sampling frequency, or the
	# points will not be equispaced between 0 and 1. We then add +2 to f_sample
	t, dt = np.linspace(0, T, f_sample + 2, endpoint=False, retstep=True)

	y = np.fft.rfft(f(t)) / t.size

	if return_complex:
	   return y
	else:
	   y *= 2
	   return y[0].real, y[1:-1].real, -y[1:-1].imag


# FFT Plotting Function
def fft_plot(f, N, T=1, mn=False, mx=False, fftplot=True, absolute=False, title=False, plotall=True):
	""" Plots the FFT of the provided function as a stem plot.

	Arguments
	----------
	f : the periodic function, a callable like f(t)
	T : the period of the function f, so that f(0)==f(T)
	N_max : the function will return the first N_max + 1 Fourier coeff.
	mn : the minimum time value of the original signal
	mx : the maximum time value of the original signal

	Returns
	-------
	if fftplot=True, the function returns:
	Plot of FFT output of function

	if mx != False, the function returns:
	Approximation of original signal from FFT results

	if absolute=True, the function will:
	Return absolute values of the coefficients
	
	if plotall=True, the function will:
	Plot each summed frequency

	"""

	# Calculate FFT and find coefficients
	a0, a, b = fft_coef(f, N, T)

	# If provided a title, add it to the title string
	tStr = ""
	if title!=False:
		tStr = title

	# Define Range values for plots
	rng = range(1,len(a)+1,1)
	xtic = range(0,len(a)+1,1)
	
	# Plot FFT results with respect to their sign
	if fftplot and not absolute:
		# Set up Arguments
		a0x = [0,0]
		a0y = [0,a0/2]

		# Plot
		plt.title("Fourier Coefficients"+tStr)
		plt.plot(a0x,a0y,'g',label="DC-Term")
		plt.stem(rng,a,'r','ro',label="A-Terms")
		plt.stem(rng,b,'b','bo',label="B-Terms")
		plt.legend()
		if(len(xtic) < 50):
			plt.xticks(xtic)
		plt.show()

	# Plot absolute value of FFT results
	if fftplot and absolute:
		# Set up Arguments
		a0x = [0,0]
		a0y = [0,abs(a0)/2]

		# Plot
		plt.title("Fourier Coefficients"+tStr)
		plt.plot(a0x,a0y,'g',label="DC-Term")
		plt.stem(rng,np.abs(a),'r','ro',label="A-Terms")
		plt.stem(rng,np.abs(b),'b','bo',label="B-Terms")
		plt.legend()
		if(len(xtic) < 50):
			plt.xticks(xtic)
		plt.show()

	# Plot original function as described by FFT results
	if mx!=False:
		# Create domain variable
		x = np.arange(mn,mx,(mx-mn)/1000)
		# Set output to DC constant
		yout = np.ones(len(x))*a0
		# Plot each iteration of the Fourier Series
		for k in range(1,N):
			if plotall:
				plt.plot(x,yout)
			yout = yout + a[k-1]*np.cos(k*2*np.pi*x/T) + b[k-1]*np.sin(k*2*np.pi*x/T)
		plt.plot(x,yout)
		plt.title("Fourier Series Summation"+tStr)
		plt.show()
