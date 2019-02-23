#################################################################################
#   FILTER.PY
#
#   This file contains a variety of functions and constants related to filtering.
#   These items will commonly be used in Electrical Engineering Applications.
#
#   February 22, 2019
#
#   Written by Joe Stanley
#   Special Thanks To and Code Support From:
#   Dr. Dennis Sullivan
#
#   Included Functions:
#   - Phase Lead System:                phase_lead
#   - Butterworth Min Order Solver:     but_minord
#   - Butterworth Filter Plotter:       plot
#   - Butterworth Filter Generator:     butter_gen
#   - Chebyshev Filter I Term Solver:   cheb_I_terms
#   - Chebyshev Filter II Term Solver:  cheb_II_terms
#   - Filter Conversion:                convert
#   - Filter Polynomial Factoring:      factor
#   - Phase Margin:                     pm
#   - Gain Margin:                      gm
#   - FIR Filter Design Assistant:      firdesign
#   - Automatic Butterworth Builder:    autobutter
#
#   Private Functions ( Those not Intended for Use Outside of Library )
#   - TF System Conditioning:           sys_condition
#################################################################################

# Import Required Libraries
import numpy as np
import scipy.signal as sig
import matplotlib.pyplot as plt

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
def factor(poly):
	""" Filter Factorization Function
	
	Purpose: accept a system polynomial and factor it as needed
	to create a set of 1st and 2nd order polynomial factors.
	
	Arguments:
	----------
	poly:   The polynomial passed to the factorization function.
	
	Returns:
	--------
	poles:	The numpy array of polynomial factors.
	
	"""
	# Find roots to evaluate expression
	expr = np.roots(poly)
	
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
def convert( sys, convn, convd=1, debug=False, TFprint=False):
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
def plot(system,mn=-1,mx=3,npts=1000,yticks=False,forceticks=False,gtitle="",
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

# Define FIR Design Function
def firdesign(f,MM,filtermode,Mstart=None,sigma=None,nc=None,dt=1e-5,
              NN=1000,window=2,dispfund=True,scalef=None,freqs=None,
              axis=None,plotfast=True,maxaproach=10,printscale=False):
    """
    FIRDESIGN Function
    
    Purpose:
    --------
    This function is designed to assist users in designing an FIR filter
    by accepting the basic parameters and plotting the specified outputs.
    
    Required Arguments:
    -------------------
    f:          The input system function
    MM:         Design parameter, number of points
    filtermode: Specifies Low-, High-, or Band-Pass as follows:
                    1 = Low-Pass
                    2 = High-Pass
                    3 = Band-Pass
    
    Optional Arguments:
    -------------------
    dt:         The time-step size; default=1e-5.
    NN:         The number of time-steps; default=1000
    window:     Specifies whether or not to use a square or Gaussian
                filter window as follows:
                    1 = Square
                    2 = Gaussian (default)
                    (function) = Uses the function to specify the window
    Mstart:     Parameter for Low-Pass Envelope; default=None
    sigma:      Parameter for Band-Pass Envelope; default=None
    nc:         Parameter for High-Pass Envelope; default=None
    dispfund:   Control argument to display the fundamental frequency;
                default=True
    scalef:     The scaling factor to set maximum of the output FFT;
                default=None
    freqs:      The set of frequencies to be used in the x-axis label;
                default=None
    axis:       The bounds of the x- and y-axes; default=None
    plotfast:   Control argument to allow the system to plot "in process";
                default=True
    maxaproach: Limit of how many recursive attempts to achieve appropriately
                scale the FFT limit; default=10
    printscale: Control value used to enable printinf of the final scalef value;
                default=False
    
    Returns:
    --------
    NONE (Plots Generated)
    """
    # Generate the required constants and arrays
    N2 = NN//2 # '//' operator is integer-division
    M = N2
    K = 1
    scale = 1
    x = np.zeros(NN)
    y = np.zeros(NN)
    H = np.zeros(NN)
    TT = np.linspace(0,dt*(NN-1),NN)
    DF = 1/(dt*NN)
    FF = np.linspace(0,DF*(NN-1),NN)
    
    # Calculate Fundamental Frequency
    if(dispfund):
        print("Fundamental Frequency:", DF,"Hz")
    
    # Append 0 to the frequencies so that it's plotted
    if(freqs!=None):
        freqs = np.append(np.array([0]),freqs)
    
    # Capture Input Function Over Range
    y[0] = 0
    for n in range(1,NN):
        x[n] = f(dt*n)
        
    # Generate FFT Decomposition
    X = (2/NN)*np.fft.fft(x)
    
    # Plot Original Function
    plt.subplot(321)
    plt.plot(TT,x,'k')
    plt.title('Input System')
    plt.ylabel('x[n]')
    plt.xlabel('T (sec)')
    plt.grid()
    # Plot FFT Decomposition of Function
    plt.subplot(322)
    plt.plot(FF,abs(X),'k')
    plt.grid()
    if(axis!=None):
        plt.axis(axis)
    if(freqs!=None):
        plt.xticks(axis)
    plt.xlabel('Freq (Hz)')
    plt.title('FFT Decomposition')
    if(plotfast):
        plt.show()
    
    # Generate Filter Window to be Convolved With
    if(filtermode==1): # Mode set for Low-Pass
        if(Mstart==None):
            raise ValueError("ERROR: Mstart Must be Specified.")
        if(window==1):
            for n in range(0,Mstart):  # Rectangular Window
                H[n] = 1
        elif(window==2):
            for n in range(Mstart,N2): # Gaussian Window   
                H[n] =np.exp(-0.5*( ((n-Mstart)/1.5)**2))
        else:
            for n in range(0,Mstart):  # Personalized Window
                H[n] = window( n ) # *window* is a function handle
    elif(filtermode==2): # Mode set for High-Pass
        if(nc==None):
            raise ValueError("ERROR: nc Must be Specified.")
        if(window==1):
            for n in range(nc,N2):  # Rectangular Window
                H[n] = 1
        elif(window==2):
            for n in range(0,N2): # Gaussian Window   
                H[n] =np.exp(-0.5*( ((n-nc)/4)**2))
        else:
            for n in range(0,Mstart):  # Personalized Window
                H[n] = window( n ) # *window* is a function handle
    elif(filtermode==3): # Mode set for Band-Pass
        if(window==1):
            raise ValueError("ERROR: Window Set Improperly.",
                             "Window may not be set to '1' (Square)",
                             "for Band-Pass filters (mode=3).")
        elif(window==2):
            if(nc==None):
                nc = 80//DF
            if(sigma==None):
                raise ValueError("ERROR: Sigma Must be Specified.")
            for n in range(0,N2): # Gaussian Window   
                H[n] =np.exp(-0.5*( ((n-nc)/sigma)**2))
        else:
            for n in range(0,N2):  # Personalized Window
                H[n] = window( n ) # *window* is a function handle
    else:
        raise ValueError("ERROR: Filter Mode Set Improperly.",
                         "Mode must be set between 1 and 3.")
    
    # Reflect around zero to get the negative frequencies.
    for n in range(1,M):
        H[NN-n] = H[n]
    
    # Multiply X and H to get Y, freq domain output    
    Y = np.copy( X )
    for n in range(0,NN):
        Y[n] = H[n]*X[n]
    
    # Use an inverse FFT to Display Filtered System
    y = (NN/2)*np.fft.ifft(Y).real
    
    # Plot Filtered System
    plt.subplot(323)
    plt.plot(y,'k')
    plt.ylabel('y[n]')
    plt.title('Filtered System')
    plt.grid()
    # Plot Filtered FFT
    plt.subplot(324)
    plt.plot(FF,H,'k--')
    plt.plot(FF,abs(Y),'k')
    if(axis!=None):
        plt.axis(axis)
    if(freqs!=None):
        plt.xticks(axis)
    plt.title('Proposed Filter in Freq. Domain')
    plt.grid()
    if(plotfast):
        plt.show()
    
    # Use an inverse FFT to capture the Filter
    h = (NN/2)*np.fft.ifft(H).real
    hmax = max(h)
    for n in range(NN):
        h[n] = h[n]/hmax
    
    #  Shift h  to hh
    hh = np.zeros(2*MM+1)
    for n in range(0,MM): # Shifted filter is MM*2 points long
        hh[n+MM] = h[n]
    hh[0] = 0
    for n in range(1,MM):    
        hh[n] = h[NN-MM+n]
    hsum = sum(abs(h))
    
    # Plot the Proposed Filter
    plt.subplot(325)
    plt.plot(h,'k')
    plt.grid()
    plt.title('Filter')
    # Plot Shifted Proposed Filter
    plt.subplot(326)
    plt.plot(hh,'ko')
    if(axis!=None):
        plt.axis(axis)
    if(freqs!=None):
        plt.xticks(axis)
    plt.title("Shifted Filter 'h'")
    plt.grid()
    if(plotfast):
        plt.show()
    
    # Determine if Specific "scalef-Factor" was given
    if(scalef!=None):
        K = scalef
        maxaproach = 1
    
    # Attempt to correct the "scalef-factor"
    for i in range(maxaproach):
        w = (K/hsum)*np.convolve(hh,x)
        z = np.zeros(NN)
        for n in range(NN):
            z[n] = w[n+MM]
        Z = (2/NN)*np.fft.fft(z)
        # Evaluate the Maximum Z
        mxZ = max(abs(Z))
        if(mxZ < 1): # Increase scalef-factor
            K += scale
        elif(mxZ > 1): # Overshoot!
            # Reset, then downsize scale
            K -= scale
            scale = scale/10
            K += scale
        else: # MaxZ equal to 1?
            break
    
    # Plot the Output
    plt.subplot(311)
    plt.plot(w,'k')
    if(axis!=None):
        plt.axis(axis)
    plt.grid()
    plt.subplot(311)
    plt.plot(z,'k')
    plt.grid()
    plt.subplot(313)
    plt.plot(FF,abs(Z),'k')
    plt.ylabel('|Z(w)|')
    plt.yticks([0,.1,.9,1.])
    if(axis!=None):
        plt.axis(axis)
    if(freqs!=None):
        plt.xticks(axis)
    plt.title('Final Filtered FFT of Output System')
    plt.grid()
    plt.show()
    
    # Print scalef Factor
    if(printscale):
        print("Scaling Value:", K)
        print("Max of FFT:", mxZ)

# Define Automatic Butterworth Filter Builder
def autobutter(wc,wp,wpm,ws,wsm,mode,W=None,n=None,mn=1e-1,mx=1e1,
               dt=0.01,respfreq=None):
    """
    AUTOBUTTER Function
    
    Purpose:
    --------
    This function is intended to automatically determine the transfer-
    function parameters required to meet given requirements of a filter.
    It is intended to design a Low-Pass, High-Pass or Band-Pass filter.
    
    Required Arguments:
    -------------------
    wc:         Center Frequency for LPF, HPF; Bandwidth for BP
    wp:         Passband Cutoff
    wpm:        Passband Magnitude at Cutoff
    ws:         Stopband Cutoff
    wsm:        Stopband Magnitude at Cutoff
    mode:       Desired Filter Control; Options are:
                    1: Low-Pass Filter
                    2: High-Pass Filter
                    3: Band-Pass Filter
    
    Optional Arguments:
    -------------------
    W:          Center Frequency of Band (only for BP); Required for BP;
                default=None
    n:          The number of poles; automatically determined if not
                specified (forced); default=None
    mn:         Frequency domain graph minimum x; default=1e-1
    mx:         Frequency domain graph maximum x; default=1e+1
    dt:         Sampling Time for Z Domain Filter, AKA: time step-size;
                default=0.01
    respfreq:   Response input Frequency; default=None
    
    Returns:
    --------
    ANY?
    """

    #mn = 1e-1     #Frequency domain graph minimum x
    #mx = 1e1      #Frequency domain graph maximum x
    w = np.linspace(mn, mx,1000)    #for Bode Plot; why not np.logspace?

    #respfreq = respfreq  #Response input Frequency; what was w_in for? Not used?

    #Poles calculation
    #Butterworth Filter Conditions
    wp = wp/wc # Passband cutoff
    ws = ws/wc  # Stopband cutoff

    p1 = np.log( (1/(wpm*wpm)) - 1)/(2*np.log(wp))
    p2 = np.log( (1/(wsm*wsm)) - 1)/(2*np.log(ws))
    temp_n = round(max(p1,p2))
    n = temp_n
    if (temp_n-max(p1,p2)<=0): #To always round up
        n=temp_n+1  
        
    #Butterworth Filter Calculation
    #Lowpass Filter Calculation
    num  = np.array([0, 0, wc**n ])

    #For Odd Pole Amounts
    if ((n%2) == 1):
        angle = np.pi/n
        den = np.array([ 1, wc])  #odd pole assigned
        int=1
        temp_angle = angle
        while(int < n):
            den = np.convolve(den,np.array([1, 2*wc*np.cos(temp_angle), wc*wc]))
            temp_angle = temp_angle+angle
            int = int+2
            
    #For Even Pole Amounts
    if((n%2) == 0):
        angle = np.pi/(n)
        temp_angle = angle/2
        den = np.array([1, 2*wc*np.cos(temp_angle), wc*wc])
        int=2
        while(int < n):
            temp_angle = temp_angle+angle
            den = np.convolve(den,np.array([1, 2*wc*np.cos(temp_angle), wc*wc]))
            int = int+2

    #Save the Lowpass filter num/den for reference
    denLP = den
    numLP = num

    #Highpass Filter Converstion

    if(mode==2):
        int = 1
        while(int<=n):
            num = np.convolve(num,np.array([1/wc, 0]))
            int = int + 1
        denHP = den
        numHP = num
        
    #Bandpass Filter Conversion
    if(mode==3):

        int = 1
        BP = np.array([1, 0, W*W])
        temp = 1
        while(int <=n):
            num = np.convolve(num, np.array([1, 0]))
            int= int + 1
        
        
        a = np.zeros((n+1,3+2*(n-1)))
        temp_arr = np.zeros(3+2*(n-1)-2)
        a[0] = np.concatenate([temp_arr,[0, 1]], axis=0)
        #Calculating the S -> S**2 + wc**2 coefficients
        for i in range(0,n):
            temp = np.convolve(temp, BP)
            temp_cat = np.concatenate((np.zeros((1,len(a.T)-len(temp))),temp), axis=None)
            a[i+1]= temp_cat
        
        temp_sum = np.zeros((n+1,1+3*(n)))
        #multiplying the temp_cat coefficients by the correct denominator value shifted
        for i in range(0,n+1):
            temp_sum[i] = np.convolve(den[len(den)-i-1]*a[i], np.concatenate(((np.zeros((1,i)), 1, np.zeros((1, n-i)))), axis=None))
       
        den = np.sum(temp_sum, axis=0)
        
        
    #------------------
    [f, h] = sig.freqs(num,den,worN=np.logspace(np.log10(mn), np.log10(mx), 10000))      
    #---------------------

    plt.subplot(311)
    plt.seminp.logx(f,20*np.log10(h))
    #plt.axis([mn, mx, -60, 10])
    plt.ylabel('|H| dB')
    plt.yticks([-40,-20,-3,0])
    plt.axvline(1,color='k')
    plt.axvline(5.5e-1,color='k')
    plt.axvline(2,color='k')
    plt.title('My_Fbode')
    plt.grid(which='both')


    #Conversion to Z Domain
    #s => (1-z**-1)/dt
    a2 = np.zeros((n+1,3*n+1))
    temp_arr2 = np.zeros(1+3*n-2)
    a2[0] = np.concatenate([temp_arr2,[0, 1]], axis=0)
    BP = [1/dt, -1/dt]
    y = np.zeros((n+1,3*n+1))
    tempZ = 1
    temp_arr = np.zeros(3*n-1)
    y[0] = np.concatenate([temp_arr,[0, 1]], axis=0)

    #Calculating S-> Z/dt-1/dt
    for i in range(0,n):
        tempZ = np.convolve(tempZ, BP)
        temp_cat = np.concatenate((np.zeros((1,len(a2.T)-len(tempZ))),tempZ), axis=None)
        a2[i+1]= temp_cat

    temp_sumd = np.zeros((n+1,4*n+1))
    temp_sumn = np.zeros((n+1,4*n+1))

    #Calculating the products of temp_cat, and their respectively shifted S values for num and den
    for i in range(0,n+1):
        temp_sumd[i] = np.convolve(den[len(den)-i-1]*a2[i], np.concatenate(((np.zeros((1,i)), 1, np.zeros((1, n-i)))), axis=None))
        temp_sumn[i] = np.convolve(num[len(num)-i-1]*a2[i], np.concatenate(((np.zeros((1,i)), 1, np.zeros((1, n-i)))), axis=None))
        
    den = np.sum(temp_sumd, axis=0)
    num = np.sum(temp_sumn, axis=0)

    #-----------------------
    [fz, hz] = sig.freqz(num,den,worN=10000) 
    #----------------------- 

    plt.subplot(313)
    plt.seminp.logx(fz,20*np.log10(hz))
    #plt.plot(fz,20*np.log10(hz))
    #plt.axis([1e-3, np.pi, -20, 0.1])
    plt.ylabel('|H| dB')
    plt.yticks([-40,-20,-3,0])
    #plt.axvline(30,color='k')
    #plt.text(50,-15,'$\phi$ = {}'.format(30,fontsize=12))
    plt.axvline(1*dt,color='k')
    plt.axvline(5.5e-1*dt,color='k')
    plt.axvline(2*dt,color='k')
    plt.title('My_Zbode')
    plt.grid(which='both')

# End of FILTER.PY