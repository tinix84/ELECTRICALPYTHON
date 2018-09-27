#################################################################################
#   EESignal.py
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
#   This library was compiled by Joe Stanley, special thanks to stackOverflow
#   user: gg349 whose work is well documented and used for the FFT calculations
#   used throughout this library.
#
#   Included Functions
#   - FFT Coefficient Calculator:		fft_coef
#   - FFT Plotting Function:			fft_plot
#   - RMS Calculator:					rms_calc
#   - State Space Simulator:			st_space
#   - Step Function						u
#   - Phase Margin:						pm
#   - Gain Margin:						gm
#   - System Response Plotter:			sys_response
#   - Multi-Argument Convolution:		convolve
#   - System Bode Plot:					bode
#
#   Private Functions ( Those not Intended for Use Outside of Library )
#   - TF System Conditioning:			sys_condition
#   - Tupple to Matrix Converter:		tuple_to_matrix
#   - Numpy Array to Matrix Converter:	nparr_to_matrix
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
#################################################################################

# Import necessary libraries
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad as integrate
from scipy import signal as sig

# Define constants
NAN = float('nan')
matrix = "<class 'numpy.matrixlib.defmatrix.matrix'>"
tuple = "<class 'tuple'>"
ndarr = "<class 'numpy.ndarray'>"
tint = "<class 'int'>"
tfloat = "<class 'float'>"
tfun = "<class 'function'>"

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

# Define System Bode Plotting Function
def bode(system,mn=-2,mx=3,npts=100,gtitle="",xlim=False,ylim=False):
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
	xlim:		Limit in x-axis for graph plot.
	ylim:		Limit in y-axis for graph plot.
	
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
	plt.title("Magnitude "+gtitle)
	plt.plot(w, mag)
	plt.xscale("log")
	plt.grid(which="both")
	plt.ylabel("Magnitude (dB)")
	plt.xlabel("Frequency (rad/sec)")
	if xlim!=False:
		plt.xlim(xlim)
	if ylim!=False:
		plt.ylim(ylim)
	plt.show()

	# Plot Angle
	plt.title("Angle "+gtitle)
	plt.plot(w, ang)
	plt.xscale("log")
	plt.grid(which="both")
	plt.ylabel("Angle (degrees)")
	plt.xlabel("Frequency (rad/sec)")
	if xlim!=False:
		plt.xlim(xlim)
	if ylim!=False:
		plt.ylim(ylim)
	plt.show()
	
	

# Define System Response Plotter function
def sys_response(system,npts=1000,dt=0.01,combine=True,gtitle="",
				stepResponse=True,rampResponse=False,parabolicResponse=False):
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
		plt.title("Step Response")
		plt.plot(TT,y1,'k--', label="Step Response")
		plt.plot(TT,step,'k', label="Step Function")
		plt.grid()
		plt.legend()
		plt.subplot(122)
		plt.title("Step Response Error")
		plt.plot(TT,errS,'k', label="Error")
		plt.grid()
		plt.legend()
		plt.subplots_adjust(wspace=0.3)
		plt.show()
	if (rampResponse):
		plt.figure()
		plt.subplot(121)
		plt.title("Ramp Response")
		plt.plot(TT,y2,'k--', label="Ramp Response")
		plt.plot(TT,ramp,'k', label="Ramp Function")
		plt.grid()
		plt.legend()
		plt.subplot(122)
		plt.title("Ramp Response Error")
		plt.plot(TT,errR,'k', label="Error")
		plt.grid()
		plt.legend()
		plt.subplots_adjust(wspace=0.3)
		plt.show()
	if (parabolicResponse):
		plt.figure()
		plt.subplot(121)
		plt.title("Parabolic Response")
		plt.plot(TT,y3,'k--', label="Parabolic Response")
		plt.plot(TT,parabola,'k', label="Parabolic Function")
		plt.grid()
		plt.legend()
		plt.subplot(122)
		plt.title("Parabolic Response Error")
		plt.plot(TT,errP,'k', label="Error")
		plt.grid()
		plt.legend()
		plt.subplots_adjust(wspace=0.3)
		plt.show()

# Define Gain Margin Calculator Function
def gm(tf,mn=-2,mx=3,npts=100,err=1e-12,printout=False,ret=True):
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
			if (ang[i] > -180):
				isP = True             		# Positive value found
				Pi = w[i]              		# Store frequency for positive value
				if abs(ang[i]+180) < err:	# if less than error, stop
					valid = False      		# stop while loop
					wg = Pi            		# store Phase Margin angle
					gm = mag[i]        		# store Phase Margin
					break              		# break out of for loop
			elif (ang[i] < -180):
				isN = True            		# Negative value found
				Ni = w[i]              		# Store frequency for negative value
				if abs(ang[i]+180) < err:   # if less than error, stop
					valid = False      		# stop while loop
					wg = Ni            		# store gain Margin angle
					gm = mag[i]        		# store gain Margin
					break              		# break out of for loop

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
def pm(tf,mn=-2,mx=3,npts=100,err=1e-12,printout=False,ret=True):
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
			if (mag[i] > 0):
				isP = True             # Positive value found
				Pi = w[i]              # Store frequency for positive value
				if abs(mag[i]) < err:  # if less than error, stop
					valid = False      # stop while loop
					wp = Pi            # store Phase Margin angle
					pm = ang[i]        # store Phase Margin
					break              # break out of for loop
			elif (mag[i] < 0):
				isN = True             # Negative value found
				Ni = w[i]              # Store frequency for negative value
				if abs(mag[i]) < err:  # if less than error, stop
					valid = False      # stop while loop
					wp = Ni            # store Phase Margin angle
					pm = ang[i]        # store Phase Margin
					break              # break out of for loop

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

# Define u(t) - Step function
def u(t):
	if (t>0):
		return(1)
	else:
		return(0)

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

def st_space(A,B,x=0,f=0,solution=2,C=False,D=False,npts=9999,NN=10000,dt=0.01,
		xlim=False,ylim=False,gtitle="",ret=False,plot=True,pltfn=False,sv=False):
	""" Plots the state-space simulation of an arbitrary set of matricies.

	Required Arguments:
	--------------------
	A :			Matrix A; if not type=numpy.matrix, converts to numpy.matrix
	B :			Matrix B; if not type=numpy.matrix, converts to numpy.matrix

	Optional Arguments:
	--------------------
	x :			Matrix x; if not type=numpy.matrix, converts to numpy.matrix
	fn :		Forcing Function; must be provided as callable function that
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

	Returns:
	--------
	If plot=True:	Generates and displays a plot of state-space simulation
	If ret=True:	Returns X-Axis variable and each state space term
					ex: ( x_axis, ( x1, x2, ... , xn ) )
					State-Space Variables return as tuple.

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
	elif (mF==tint) or (mF==tfloat): # If function returns int or float
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
		#print(fn_arr)
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
def rms_calc(f, T):
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
def fft_coef(f, T, N, return_complex=False):
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
def fft_plot(f, T, N, mn=False, mx=False, fftplot=True, absolute=False, title=False):
	""" Plots the FFT of the provided function as a stem plot.

	Arguments
	----------
	f : the periodic function, a callable like f(t)
	T : the period of the function f, so that f(0)==f(T)
	N_max : the function will return the first N_max + 1 Fourier coeff.
	mn : the minimum of the original signal
	mx : the maximum of the original signal

	Returns
	-------
	if fftplot=True, the function returns:
	Plot of FFT output of function

	if genSignal=True, the function returns:
	Approximation of original signal from FFT results

	if absolute=True, the function will:
	Return absolute values of the coefficients

	"""

	# Calculate FFT and find coefficients
	a0, a, b = fft_coef(f, T, N)

	# If provided a title, add it to the title string
	tStr = ""
	if title!=False:
		tStr = title

	# Plot FFT results with respect to their sign
	if fftplot and not absolute:
		# Set up Arguments
		rng = range(1,len(a)+1,1)
		xtic = range(0,len(a)+1,1)
		a0x = [0,0]
		a0y = [0,a0/2]

		# Plot
		plt.title("Fourier Coefficients"+tStr)
		plt.plot(a0x,a0y,'g')
		plt.stem(rng,a,'r','ro')
		plt.stem(rng,b,'b','bo')
		plt.xticks(xtic)
		plt.show()

	# Plot absolute value of FFT results
	if fftplot and absolute:
		# Set up Arguments
		rng = range(1,len(a)+1,1)
		xtic = range(0,len(a)+1,1)
		a0x = [0,0]
		a0y = [0,abs(a0)/2]

		# Plot
		plt.title("Fourier Coefficients"+tStr)
		plt.plot(a0x,a0y,'g')
		plt.stem(rng,np.abs(a),'r','ro')
		plt.stem(rng,np.abs(b),'b','bo')
		plt.xticks(xtic)
		plt.show()

	# Plot original function as described by FFT results
	if mx!=False:
		# Create domain variable
		x = np.arange(mn,mx,(mx-mn)/1000)
		# Set output to DC constant
		yout = a0
		# Plot each iteration of the Fourier Series
		for k in range(1,N):
			yout = yout + a[k-1]*np.cos(k*2*np.pi*x/T) + b[k-1]*np.sin(k*2*np.pi*x/T)
			plt.plot(x,yout)
		plt.title("Fourier Series Summation"+tStr)
		plt.show()
