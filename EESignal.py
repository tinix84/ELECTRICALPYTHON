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
#
#   This library was compiled by Joe Stanley, special thanks to stackOverflow
#   user: gg349 whose work is well documented and used here.
#
#   Included Functions
#   - FFT Coefficient Calculator:		fft_coef
#   - FFT Plotting Function:			fft_plot
#   - RMS Calculator:					rms_calc
#   - State Space Simulator:			st_space
#################################################################################

# Import necessary libraries
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad as integrate

def st_space(A,B,x,fn=0,solution=2,nsteps=6000,axis=[0,10,-1,2],gtitle=""):
    """ Plots the state-space simulation of an arbitrary set of matricies.
	
	Parameters:
	-----------
	A :			Matrix A
	B :			Matrix B
	x :			Matrix x
	fn :		Forcing Function
	solution:	Determines What Type of Solution Simulation is Required;
				Default of solution is 2
				0=zero-input
				1=zero-state
				2=total
	nsteps: 	Changes the range of simulation; defualt=6000
	axis:		Defines Plot Axes; [ x-min, x-max, y-min, y-max ]
				default=[0,10,-1,2]
	gtitle:		Additional String for Plot Title
	
	Returns:
	--------
	Generates and displays a plot of state-space simulation
	
	"""
    
    # Test for matrix inputs
    matrix = "<class 'numpy.matrixlib.defmatrix.matrix'>"
    mA = str(type(A))
    mB = str(type(B))
    mx = str(type(x))
    if (mA!=matrix) or (mB!=matrix) or (mx!=matrix):
        raise ValueError("ERROR: All Arguments must be Numpy Matrix of type: "+
                        matrix+"\nOne or more arguments do not match criteria."+
                        "\n\nArguments:\nA: "+mA+"\nB: "+mB+"\nx: "+mx)
	# Test for inputs
	if (fn==0) and (solution!=0):
		raise ValueError("ERROR: Arguments indicate no forcing function specified"+
						"yet solution has been requested in a form where forcing function"+
						"is required. Please check arguments.")
    
    # Start by defining Constants
    NN = 10000
    #change nsteps this to change the range of simulation
    dt = 0.01
    T = 0
    TT = np.arange(0,(dt*(NN)),dt)
    
    # Create a keyed list of state-space variables
    xtim = {}
    xtim_len = x.shape[0]
    for n in range(xtim_len):
        key = n #Each key should be the iterative variable
        xtim_init = np.zeros(NN) #Define the initial array
        xtim[key] = xtim_init #Create each xtim
    
    # When asked to find zero-state, set all ICs to zero
    if solution == 1:
        for n in range(xtim_len):
            x[n] = 0 #Set each value to zero
    
    # Finite-Difference Simulation
    for i in range(0,nsteps):
        for n in range(xtim_len):
            xtim[n][i] = x[n] #xtim[state-space][domain] = x[state-space]

        if solution == 0: #Zero-input, no added function input
            x = x + dt*A*x
        else: #Zero-state or Total, add function input
            x = x + dt*A*x + dt*B*fn(T)
        
        T = T+dt #Add discrete increment to T
        
    # Plot each state-variable over time
    for x in range(xtim_len):
        plt.plot(TT,xtim[x],label="x"+str(x+1))
    plt.axis(axis)
    plt.title("Simulated Output"+gtitle)
    plt.xlabel("Time (seconds)")
    plt.legend(title="State Space")
    plt.show()

# RMS Calculating Function
def rms_calc(f, T):
	""" Calculates the RMS value of the provided function.
	
	Parameters
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

    Parameters
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
	
	Parameters
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
		# Set up parameters
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
		# Set up parameters
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

	