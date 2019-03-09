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
#   Jeremy Perhac
#   Daniel Allen
#
#   Special thanks to stackOverflow user: gg349 whose work is well documented
#   and used for the FFT calculations used throughout this library.
#
#   Included Functions
#   - FFT Coefficient Calculator:          fft_coef
#   - FFT Plotting Function:               fft_plot
#   - RMS Calculator:                      rms
#   - Step Function:                       step
#   - Multi-Argument Convolution:          convolve
#   - Convolution Bar Graph Visualizer:    convbar
#
#   Private Functions ( Those not Intended for Use Outside of Library )
#   - Tupple to Matrix Converter:       tuple_to_matrix
#   - Numpy Array to Matrix Converter:  nparr_to_matrix
#
#   Private Classes ( Those not Intended for Use Outside of Library )
#   - Function Concatinator:            c_func_concat
#
#   Constants
#   - NaN (Not A Number):                NAN
#   - Type: Numpy Matrix:                matrix
#   - Type: Tuple:                       tuple
#   - Type: Numpy Array:                 ndarr
#   - Type: Integer:                     tint
#   - Type: Float:                       tfloat
#   - Type: Function Handle:             tfun
#
#   Submodules
#   - Bode Plot Generator               BODE.PY         Imported as: *
#   - Filter Simulations                FILTERSIM.PY    Imported as: *
#   - Filter Operations/Tools           FILTER.PY       Imported as: filter
#################################################################################
name = "eesignal"
ver = "2.11.1"

# Import Submodules as Internal Functions
from .bode import *
from .filtersim import *
# Import Submodules as External Functions
from . import filter

# Import necessary libraries
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad as integrate
from scipy import signal as sig

# Define constants
u = 1e-6 # Micro (mu) Multiple
m = 1e-3 # Mili Multiple
k = 1e+3 # Kili Multiple
M = 1e+6 # Mega Multiple
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
    tuple:        Tuple of terms to be convolved.
                i.e. ( [1, 2], [3, 4], ..., [n-1, n] )
    
    Returns
    -------
    c:            The convolved set of the individual terms.
                i.e. np.array([ x1, x2, x3, ..., xn ])
    """
    c = sig.convolve(tuple[0],tuple[1])
    if (len(tuple) > 2):
        # Iterate starting with second element and continuing
        for i in range(2,len(tuple)):
            c = sig.convolve(c,tuple[i])
    return(c)

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


# End of __INIT__.PY