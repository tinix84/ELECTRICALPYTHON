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
#   - Gaussian Function:                   gaussian
#   - Gaussian Distribution Calculator:    gausdist
#   - Probability Density Calculator:      probdensity
#   - Real FFT Evaluator:                  rfft
#   - Normalized Power Spectrum:           wrms
#   - Hartley's Data Capacity Equation:    hartleydata
#   - Shannon's Data Capacity Equation:    shannondata
#   - String to Bit-String Converter:      string_to_bits
#   - CRC Message Generator:               crcsender
#   - CRC Remainder Calculator:            crcremainder
#
#   Private Functions ( Those not Intended for Use Outside of Library )
#   - Tupple to Matrix Converter:       tuple_to_matrix
#   - Numpy Array to Matrix Converter:  nparr_to_matrix
#
#   Private Classes ( Those not Intended for Use Outside of Library )
#   - Function Concatinator:            c_func_concat
#
#   Constants
#   - Micro (mu) Multiple:               u
#   - Mili Multiple:                     m
#   - Kilo Multiple:                     k
#   - Mega Multiple:                     M
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
ver = "2.13.10"

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
p = 1e-12 # Pico Multiple
n = 1e-9 # Nano Multiple
u = 1e-6 # Micro (mu) Multiple
m = 1e-3 # Mili Multiple
k = 1e+3 # Kili Multiple
M = 1e+6 # Mega Multiple
NAN = float('nan')
matrix = "<class 'numpy.matrixlib.defmatrix.matrix'>"
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
    
    Parameters
    ---------
    tuple:      Tuple of terms to be convolved.
                i.e. ( [1, 2], [3, 4], ..., [n-1, n] )
    
    Returns
    -------
    c:          The convolved set of the individual terms.
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

# Define Gaussian Function
def gaussian(x,mu=0,sigma=1):
    """
    GAUSSIAN Function:
    
    Purpose:
    --------
    This function is designed to generate the gaussian
    distribution curve with configuration mu and sigma.
    
    Required Arguments:
    -------------------
    x:       The input (array) x.
    
    Optional Arguments:
    -------------------
    mu:      Optional control argument, default=0
    sigma:   Optional control argument, default=1
    
    Returns:
    --------
    Computed gaussian (array) of the input x
    """
    return( 1/(sigma * np.sqrt(2 * np.pi)) *
            np.exp(-(x - mu)**2 / (2 * sigma**2)) )

# Define Gaussian Distribution Function
def gausdist(x,mu=0,sigma=1):
    """
    GAUSSDIST Function:
    
    Purpose:
    --------
    This function is designed to calculate the generic
    distribution of a gaussian function with controls
    for mu and sigma.
    
    Required Arguments:
    -------------------
    x:       The input (array) x
    
    Optional Arguments:
    -------------------
    mu:      Optional control argument, default=0
    sigma:   Optional control argument, default=1
    
    Returns:
    --------
    Computed distribution of the gausian function at the
    points specified by (array) x
    """
    F = np.array([])
    try:
        lx = len(x) # Find length of Input
    except:
        lx = 1 # Length 1
        x = [x] # Pack into list
    for i in range(lx):
        x_tmp = x[i]
        # Evaluate X (altered by mu and sigma)
        X = (x_tmp-mu) / sigma
        # Define Integrand
        def integrand(sq):
            return( np.exp(-sq**2/2) )
        integral = integrate(integrand,np.NINF,X) # Integrate
        result = 1/np.sqrt(2*np.pi) * integral[0] # Evaluate Result
        F = np.append(F, result) # Append to output list
    # Return only the 0-th value if there's only 1 value available
    if(len(F)==1):
        F = F[0]
    return(F)

# Define Probability Density Function
def probdensity(func,x,x0=0,scale=True):
    """
    PROBDENSITY Function:
    
    Purpose:
    --------
    This function uses an integral to compute the probability
    density of a given function.
    
    Required Arguments:
    -------------------
    func:    The function for which to calculate the PDF
    x:       The (array of) value(s) at which to calculate
             the PDF
    
    Optional Arguments:
    -------------------
    x0:      The lower-bound of the integral, starting point
             for the PDF to be calculated over, default=0
    scale:   The scaling to be applied to the output,
             default=True
    
    Returns:
    --------
    sumx:    The (array of) value(s) computed as the PDF at
             point(s) x
    """
    sumx = np.array([])
    try:
        lx = len(x) # Find length of Input
    except:
        lx = 1 # Length 1
        x = [x] # Pack into list
    # Recursively Find Probability Density
    for i in range(lx):
        sumx = np.append(sumx,integrate(func,x0,x[i])[0])
    # Return only the 0-th value if there's only 1 value available
    if(len(sumx)==1):
        sumx = sumx[0]
    else:
        if(scale==True):
            mx = sumx.max()
            sumx /= mx
        elif(scale!=False):
            sumx /= scale
    return(sumx)

# Define Real FFT Evaluation Function
def rfft(arr,dt=0.01,absolute=True,resample=True):
    """
    RFFT Function
    
    Purpose:
    --------
    This function is designed to evaluat the real FFT
    of a input signal in the form of an array or list.
    
    Required Arguments:
    -------------------
    arr:        The input array representing the signal
    
    Optional Arguments:
    -------------------
    dt:         The time-step used for the array,
                default=0.01
    absolute:   Control argument to force absolute
                values, default=True
    resample:   Control argument specifying whether
                the FFT output should be resampled,
                or if it should have a specific
                resampling rate, default=True
    
    Returns:
    --------
    FFT Array
    """
    # Calculate with Absolute Values
    if absolute:
        fourier = abs(np.fft.rfft(arr))
    else:
        foruier = np.fft.rfft(arr)
    if resample==True:
        # Evaluate the Downsampling Ratio
        dn = int(dt*len(arr))
        # Downsample to remove unnecessary points
        fixedfft = filter.dnsample(fourier,dn)
        return(fixedfft)
    elif resample==False:
        return(fourier)
    else:
        # Condition Resample Value
        resample = int(resample)
        # Downsample to remove unnecessary points
        fixedfft = filter.dnsample(fourier,resample)
        return(fixedfft)

# Define Normalized Power Spectrum Function
def wrms(func,dw=0.1,NN=100,quad=False,plot=True,
         title="Power Density Spectrum",round=3):
    """
    WRMS Function:
    
    Purpose:
    --------
    This function is designed to calculate the RMS
    bandwidth (Wrms) using a numerical process.
    
    Required Arguments:
    -------------------
    func:      The callable function to use for evaluation
    
    Optional Arguments:
    -------------------
    dw:        The delta-omega to be used, default=0.1
    NN:        The total number of points, default=100
    quad:      Control value to enable use of integrals
               default=False
    plot:      Control to enable plotting, default=True
    title:     Title displayed with plot,
               default="Power Density Spectrum"
    round:     Control to round the Wrms on plot,
               default=3
    
    Returns:
    --------
    W:         Calculated RMS Bandwidth (rad/sec)
    """
    # Define omega
    omega = np.linspace(0,(NN-1)*del_w,NN)
    # Initialize Fraction Terms
    Stot = Sw2 = 0
    # Power Density Spectrum
    Sxx = np.array([])
    for n in range(NN):
        # Calculate Power Density Spectrum
        Sxx = np.append(Sxx,func(omega[n]))
        Stot = Stot + Sxx[n]
        Sw2 = Sw2 + (omega[n]**2)*Sxx[n]
    if(quad):
        def intf(w):
            return(w**2*func(w))
        num = integrate(intf,0,np.inf)[0]
        den = integrate(func,0,np.inf)[0]
        # Calculate W
        W = np.sqrt(num/den)
    else:
        # Calculate W
        W = np.sqrt(Sw2/Stot)
    Wr = np.around(W,round)
    # Plot Upon Request
    if(plot):
        plt.plot(omega,Sxx)
        plt.title(title)
        # Evaluate Text Location
        x = 0.65*max(omega)
        y = 0.80*max(Sxx)
        plt.text(x,y,"Wrms: "+str(Wr))
        plt.show()
    # Return Calculated RMS Bandwidth
    return(W)
        
# Define Hartley's Equation for Data Capacity
def hartleydata(BW,M):
    """
    hartleydata Function
    
    Function to calculate Hartley's Law,
    the maximum data rate achievable for
    a given noiseless channel.
    
    Parameters
    ----------
    BW:         float
                Bandwidth of the data channel.
    M:          float
                Number of signal levels.
    
    Returns:
    --------
    C:          float
                Capacity of channel (in bits per second)
    """
    C = 2*BW*np.log2(M)
    return(C)

# Define Shannon's Equation For Data Capacity
def shannondata(BW,S,N):
    """
    shannondata Function
    
    Function to calculate the maximum data
    rate that may be achieved given a data
    channel and signal/noise characteristics
    using Shannon's equation.
    
    Parameters
    ----------
    BW:         float
                Bandwidth of the data channel.
    S:          float
                Signal strength (in Watts).
    N:          float
                Noise strength (in Watts).
    
    Returns
    -------
    C:          float
                Capacity of channel (in bits per second)
    """
    C = BW*np.log2(1+S/N)
    return(C)

# Define CRC Generator (Sender Side)
def crcsender(data, key):
    """
    crcsender Function
    
    Function to generate a CRC-embedded
    message ready for transmission.
    
    Contributing Author Credit:
    Shaurya Uppal
    Available from: geeksforgeeks.org
    
    Parameters
    ----------
    data:       string of bits
                The bit-string to be encoded.
    key:        string of bits
                Bit-string representing key.
    
    Returns
    -------
    codeword:   string of bits
                Bit-string representation of
                encoded message.
    """
    # Define Sub-Functions
    def xor(a, b): 
        # initialize result 
        result = [] 

        # Traverse all bits, if bits are 
        # same, then XOR is 0, else 1 
        for i in range(1, len(b)): 
            if a[i] == b[i]: 
                result.append('0') 
            else: 
                result.append('1') 

        return(''.join(result))

    # Performs Modulo-2 division 
    def mod2div(divident, divisor):
        # Number of bits to be XORed at a time. 
        pick = len(divisor) 

        # Slicing the divident to appropriate 
        # length for particular step 
        tmp = divident[0 : pick] 

        while pick < len(divident): 

            if tmp[0] == '1': 

                # replace the divident by the result 
                # of XOR and pull 1 bit down 
                tmp = xor(divisor, tmp) + divident[pick] 

            else:   # If leftmost bit is '0' 

                # If the leftmost bit of the dividend (or the 
                # part used in each step) is 0, the step cannot 
                # use the regular divisor; we need to use an 
                # all-0s divisor. 
                tmp = xor('0'*pick, tmp) + divident[pick] 

            # increment pick to move further 
            pick += 1

        # For the last n bits, we have to carry it out 
        # normally as increased value of pick will cause 
        # Index Out of Bounds. 
        if tmp[0] == '1': 
            tmp = xor(divisor, tmp) 
        else: 
            tmp = xor('0'*pick, tmp) 

        checkword = tmp 
        return(checkword)
    
    # Condition data
    data = str(data)
    # Condition Key
    key = str(key)
    l_key = len(key)
   
    # Appends n-1 zeroes at end of data 
    appended_data = data + '0'*(l_key-1) 
    remainder = mod2div(appended_data, key) 
   
    # Append remainder in the original data 
    codeword = data + remainder 
    return(codeword)

# Define CRC Generator (Sender Side)
def crcremainder(data, key):
    """
    crcremainder Function
    
    Function to calculate the CRC
    remainder of a CRC message.
    
    Contributing Author Credit:
    Shaurya Uppal
    Available from: geeksforgeeks.org
    
    Parameters
    ----------
    data:       string of bits
                The bit-string to be decoded.
    key:        string of bits
                Bit-string representing key.
    
    Returns
    -------
    remainder: string of bits
                Bit-string representation of
                encoded message.
    """
    # Define Sub-Functions
    def xor(a, b): 
        # initialize result 
        result = [] 

        # Traverse all bits, if bits are 
        # same, then XOR is 0, else 1 
        for i in range(1, len(b)): 
            if a[i] == b[i]: 
                result.append('0') 
            else: 
                result.append('1') 

        return(''.join(result))

    # Performs Modulo-2 division 
    def mod2div(divident, divisor):
        # Number of bits to be XORed at a time. 
        pick = len(divisor) 

        # Slicing the divident to appropriate 
        # length for particular step 
        tmp = divident[0 : pick] 

        while pick < len(divident): 

            if tmp[0] == '1': 

                # replace the divident by the result 
                # of XOR and pull 1 bit down 
                tmp = xor(divisor, tmp) + divident[pick] 

            else:   # If leftmost bit is '0' 

                # If the leftmost bit of the dividend (or the 
                # part used in each step) is 0, the step cannot 
                # use the regular divisor; we need to use an 
                # all-0s divisor. 
                tmp = xor('0'*pick, tmp) + divident[pick] 

            # increment pick to move further 
            pick += 1

        # For the last n bits, we have to carry it out 
        # normally as increased value of pick will cause 
        # Index Out of Bounds. 
        if tmp[0] == '1': 
            tmp = xor(divisor, tmp) 
        else: 
            tmp = xor('0'*pick, tmp) 

        checkword = tmp 
        return(checkword)
    
    # Condition data
    data = str(data)
    # Condition Key
    key = str(key)
    l_key = len(key)
   
    # Appends n-1 zeroes at end of data 
    appended_data = data + '0'*(l_key-1) 
    remainder = mod2div(appended_data, key) 
   
    return(remainder)

def string_to_bits(str):
    """
    string_to_bits Function
    
    Converts a Pythonic string to the string's
    binary representation.
    
    Parameters
    ----------
    str:        string
                The string to be converted.
    
    Returns
    -------
    data:       string
                The binary representation of the
                input string.
    """
    data = (''.join(format(ord(x), 'b') for x in str))
    return(data)

# End of __INIT__.PY