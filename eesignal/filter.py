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
#   Steven Weeks
#   Jeremy Perhac
#   Daniel Allen
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
#   - Complex Complete-Square Terms     completesquare
#   - Down-Sampler                      dnsample
#   - Up-Sampler                        upsample
#   - Nth-Order Quadrature-Mirror       quadmirror
#   - Quad-Mirror Transfer Builder      quadtransfers
#   - Autocorrelation Function:         autocorr
#
#   Private Functions ( Those not Intended for Use Outside of Library )
#   - TF System Conditioning:           sys_condition
#################################################################################

# Import Required Libraries
import numpy as np
import scipy.signal as sig
import sympy as sym
from sympy.abc import s as s
import matplotlib.pyplot as plt

# Define Butterworth Minimum Order Solver:
def but_minord(mxDev, w, wc=1):
    """ Minimum Order Solving Function
    
    Finds the minimum order allowable to create the butterworth
    filter to match the specified criteria.
    
    Arguments:
    ----------
    mxDev:        The maximum allowable deviation.
    w:            The frequency at which mxDev occurs.
    wc:            The cuttoff frequency (default is 1).
    
    Returns:
    --------
    order:        An integer representing the minimum order."""
    
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
    n:        Filter Order
    a:        Calculated descriptive term "a"
    b:        Calculated descriptive term "b"
    
    Optional Arguments:
    -------------------
    type:    The Filter type, either 1 or 2 denoting
            Chebyshev type I or type II.
            
    Returns:
    --------
    totPole:    The complete system of poles returned as a
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
    n:        Filter Order
            
    Returns:
    --------
    wk:        Each omega returned as a list
    zeros:    The complete system of zeros returned as a
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
    ws:        Stop-Band Frequency
    Hs:        Stop-Band Magnitude
    Hp:        Pass-Band Magnitude
    
    Optional Arguments:
    -------------------
    n:        Filter Order, used to force and bypass calculation of n.
    
    Returns:
    ep:        Epsilon of system
    n:        System Order
    alpha:    Alpha of system
    a:        A constant of system
    b:        B constant of system
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
    wp:        Pass-Band Frequency
    Hs:        Stop-Band Magnitude
    Hp:        Pass-Band Magnitude
    
    Optional Arguments:
    -------------------
    n:        Filter Order, used to force and bypass calculation of n.
    
    Returns:
    ep:        Epsilon of system
    n:        System Order
    alpha:    Alpha of system
    a:        A constant of system
    b:        B constant of system
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
    poles:    The numpy array of polynomial factors.
    
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
    if ( len(system) == 2 ):        # System found to be num and den
        num = system[0]
        den = system[1]
        # Convolve numerator or denominator as needed
        if (str(type(num)) == tuple):
            num = convolve(num)        # Convolve terms in numerator
        if (str(type(den)) == tuple):
            den = convolve(den)        # Convolve terms in denominator
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
                num = num[i:]        # Slice zeros off the front of the numerator
                break                 # Break out of for loop
        for i in range( len( den ) ):
            if (den[i] != 0):
                den = den[i:]        # Slice zeros off the front of the denominator
                break                 # Break out of for loop
        system = (num,den)  # Repack system
    return(system) # Return the conditioned system

# Define Filter to band-pass function
def convert( sys, son, sod=None, debug=False, TFprint=False):
    """ Filter Conversion Function
    
    Purpose: This function is developed to perform the polynomial
    shift and conversion as dictated by the inputs. This function
    is to be used to shift transfer function systems or to convert
    the system as prescribed.
    
    Required Arguments:
    -------------------
    sys:        The tuple of the system, (numerator, denominator)
    son:        The numerator of the conversion factor, provided
                as an array or list. (i.e. [ s^n, ... s^1, s^0 ])
    
    Optional Arguments:
    -------------------
    sod:        The denominator of the conversion factor, must be
                provided as an array or list,
                (i.e. [ s^n, ... s^1, s^0 ]), default=None
    debug:      Print debugging information, default=False
    TFprint:    Print the resulting transfer function, default=False
    
    Returns:
    --------
    num:        The newly updated numerator polynomial
    den:        The newly updated denominator polynomial
    
    """
    # Evaluate Symbolic Representation of Polynomials
    convn = 0
    for i in range(len(son)):
        convn = convn + s**i * son[-(i+1)] # Add symbolic terms, low-ord first
    if(sod!=None):
        convd = 0
        for i in range(len(sod)):
            convd = convd + s**i * sod[-(i+1)] # Add symbolic terms, low-ord first
    else: convd = 1*s**0 # No denominator conversion term provided
    
    # Condition Symbolic Conversion terms
    convn = sym.expand( convn )
    convd = sym.expand( convd )
    
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
    # Print debug information if needed
    if debug:
        print("Input Conversion Factors:", son, sod)
        print(len(son),len(sod))
        print("Polynomial Conversion Factors:", convn, convd)
        print("Numerator Conversion Process:", num1,num2)
        print("Denominator Conversion Process:", den1,den2)
        print(m)
        print(LC)
        print(final_den)
    den = np.asarray(final_den.all_coeffs()).astype(np.double)
    try: # If the length is greater than 1, no exception
        final_num = sym.poly(num3)
        num = np.asarray(final_num.all_coeffs()).astype(np.double)
    except: # Numerator only 1-element long
        num = np.asarray(num3).astype(np.double)

    # Print debug information if needed
    if debug:
        print(final_num)
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
    system:        The tupled system of (numerator, denominator)
    
    Optional Arguments:
    -------------------
    mn:            The minimum value (10^mn) to be calculated; default=-1.
    mx:            The maximum value (10^mx) to be calculated; default=3.
    npts:        The number of points to calculate over; default=1000.
    yticks:        An array of the points to plot as tickmarks; default=False.
    forceticks:    A value that will allow forced elimination of preexisting
                tickmarks and only leave those provided as yticks or xticks;
                default=False.
    gtitle:        A string to be presented as the plot title; default=""
    xlim:        The list representing the minimum and maximum values over
                which to plot in x-axis; default=False
    ylim:        Same as xlim, but for y-axis.
    ysize:        The font size for y-tickmarks; default is 10
    xticks:        Same as yticks, but for x-axis.
    xsize:        Same as ysize, but for x-axis.
    sv:            A value that will allow the saving of the plotted figure as
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
    tf:            The Transfer Function; can be provided as the following:
                - 1 (instance of lti)
                - 2 (num, den)
                - 3 (zeros, poles, gain)
                - 4 (A, B, C, D)
    
    Optional Arguments:
    -------------------
    mn:            The minimum frequency (as an exponent to 10, e.g. 10^mn)
                to be calculated for. Default is -2.
    mx:            The maximum frequency (as an exponent to 10, e.g. 10^mx)
                to be calculated for. Default is 3.
    npts:         The number of points over which to calculate the system.
                Default is 100.
    err:        The maximum allowable error for an aproximation of zero
                (i.e. the difference between found value and zero).
                Default is 1e-12.
    printout:    If set to true, will automatically print the values.
                Default is False.
    ret:        If set to true, will return the gain margin and frequency.
                Default is True.
    find:        The value to be searched for. For gain margin, the default
                is -180 (degrees).
    
    Returns:
    --------
    wg:            First returned argument; the frequency at which gain margin
                occurs (in radians per second e.g. rad/s).
    gm:            Second and final returned argument; Gain Margin (in dB).
    
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
                isP = True                         # Positive value found
                Pi = w[i]                          # Store frequency for positive value
                if abs(ang[i]-find) < err:        # if less than error, stop
                    valid = False                  # stop while loop
                    wg = Pi                        # store Phase Margin angle
                    gm = mag[i]                    # store Phase Margin
                    break                          # break out of for loop
            elif (ang[i] < find):
                isN = True                        # Negative value found
                Ni = w[i]                          # Store frequency for negative value
                if abs(ang[i]-find) < err:        # if less than error, stop
                    valid = False                  # stop while loop
                    wg = Ni                        # store gain Margin angle
                    gm = mag[i]                    # store gain Margin
                    break                          # break out of for loop

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
    tf:            The Transfer Function; can be provided as the following:
                - 1 (instance of lti)
                - 2 (num, den)
                - 3 (zeros, poles, gain)
                - 4 (A, B, C, D)
    
    Optional Arguments:
    -------------------
    mn:            The minimum frequency (as an exponent to 10, e.g. 10^mn)
                to be calculated for. Default is -2.
    mx:            The maximum frequency (as an exponent to 10, e.g. 10^mx)
                to be calculated for. Default is 3.
    npts:         The number of points over which to calculate the system.
                Default is 100.
    err:        The maximum allowable error for an aproximation of zero
                (i.e. the difference between found value and zero).
                Default is 1e-12.
    printout:    If set to true, will automatically print the values.
                Default is False.
    ret:        If set to true, will return the phase margin and frequency.
                Default is True.
    find:        The value to be searched for. Default for phase margin is
                0 (dB).
    
    Returns:
    --------
    wp:            First returned argument; the frequency at which phase margin
                occurs (in radians per second e.g. rad/s).
    pm:            Second and final returned argument; Phase Margin (in degrees).
    
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
                isP = True                         # Positive value found
                Pi = w[i]                          # Store frequency for positive value
                if abs(mag[i]-find) < err:        # if less than error, stop
                    valid = False                  # stop while loop
                    wp = Pi                        # store Phase Margin angle
                    pm = ang[i]                    # store Phase Margin
                    break                          # break out of for loop
            elif (mag[i] < find):
                isN = True                         # Negative value found
                Ni = w[i]                          # Store frequency for negative value
                if abs(mag[i]-find) < err:        # if less than error, stop
                    valid = False                  # stop while loop
                    wp = Ni                        # store Phase Margin angle
                    pm = ang[i]                    # store Phase Margin
                    break                          # break out of for loop

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
    system:                The Transfer Function; can be provided as the following:
                        - 1 (instance of lti)
                        - 2 (num, den)
                        - 3 (zeros, poles, gain)
                        - 4 (A, B, C, D)
    desired:            The goal Phase-Margin.
    
    Optional Arguments:
    -------------------
    tolerance:            The additional phase given to make bring
                        the output closer to the desired result.
                        Default is 5.
    printout:            Print out the all values in order of calculation.
                        Default is False.
    ret:                Determines if there are values returned.
                        Default is True.
    plot:                Determines if plots are generated, will generate
                        two graphs, original and corrected. Default is False.
                        
    Return:
    -------
    wp:                    The Pole value of Phase Lead circuit { G(s) }
    wz:                    The Zero value of Phase Lead circuit { G(s) }
    pm:                    The phase margin of the new system.
    
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
                default=None; example: freqs=[100,200,300]
    axis:       The bounds of the x- and y-axes; default=None;
                example: axis=[0,6500,-0.2,1.1]
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
        freqs.append(0)
    
    # Capture Input Function Over Range
    y[0] = 0
    for n in range(1,NN):
        x[n] = f(dt*n)
        
    # Generate FFT Decomposition
    X = (2/NN)*np.fft.fft(x)
    
    # Plot Original Function
    plt.figure()
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
        plt.xticks(freqs)
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
            else:
                nc = nc//DF
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
    plt.tight_layout()
    plt.grid()
    # Plot Filtered FFT
    plt.subplot(324)
    plt.plot(FF,H,'k--')
    plt.plot(FF,abs(Y),'k')
    if(axis!=None):
        plt.axis(axis)
    if(freqs!=None):
        plt.xticks(freqs)
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
    plt.axis([0,2*MM,-1,1.2])
    plt.xticks([0,MM,int(2*MM)])
    plt.title("Shifted Filter 'h'")
    plt.tight_layout()
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
    plt.figure()
    plt.subplot(311)
    plt.plot(w,'k')
    plt.grid()
    plt.subplot(312)
    plt.plot(z,'k')
    plt.grid()
    plt.subplot(313)
    plt.plot(FF,abs(Z),'k')
    plt.ylabel('|Z(w)|')
    plt.yticks([0,.1,.9,1.])
    plt.tight_layout()
    if(axis!=None):
        plt.axis(axis)
    if(freqs!=None):
        plt.xticks(freqs)
    plt.title('Final Filtered FFT of Output System')
    plt.grid()
    plt.subplots_adjust(hspace=0.8)
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

# Define Complex Completing-The-Square Terms Finder
def completesquare(system):
    """
    COMPLETESQUARE Function
    
    Purpose:
    --------
    This function is designed to decompose a complex
    laplacian-domain system into K1, K2, alpha, omega
    in the following structure:
    
       a1*s + a0              o                s + a
    ---------------  = K1------------- + K2--------------
    s^2 + b1*s + b0      (s+a)^2 + o^2     (s+a)^2 + o^2
    
    where: a = alpha, and o = omega.
    
    Required Arguments:
    -------------------
    system:    The Transfer Function; can be provided as the following:
                - 1 (instance of lti)
                - 2 (num, den)
                - 3 (zeros, poles, gain)
                - 4 (A, B, C, D)
                
    Returns:
    --------
    (K1,K2,alpha,omega): Set of constants as defined in Purpose
                         section above.
    """
    # Condition system
    system = sys_condition(system,False)
    
    # Split into Numerator and Denominator
    num, den = system
    
    # Determine if the numerator and denominator
    # meet the specified criteria.
    num_sz = num.size
    den_sz = den.size
    if( num_sz!=2 or den_sz!=3):
        raise ValueError("ERROR: Improper input system size."+
                         " Numerator:"+str(num_sz)+
                         " Denominator:"+str(den_sz))
    
    # Calculate Alpha and Omega
    alpha = den[1]/2
    omega = np.sqrt(den[2]-alpha**2)
    
    # Calculate K1 and K2
    K1 = (num[1]-num[0]*alpha)/omega
    K2 = num[0]
    
    # Return Values
    return(K1,K2,alpha,omega)

# Define Down-Sampler
def dnsample(inarr,n=2):
    """
    DNSAMPLE Function
    
    Purpose:
    --------
    This function is designed to perform the operation
    of down-sampling an array or list input for the
    purposes of filtering. Down-Sampling will remove
    every n-th term.
    
    Required Arguments:
    -------------------
    inarr:  The input array, must be list or array.
    
    Optional Arguments:
    -------------------
    n:      The down-sample rate, default=2
    
    Returns:
    --------
    outarr: The down-sampled array
    """
    # Condition Input
    inarr = np.asarray(inarr)
    # Down-Sample
    outarr = np.copy(inarr[::n])
    return(outarr)

# Define Up-Sampler
def upsample(inarr,n=2,trailing=False):
    """
    UPSAMPLE Function
    
    Purpose:
    --------
    This function is designed to perform the operation
    of up-sampling an array or list input for the
    purposes of filtering. Up-Sampling will add a n-1
    zeros (0) between every other term.
    
    Required Arguments:
    -------------------
    inarr:  The input array, must be list or array.
    
    Optional Arguments:
    -------------------
    n:      The up-sample rate, default=2
    
    Returns:
    --------
    outarr: The up-sampled array
    """
    # Condition Input
    inarr = np.asarray(inarr)
    # Up-Sample
    arrlen = len(inarr)
    zero = np.zeros(n-1)
    outarr = np.array([])
    for i in range(arrlen):
        if(i==(arrlen-1) and not trailing):
            outarr = np.append(outarr,inarr[i])
        else:
            ap = np.append(inarr[i],zero)
            outarr = np.append(outarr,ap)
    return(outarr)
    

# Define Quadrature Mirror Filter
def quadmirror(farray,filterset=None,showall=False,pltinout=False,ord=1,
               stem=True,ret=True,trailing=False,figsize=None,reduce=False):
    """
    QUADMIRROR Function
    
    Purpose:
    --------
    This function is designed to provide a 1st-order basic quadrature
    mirror application for the purposes of filtering and compression.
    
    Required Arguments:
    -------------------
    farray:     The Function Array to be utilized in the Quad-Mirror
    
    Optional Arguments:
    -------------------
    filterset:  Specifies the four filters to be used (h0, f0, h1, f1),
                must be specified as a tupple in the before mentioned
                order, default=None
    showall:    Control argument to enable plotting intermediate steps,
                default=False
    pltinout:   Control argument to enable plotting in- and out-put arrays,
                default=False
    stem:       Control argument to allow plotting points as stem-function,
                default=True
    ret:        Control argument to enable return of array, default=True
    figsize:    Control argument to force size of subplot sizes,
                default=None
    ord:        Order of Quad-Mirror Filter, default=2, must be >= 1
    reduce:     Control argument to eliminate (zero-out) the high-frequency
                datapoints (C1) from the system, used for simulation of
                data compression, default=False
    
    Returns:
    --------
    Function returns the computed output array, commonly reffered to as x-hat.
    Function may also return the down-sampled stage c0 and c1 for the
    highest order when the *reduce* argument is set to True. [y, c0, c1]
    """
    # Order must be greater than or equal to 1
    if( ord<1 ):
        print("WARNING: Order cannot be less than 1.")
        print("Setting *ord* to 1.")
        ord = 1
    # Plot input array
    if pltinout:
        plt.figure()
        if stem: # Plot as stem values
            plt.stem(farray,linefmt='k',
                     basefmt='k',markerfmt='k.')
        else: plt.plot(farray) # Plot as standard function
        plt.title("Input")
        plt.show()
    
    # Define default filter elements if None provided
    if(filterset==None):
        h0 = [1,1]
        h1 = [1,-1]
        f0 = [1,1]
        f1 = [-1,1]
    else:
        h0,f0,h1,f1 = filterset
    
    # Calculate the shifting terms
    nshift = np.zeros(ord+1)
    nshift[1] = 0
    for n in range(2,ord+1):
        nshift[n] = 2*nshift[n-1] + 1
    # Calculate the size of the filters
    N = len(h0)
    M = len(h1)
    ishift = int((N + M)/2 - 1)
    # Generate the Necessary Number of zeros
    shift = np.zeros(ord)
    for i in range(ord):
        shift[i] = int(nshift[i+1]*ishift)
    shift = np.flip(shift,0)
    
    # Convolve to generate intermediate terms
    theta0 = [np.convolve(farray,h0)]
    theta1 = [np.convolve(farray,h1)]
    # Define Empty Lists
    c0 = []
    c1 = []
    y0 = []
    y1 = []
    y  = []
    # Ascend the filter
    for branch in range(ord):
        # Downsample
        c0.append(dnsample(theta0[branch]))
        c1t = dnsample(theta1[branch])
        zer = int(shift[branch])
        c1.append(np.append(np.zeros(zer),c1t))
        if(branch<ord-1):
            # Re-Evaluate Theta by Convolving
            theta0.append(np.convolve(c0[branch],h0))
            theta1.append(np.convolve(c0[branch],h1))
    # Upsample
    u0 = [upsample(c0[ord-1],trailing=trailing)]
    # When Data Compression Requested, Zero-Out C1 Terms
    if reduce:
        u1 = [upsample(np.zeros(len(c1[ord-1])),trailing=trailing)]
    else:
        u1 = [upsample(c1[ord-1],trailing=trailing)]
    # Descend the filter
    for i in range(ord):
        # Determine Branch
        branch = ord-(2+i)
        # Convolve to generate final term set
        y0.append(np.convolve(u0[i],f0))
        y1.append(np.convolve(u1[i],f1))
        # Ensure Equal Length of Arrays
        zerolen = len(y0[i])-len(y1[i])
        y1[i] = np.append(y1[i],np.zeros(zerolen))
        # Sum output
        y.append((y0[i]+y1[i])/2)
        if(branch>=0):
            # Re-Evaluate U by Upsampling
            u0.append(upsample(y[i],trailing=trailing))
            # When Data Compression Requested, Zero-Out C1 Terms
            if reduce:
                u1.append(upsample(np.zeros(len(c1[branch])),trailing=trailing))
            else:
                u1.append(upsample(c1[branch],trailing=trailing))
    u0.reverse()
    u1.reverse()
    # Plot all intermediate steps
    if(showall):
        plt.figure()
        if(ord==1):
            if(figsize!=None): plt.figure(figsize=figsize)
            plots = np.array([theta0,theta1,c0,c1,u0,u1,y0,y1])
            label = np.array(["θ-0[n]","θ-1[n]","C-0","C-1",
                              "U-0","U-1","Y-0","Y-1"])
            # Iteratively generate subplots
            for plot in range(len(plots)):
                plt.subplot(4,2,plot+1)
                plt.grid(True)
                plt.title(label[plot])
                vals = plots[plot]
                if stem: # Plot as stem values
                    plt.stem(vals[0],linefmt='k',
                             basefmt='k',markerfmt='k.')
                else: plt.plot(vals[0]) # Plot as standard function
                plt.tight_layout()
        else:
            if(figsize!=None): plt.figure(figsize=figsize)
            # Iteratively Generate Subplots
            for n in range(ord):
                plt.subplot(ord,2,n*2+1)
                plt.grid(True)
                plt.title("Level "+str(n+2)+" Downsampled Stage")
                if stem: # Plot as stem values
                    plt.stem(c0[n],linefmt='r',
                             basefmt='k',markerfmt='r.')
                    plt.stem(c1[n],linefmt='b',
                             basefmt='k',markerfmt='b.')
                else:
                    plt.plot(c0[n])
                    plt.plot(c1[n])
                plt.subplot(ord,2,n*2+2)
                plt.grid(True)
                plt.title("Level "+str(n+2)+" Upsampled Stage")
                if stem: # Plot as stem values
                    plt.stem(u0[n],linefmt='r',label='Top Branch',
                             basefmt='k',markerfmt='r.')
                    plt.stem(u1[n],linefmt='b',label='Bottom Branch',
                             basefmt='k',markerfmt='b.')
                else:
                    plt.plot(u0[n],label='Top Branch')
                    plt.plot(u1[n],label='Bottom Branch')
                plt.tight_layout()
                plt.legend()
        plt.show()
    # Plot the final output of the quadrature-mirror
    if pltinout:
        yout = y[ord-1]
        plt.figure()
        if stem: # Plot as stem values
            plt.stem(yout,linefmt='k',basefmt='k',
                     markerfmt='k.')
        else: plt.plot(yout) # Plot as standard function
        plt.title("Output")
        plt.show()
    # Return the computed output array
    if ret:
        if(ord==1):
            return(y[0])
        else:
            if reduce:
                return(y[ord-1],c0[ord-1],c1[ord-1])
            else:
                return(y[ord-1])


# Define Quad-Mirror Transfer Function Builder
def quadtransfers(p,offset=0,complex=True,round=None):
    """
    QUADTRANSFERS Function
    
    Purpose:
    --------
    This function is designed to generate the necessary
    transfer function arrays for a quadrature mirror filter.
    
    Required Arguments:
    -------------------
    p:          The order of the filter arrays
    
    Optional Arguments:
    -------------------
    offset:     A setting value to arrange the poles into the two
                filters (i.e. H0/H1) a positive value will shift
                more poles from F0 to H0, a negative value will
                shift poles from H0 to F0; default=0
    complex:    Control argument to allow filters to be returned
                as complex values, default=True
    round:      Rounding value to be used to simplify filters,
                default=None
    
    Returns:
    --------
    [h0,f0,h1,f1] : Tuple of Numpy Arrays Representing Transfer Functions
    """
    # Corner case where p=0
    if p==0:
        return([1],[1],[1],[1])
    # Corner case where p=1
    if p==1:
        return([1,1],[1,1],[1,-1],[-1,1])
    # Generate Transfer Functions for any p>1
    s = [1,1]
    s2 = s
    # Generate Pascal's Triangle by Convolution
    for i in range(1,2*p):
        s = np.convolve(s,s2)
    # Generate the Matrix To Solve for the Filters
    b = np.zeros(p-1)
    for i in range(0,p-1):
        b[i] = -s[2*i+1]
    a = [[0 for i in range(p-1)] for i in range(p-1)]
    for m in range(1,p):
        for n in range(2*m-1):
            if(n >= p-1):
                temp = 2*p-n-4
                a[m-1][temp] = a[m-1][temp]+s[2*m-n-2]
            else:
                a[m-1][n] = s[2*m-n-2]
    # Solve Matrix to Generate Filters
    c = np.linalg.solve(a,b)
    q = np.concatenate(([1],c),axis=None)
    h0 = np.zeros(2*len(q)-1)
    for i in range(len(q)):
        h0[i] = q[i]
    for i in range(len(q)-1):
        h0[2*len(q)-i-2] = q[i]
    h0 = h0*(-1)**(p)
    f0 = s
    P0 = np.convolve(h0,f0)
    for i in range(0,len(P0)-1):
        if abs(P0[i]) < 10**-10:
            P0[i] = 0
    H = np.roots(h0)
    F = np.roots(f0)
    # Alter the arangment of the filters to meet offset
    if(offset > 0):     # Shift poles from F0 to H0 (increase H0 length)
        SPLIT = offset
        temp_H = np.zeros(SPLIT,dtype=np.complex)
        temp_F = np.zeros(len(F)-SPLIT,dtype=np.complex)
        for i in range(0,SPLIT):
            temp_H[i] = F[i]
        for i in range(0,len(F)-SPLIT):
            temp_F[i] = F[i+SPLIT]
        H2 = np.concatenate((H,temp_H))
        F2 = temp_F
    elif(offset < 0):   # Shift poles from H0 to F0 (decrease H0 length)
        SPLIT = abs(offset)
        temp_F = np.zeros(SPLIT,dtype=np.complex)
        temp_H = np.zeros(len(H)-SPLIT,dtype=np.complex)
        for i in range(0,SPLIT):
            temp_F[i] = H[i]
        for i in range(0,len(H)-SPLIT):
            temp_H[i] = H[i+SPLIT]
        F2 = np.concatenate((F,temp_F))
        H2 = temp_H
    else:               # No offset required
        H2 = H
        F2 = F
    #Using the new roots to calculate the new filter coefficents
    h0 = np.polynomial.polynomial.polyfromroots(H2)
    f0 = np.polynomial.polynomial.polyfromroots(F2)
    f0 = f0*(-1)**(p-1)
    # Negate Filter Arrays
    h0 = -1*h0
    f0 = -1*f0
    # Normalize Filter Arrays
    h0 = h0 / np.sum(h0)*2
    f0 = f0 / np.sum(f0)*2
    # Generate H1 and F1 from H0 and F0
    N = len(h0)
    f1 = np.zeros(N)
    for i in range(0,N):
        f1[i] = ((-1)**(i-1))*h0[i]
    M = len(f0)
    h1 = np.zeros(M)
    for i in range(0,M):
        h1[i] = (-(-1)**(i-1))*f0[i]
    # Condition the Filter Arrays Before Returning
    if not complex: # Real Part Only Requested
        h0 = h0.real
        h1 = h1.real
        f0 = f0.real
        f1 = f1.real
    if round!=None: # Requested as Rounded Values
        h0 = np.around(h0,decimals=round)
        h1 = np.around(h1,decimals=round)
        f0 = np.around(f0,decimals=round)
        f1 = np.around(f1,decimals=round)
    # Return the Final Filter Arrays as Tuple
    return(h0,f0,h1,f1)


# Define Autocorrelation Function
def autocorr(arr,dt=0.01,tau=True,freqs=False):
    """
    AUTOCORR Function
    
    Purpose:
    --------
    This function is designed to perform the autocorrelation
    of an array input and return that correlated array.
    
    Required Arguments:
    -------------------
    arr:        The input array to be autocorrelated
    
    Optional Arguments:
    -------------------
    dt:         The time-step used for the array,
                default=0.01
    tau:        Control argument to force function to return
                the corresponding tau value along with the
                autocorrelated, default=True
    freqs:      Control argument to force function to return
                FFT frequency domain array, default=False
    
    Returns:
    --------
    out:        List of arrays in the form:
                [ autocorr, tau (option), freqs (option) ]
    """
    MM = len(arr)//2 - 1 # Evaluate 1/2 array length
    Rxx = np.zeros(MM)
    for m in range(0,MM):
        for n in range(0,MM):
            Rxx[m] = Rxx[m] + arr[n+m]*arr[n]
    out = [Rxx]
    if tau:
        # User desires tau to be returned
        out.append(np.linspace(dt,dt*MM,MM)) # Evaluate tau
    if freqs:
        # User desires FFT frequencies to be returned
        MM = len(arr)//4
        del_FR = 1/(MM*dt)
        out.append(np.linspace(0,(MM-1)*del_FR//2,MM))
    return(out)

# End of FILTER.PY