###################################################################
#   EEPOWER.PY
#
#   A library of functions, constants and more
#   that are related to Power in Electrical Engineering.
#
#   Written by Joe Stanley
#
#   Special Thanks To:
#   Paul Ortmann - Idaho Power
#
#   Included Constants:
#   - Micro (mu) Multiple:                      u
#   - Mili Multiple:                            m
#   - Kilo Multiple:                            k
#   - Mega Multiple:                            M
#   - 'A' Operator for Symmetrical Components:  a
#   - Not a Number value (NaN):                 NAN
#
#   Symmetrical Components Matricies:
#   - ABC to 012 Conversion:        Aabc
#   - 012 to ABC Conversion:        A012
#
#   Included Functions
#   - Phasor V/I Generator:         phasor
#   - Phasor Impedance Generator:   phasorz
#   - Complex Display Function:     cprint
#   - Parallel Impedance Adder:     parallelz
#   - V/I Line/Phase Converter:     phaseline
#   - Power Set Values:             powerset
#   - Power Triangle Function:      powertriangle
#   - Transformer SC OC Tests:      transformertest
#   - Phasor Plot Generator:        phasorplot
#   - Total Harmonic Distortion:    thd
#   - Total Demand Distortion:      tdd
#   - Reactance Calculator:         reactance
#   - Non-Linear PF Calc:           nlinpf
#   - Harmonic Limit Calculator:    harmoniclimit
#   - Power Factor Distiortion:     pfdist
#   - Short-Circuit RL Current:     iscrl
#   - Voltage Divider:              voltdiv
#   - Current Divider:              curdiv
#   - Instantaneous Power Calc.:    instpower
#   - Delta-Wye Network Converter:  dynetz
#   - Single Line Power Flow:       powerflow
#
#   Additional functions available in sub-modules:
#   - passives.py (renamed from capacitor.py)
#   - fault.py
#   - electronics.py
#   - perunit.py
#   - systemsolution.py
###################################################################
name = "eepower"
ver = "2.2.1"

# Import Submodules
from .passives import *
from .perunit import *
from .systemsolution import *
# Import Submodules as External Functions
from . import fault
from . import electronics

# Import libraries as needed:
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import cmath as c

# Define constants
a = c.rect(1,np.radians(120)) # A Operator for Sym. Components
p = 1e-12 # Pico Multiple
n = 1e-9 # Nano Multiple
u = 1e-6 # Micro (mu) Multiple
m = 1e-3 # Mili Multiple
k = 1e+3 # Kili Multiple
M = 1e+6 # Mega Multiple
NAN = float('nan')
VLLcVLN = c.rect(np.sqrt(3),np.radians(30)) # Conversion Operator
ILcIP = c.rect(np.sqrt(3),np.radians(-30)) # Conversion Operator

# Define symmetrical components matricies
Aabc = 1/3 * np.array([[ 1, 1, 1    ],  # Convert ABC to 012
                       [ 1, a, a**2 ],  # (i.e. phase to sequence)
                       [ 1, a**2, a ]])
A012 = np.array([[ 1, 1, 1    ],        # Convert 012 to ABC
                 [ 1, a**2, a ],        # (i.e. sequence to phase)
                 [ 1, a, a**2 ]])

# Define type constants
matrix = "<class 'numpy.matrixlib.defmatrix.matrix'>"
ndarr = "<class 'numpy.ndarray'>"

# Define Phasor Generator
def phasor( mag, ang ):
    """
    phasor Function:
    
    Generates the standard Pythonic complex representation
    of a phasor voltage or current when given the magnitude
    and angle of the specific voltage or current.
    
    Parameters
    ----------
    mag:        float
                The Magnitude of the Voltage/Current
    ang:        float
                The Angle (in degrees) of the Voltage/Current
    
    Returns
    -------
    phasor:     complex
                Standard Pythonic Complex Representation of
                the specified voltage or current.
    """
    return( c.rect( mag, np.radians( ang ) ) )

# Define Reactance Calculator
def reactance(z,f=60,sensetivity=1e-12):
    """
    reactance Function:
    
    Calculates the Capacitance or Inductance in Farads or Henreys
    (respectively) provided the impedance of an element.
    Will return capacitance (in Farads) if ohmic impedance is
    negative, or inductance (in Henrys) if ohmic impedance is
    positive. If imaginary: calculate with j factor (imaginary number).
    
    Parameters
    ----------
    z:              complex
                    The Impedance Provided, may be complex (R+jI)
    f:              float, optional
                    The Frequency Base for Provided Impedance, default=60
    sensetivity:    float, optional
                    The sensetivity used to check if a resistance was
                    provided, default=1e-12
    
    Returns
    -------
    out:            float
                    Capacitance or Inductance of Impedance
    """
    # Evaluate Omega
    w = 2*np.pi*f
    # Input is Complex
    if isinstance(z, complex):
        # Test for Resistance
        if(abs(z.real) > sensetivity):
            R = z.real
        else:
            R = 0
        if (z.imag > 0):
            out = z/(w*1j)
        else:
            out = 1/(w*1j*z)
        out = abs(out)
        # Combine with resistance if present
        if(R!=0): out = (R, out)
    else:
        if (z > 0):
            out = z/(w)
        else:
            out = 1/(w*z)
        out = abs(out)
    # Return Output
    return(out)

# Define display function
def cprint(val,unit="",label="",printval=True,ret=False,round=3):
    """
    cprint Function
    
    This function is designed to accept a complex value (val) and print
    the value in the standard electrical engineering notation:
    
    **magnitude ∠ angle °**
    
    This function will print the magnitude in degrees, and can print
    a unit and label in addition to the value itself.
    
    Parameters
    ----------
    val:        complex
                The Complex Value to be Printed, may be singular value,
                tuple of values, or list/array.
    unit:       string, optional
                The string to be printed corresponding to the unit mark.
                default=""
    label:      string, optional
                The pre-pended string used as a descriptive labeling string.
                default=""
    printval:   bool, optional
                Control argument enabling/disabling printing of the string.
                default=True
    ret:        bool, optional
                Control argument allowing the evaluated value to be returned.
                default=False
    round:      int, optional
                Control argument specifying how many decimals of the complex
                value to be printed. May be negative to round to spaces
                to the left of the decimal place (follows standard round()
                functionality). default=3
    
    Returns
    -------
    numarr:     numpy.ndarray
                The array of values corresponding to the magnitude and angle,
                values are returned in the form: [[ mag, ang ],...,[ mag, ang ]]
                where the angles are evaluated in degrees.
    """
    printarr = np.array([]) # Empty array
    numarr = np.array([]) # Empty array
    # Find length of the input array
    try:
        len(val) # Test Case, For more than one Input
        val = np.asarray(val) # Ensure that input is array
        shp = val.shape
        if(len(shp) > 1):
            row, col = shp
        else:
            col = shp[0]
            row = 1
            val = val.reshape((col,row))
            col = row
            row = shp[0]
        sz = val.size
        mult = True
        # Handle Label for Each Element
        if label=="":
            label = np.array([])
            for _ in range(sz):
                label = np.append(label,[""])
        elif len(label)==1 or str(type(label))==tstr:
            tmp = label
            for _ in range(sz):
                label = np.append(label,[tmp])
        # Handle Unit for Each Element
        if unit=="":
            unit = np.array([])
            for _ in range(sz):
                unit = np.append(unit,[""])
        elif len(unit)==1 or str(type(unit))==tstr:
            tmp = unit
            for _ in range(sz):
                unit = np.append(unit,[tmp])
    except:
        row = 1
        col = 1
        sz = 1
        mult = False
        _label = label
        _unit = unit
    # For each value in the input (array)
    for i in range(row):
        if mult:
            _val = val[i]
            _label = label[i]
            _unit = unit[i]
        else:
            _val = val
            _label = label
            _unit = unit
        mag, ang_r = c.polar(_val) #Convert to polar form
        ang = np.degrees(ang_r) #Convert to degrees
        mag = np.around( mag, round ) #Round
        ang = np.around( ang, round ) #Round
        strg = _label+" "+str(mag)+" ∠ "+str(ang)+"° "+_unit
        printarr = np.append(printarr, strg)
        numarr = np.append(numarr, [mag, ang])
    # Reshape the array to match input
    printarr = np.reshape(printarr, (row,col))
    numarr = np.reshape(numarr, (sz, 2))
    # Print values (by default)
    if printval and row==1:
        print(strg)
    elif printval:
        print(printarr)
    # Return values when requested
    if ret:
        return(numarr)

# Define Impedance Conversion function
def phasorz(C=None,L=None,f=60,complex=True):
    """
    phasorz Function:
    
    This function's purpose is to generate the phasor-based
    impedance of the specified input given as either the
    capacitance (in Farads) or the inductance (in Henreys).
    The function will return the phasor value (in Ohms).
    
    Parameters
    ----------
    C:          float, optional
                The capacitance value (specified in Farads),
                default=None
    L:          float, optional
                The inductance value (specified in Henreys),
                default=None
    f:          float, optional
                The system frequency to be calculated upon, default=60
    complex:    bool, optional
                Control argument to specify whether the returned
                value should be returned as a complex value.
                default=True
    
    Returns
    -------
    Z:      complex
            The ohmic impedance of either C or L (respectively).
    """
    w = 2*np.pi*f
    #C Given in ohms, return as Z
    if (C!=None):
        Z = -1/(w*C)
    #L Given in ohms, return as Z
    if (L!=None):
        Z = w*L
    #If asked for imaginary number
    if (complex):
        Z *= 1j
    return(Z)

# Define Parallel Impedance Adder
def parallelz(*args):
    """
    parallelz Function:
    
    This function is designed to generate the total parallel
    impedance of a set (tuple) of impedances specified as real
    or complex values.
    
    Parameters
    ----------
    Z:      tuple of complex
            The tupled input set of impedances, may be a tuple
            of any size greater than 2. May be real, complex, or
            a combination of the two.
    
    Returns
    -------
    Zp:     complex
            The calculated parallel impedance of the input tuple.
    """
    # Gather length (number of elements in tuple)
    L = len(args)
    if L==1:
        Z = args[0] # Only One Tuple Provided
        try:
            L = len(Z)
            if(L==1):
                Zp = Z[0] # Only one impedance, burried in tuple
            else:
                # Inversely add the first two elements in tuple
                Zp = (1/Z[0]+1/Z[1])**(-1)
                # If there are more than two elements, add them all inversely
                if(L > 2):
                    for i in range(2,L):
                        Zp = (1/Zp+1/Z[i])**(-1)
        except:
            Zp = Z # Only one impedance
    else:
        Z = args # Set of Args acts as Tuple
        # Inversely add the first two elements in tuple
        Zp = (1/Z[0]+1/Z[1])**(-1)
        # If there are more than two elements, add them all inversely
        if(L > 2):
            for i in range(2,L):
                Zp = (1/Zp+1/Z[i])**(-1)
    return(Zp)

# Define Phase/Line Converter
def phaseline(VLL=None,VLN=None,Iline=None,Iphase=None,complex=False):
    """
    phaseline Function
    
    This function is designed to return the phase- or line-equivalent
    of the voltage/current provided. It is designed to be used when
    converting delta- to wye-connections and vice-versa.
    Given a voltage of one type, this function will return the
    voltage of the opposite type. The same is true for current.
    
    Parameters
    ----------
    VLL:        float, optional
                The Line-to-Line Voltage; default=None
    VLN:        float, optional
                The Line-to-Neutral Voltage; default=None
    Iline:      float, optional
                The Line-Current; default=None
    Iphase:     float, optional
                The Phase-Current; default=None
    complex:    bool, optional
                Control to return value in complex form; default=False
    
    ======  =======
    Inputs  Outputs
    ======  =======
    VLL     VLN
    VLN     VLL
    Iline   Iphase
    Iphase  Iline
    ======  =======
    """
    output = 0
    #Given VLL, convert to VLN
    if (VLL!=None):
        VLN = VLL/(VLLcVLN)
        output = VLN
    #Given VLN, convert to VLL
    elif (VLN!=None):
        VLL = VLN*VLLcVLN
        output = VLL
    #Given Iphase, convert to Iline
    elif (Iphase!=None):
        Iline = Iphase*ILcIP
        output = Iline
    #Given Iline, convert to Iphase
    elif (Iline!=None):
        Iphase = Iline/ILcIP
        output = Iphase
    #None given, error encountered
    else:
        print("ERROR: No value given"+
                "or innapropriate value"+
                "given.")
        return(0)
    #Return as complex only when requested
    if complex:
        return( output )
    return(abs( output ))

# Define Power Set Function
def powerset(P=None,Q=None,S=None,PF=None):
    """
    powerset Function
    
    This function is designed to calculate all values
    in the set { P, Q, S, PF } when two (2) of the
    values are provided. The equations in this
    function are prepared for AC values, that is:
    real and reactive power, apparent power, and power
    factor.
    
    Parameters
    ----------
    P:      float, optional
            Real Power, unitless; default=None
    Q:      float, optional
            Reactive Power, unitless; default=None
    S:      float, optional
            Apparent Power, unitless; default=None
    PF:     float, optional
            Power Factor, unitless, provided as a
            decimal value, lagging is positive,
            leading is negative; default=None
    
    Returns
    -------
    P:      float
            Calculated Real Power Magnitude
    Q:      float
            Calculated Reactive Power Magnitude
    S:      float
            Calculated Apparent Power Magnitude
    PF:     float
            Calculated Power Factor
    """
    #Given P and Q
    if (P!=None) and (Q!=None):
        S = np.sqrt(P**2+Q**2)
        PF = P/S
        if Q<0:
            PF=-PF
    #Given S and PF
    elif (S!=None) and (PF!=None):
        P = abs(S*PF)
        Q = np.sqrt(S**2-P**2)
        if PF<0:
            Q=-Q
    #Given P and PF
    elif (P!=None) and (PF!=None):
        S = P/PF
        Q = Q = np.sqrt(S**2-P**2)
        if PF<0:
            Q=-Q
    else:
        raise ValueError("ERROR: Invalid Parameters or too few"+
                        " parameters given to calculate.")
    
    # Return Values!
    return(P,Q,S,PF)

# Define Power Triangle Function
def powertriangle(P=None,Q=None,S=None,PF=None,color="red",
                  text="Power Triangle",printval=False):
    """
    powertriangle Function
    
    This function is designed to draw a power triangle given
    values for the complex power system.
    
    Parameters
    ----------
    P:          float
                Real Power, unitless; default=None
    Q:          float
                Reactive Power, unitless; default=None
    S:          float
                Apparent Power, unitless; default=None
    PF:         float
                Power Factor, unitless, provided as a
                decimal value, lagging is positive,
                leading is negative; default=None
    color:      string, optional
                The color of the power triangle lines;
                default="red"
    text:       string, optional
                The title of the power triangle plot,
                default="Power Triangle"
    printval:   bool, optional
                Control argument to allow the numeric
                values to be printed on the plot,
                default="False"
    """
    # Calculate all values if not all are provided
    if( P==None or Q==None or S==None or PF==None):
        P,Q,S,PF = powerset(P,Q,S,PF)

    #Generate Lines
    Plnx = [0,P]
    Plny = [0,0]
    Qlnx = [P,P]
    Qlny = [0,Q]
    Slnx = [0,P]
    Slny = [0,Q]

    # Plot Power Triangle
    plt.figure(1)
    plt.title(text)
    plt.plot(Plnx,Plny,color=color)
    plt.plot(Qlnx,Qlny,color=color)
    plt.plot(Slnx,Slny,color=color)
    plt.xlabel("Real Power (W)")
    plt.ylabel("Reactive Power (VAR)")
    mx = max(abs(P),abs(Q))

    if P>0:
        plt.xlim(0,mx*1.1)
        x=mx
    else:
        plt.xlim(-mx*1.1,0)
        x=-mx
    if Q>0:
        plt.ylim(0,mx*1.1)
        y=mx
    else:
        plt.ylim(-mx*1.1,0)
        y=-mx
    if PF > 0:
        PFtext = " Lagging"
    else:
        PFtext = " Leading"
    text = "P:   "+str(P)+" W\n"
    text = text+"Q:   "+str(Q)+" VAR\n"
    text = text+"S:   "+str(S)+" VA\n"
    text = text+"PF:  "+str(abs(PF))+PFtext+"\n"
    text = text+"ΘPF: "+str(np.degrees(np.arccos(PF)))+"°"+PFtext
    # Print all values if asked to
    if printval:
         plt.text(x/20,y*4/5,text,color=color)
    plt.show()

# Define Transformer Short-Circuit/Open-Circuit Function
def transformertest(Poc=False,Voc=False,Ioc=False,Psc=False,Vsc=False,
               Isc=False):
    """
    transformertest Function
    
    This function will determine the non-ideal circuit components of
    a transformer (Req and Xeq, or Rc and Xm) given the test-case
    parameters for the open-circuit test and/or the closed-circuit
    test. Requires one or both of two sets: { Poc, Voc, Ioc }, or
    { Psc, Vsc, Isc }.
    All values given must be given as absolute value, not complex.
    All values returned are given with respect to primary.
    
    Parameters
    ----------
    Poc:    float, optional
            The open-circuit measured power (real power), default=None
    Voc:    float, optional
            The open-circuit measured voltage (measured on X),
            default=None
    Ioc:    float, optional
            The open-circuit measured current (measured on primary),
            default=None
    Psc:    float, optional
            The short-circuit measured power (real power), default=None
    Vsc:    float, optional
            The short-circuit measured voltage (measured on X),
            default=None
    Isc:    float, optional
            The short-circuit measured current (measured on X),
            default=None
    
    Returns
    -------
    {Req,Xeq,Rc,Xm}:    Given all optional args
    {Rc, Xm}:           Given open-circuit parameters
    {Req, Xeq}:         Given short-circuit parameters
    """
    SC = False
    OC = False
    # Given Open-Circuit Values
    if (Poc!=None) and (Voc!=None) and (Ioc!=None):
        PF = Poc/(Voc*Ioc)
        Y = c.rect(Ioc/Voc,-np.arccos(PF))
        Rc = 1/Y.real
        Xm = -1/Y.imag
        OC = True
    # Given Short-Circuit Values
    if (Psc!=None) and (Vsc!=None) and (Isc!=None):
        PF = Psc/(Vsc*Isc)
        Zeq = c.rect(Vsc/Isc,np.arccos(PF))
        Req = Zeq.real
        Xeq = Zeq.imag
        SC = True
    # Return All if Found
    if OC and SC:
        return(Req,Xeq,Rc,Xm)
    elif OC:
        return(Rc,Xm)
    elif SC:
        return(Req,Xeq)
    else:
        print("An Error Was Encountered.\n"+
                "Not enough arguments were provided.")

# Define Phasor Plot Generator
def phasorplot(phasor,title="Phasor Diagram",legend=False,bg="#d5de9c",radius=1.2,
               colors = ["#FF0000","#800000","#FFFF00","#808000","#00ff00","#008000",
                         "#00ffff","#008080","#0000ff","#000080","#ff00ff","#800080"]):
    """
    phasorplot Function
    
    This function is designed to plot a phasor-diagram with angles in degrees
    for up to 12 phasor sets. Phasors must be passed as a complex number set,
    (e.g. [ m+ja, m+ja, m+ja, ... , m+ja ] ).
    
    Parameters
    ----------
    phasor:     list of complex
                The set of phasors to be plotted.
    title:      string, optional
                The Plot Title, default="Phasor Diagram"
    legend:     bool, optional
                Control argument to enable displaying the legend, must be passed
                as an array or list of strings, default=False
    bg:         string, optional
                Background-Color control, default="#d5de9c"
    radius:     float, optional
                The diagram radius, default=1.2
    colors:     list of str, optional
                List of hexidecimal color strings denoting the line colors to use.
    """
    # Check for more phasors than colors
    numphs = len(phasor)
    numclr = len(colors)
    if numphs > numclr:
        raise ValueError("ERROR: Too many phasors provided. Specify more line colors.")
    
    # Force square figure and square axes looks better for polar, IMO
    width, height = matplotlib.rcParams['figure.figsize']
    size = min(width, height)
    # Make a square figure
    fig = plt.figure(figsize=(size, size))
    ax = fig.add_axes([0.1, 0.1, 0.8, 0.8], polar=True, facecolor=bg)
    ax.set_rmax(radius)
    plt.grid(True)
    
    # Plot the diagram
    plt.title(title+"\n")
    handles=np.array([]) # Empty array for plot handles
    for i in range(numphs):
        mag, ang_r = c.polar(phasor[i])
        if legend!=False:
            hand = plt.arrow(0,0,ang_r,mag,color=colors[i],label=legend[i])
            handles = np.append(handles,[hand])
        else: plt.arrow(0,0,ang_r,mag,color=colors[i])
    if legend!=False: plt.legend((handles),legend)
    plt.show()

###################################################################
#   Define Total Demand Distortion function
#
#   Returns the Total demand distortion given an array of the
#   distortion factors.
#
#   Harmonic array should contain the fundamental frequency h1
#   IL: Peak Demand load current (RMS) at fundamental frequency
###################################################################
def tdd(harmonics,IL):
    # Sum all terms of 2*fundamental and up
    sum = 0
    for h in range(1,len(harmonics)):
        sum += harmonics[h]**2
    
    # Take square-root of sum and divide by IL
    TDD = np.sqrt(sum)/IL
    return(TDD)

###################################################################
#   Define Total Harmonic Distortion function
#
#   Returns the Total harmonic distortion given an array of the
#   distortion factors.
#
#   Harmonic array should contain the fundamental frequency h1
#   PFdist: the distorted power factor, can be used to find thd
###################################################################
def thd(harmonics=False,PFdist=False):
    if(PFdist != False):
        # Use PFdistortion equation to find THD
        THD = np.sqrt(1/(PFdist**2)-1)
    else:
        # Sum all terms of 2*fundamental and up
        sum = 0
        for h in range(1,len(harmonics)):
            sum += harmonics[h]**2
        # Take Square Root of Sum
        sum = np.sqrt(sum)
        # Divide by magnitude of fundamental frequency
        THD = sum/harmonics[0]
    return(THD)

###################################################################
#   Define Distorted Power Factor function
#
#   Returns the distorted power factor value given I1 and IRMS or
#   the set of harmonic current magnitudes.
#
#   Ih array should contain the fundamental frequency h1
###################################################################
def pfdist(I1=False,IRMS=False,Ih=False):
    if (I1 != False and IRMS != False):
        # Find PFdist by using fundamental and RMS current
        PFdist = I1/IRMS
    else:
        # Find PFdist by using THD
        THD = thd(Ih) # Calculate THD
        PFdist = 1/np.sqrt(1+THD**2)
    
    return(PFdist)

# Define Non-Linear Power Factor Calculator
def nlinpf(PFtrue=False,PFdist=False,PFdisp=False):
    """
    nlinpf Function
    
    This function is designed to evaluate one of three unknowns
    given the other two. These particular unknowns are the arguments
    and as such, they are described in the representative sections
    below.
    
    Parameters
    ----------
    PFtrue:     float, exclusive
                The "True" power-factor, default=None
    PFdist:     float, exclusive
                The "Distorted" power-factor, default=None
    PFdisp:     float, exclusive
                The "Displacement" power-factor, default=None
    
    Returns
    -------
    {unknown}:  This function will return the unknown variable from
                the previously described set of variables.
    """
    if(PFtrue!=None and PFdist!=None and PFdisp!=None):
        raise ValueError("ERROR: Too many constraints, no solution.") 
    elif ( PFdist!=None and PFdisp!=None ):
        return( PFdist * PFdisp )
    elif ( PFtrue!=None and PFdisp!=None ):
        return( PFtrue / PFdisp )
    elif ( PFtrue!=None and PFdist!=None ):
        return( PFtrue / PFdist )
    else:
        raise ValueError("ERROR: Function requires at least two arguments.") 

###################################################################
#   Define Harmonic Current Limit function
#
#   Returns the limits of harmonic currents given the load
#   characteristics (Short-Circuit Current, Peak Demand Current).
#   Prints results when printout=True.
#
#   By default, prints the results, does not return values as an
#   numpy array. Returns this array when ret=True.
#
#   Compares to measured harmonic currents if Ih is provided.
#
#   N is the maximum harmonic term.
###################################################################
def harmoniclimit(Isc,IL,N=0,Ih=0,printout=True,ret=False):
    percent = 1/100 # Use for scaling
    Ir = Isc/IL # compute ratio
    if(Ir < 20):
        # Generate Harmonic Factors
        f1o = 4.0*percent
        f2o = 2.0*percent
        f3o = 1.5*percent
        f4o = 0.6*percent
        f5o = 0.3*percent
        tddL = 5.0*percent
        f1e = f1o * 25*percent
        f2e = f2o * 25*percent
        f3e = f3o * 25*percent
        f4e = f4o * 25*percent
        f5e = f5o * 25*percent
        
    elif(20 <= Ir and Ir < 50):
        # Generate Harmonic Factors
        f1o = 7.0*percent
        f2o = 3.5*percent
        f3o = 2.5*percent
        f4o = 1.0*percent
        f5o = 0.5*percent
        tddL = 8.0*percent
        f1e = f1o * 25*percent
        f2e = f2o * 25*percent
        f3e = f3o * 25*percent
        f4e = f4o * 25*percent
        f5e = f5o * 25*percent
        
    elif(50 <= Ir and Ir < 100):
        # Generate Harmonic Factors
        f1o = 10.0*percent
        f2o = 4.5*percent
        f3o = 4.0*percent
        f4o = 1.5*percent
        f5o = 0.7*percent
        tddL = 12.0*percent
        f1e = f1o * 25*percent
        f2e = f2o * 25*percent
        f3e = f3o * 25*percent
        f4e = f4o * 25*percent
        f5e = f5o * 25*percent
        
    elif(100 <= Ir and Ir < 1000):
        # Generate Harmonic Factors
        f1o = 12.0*percent
        f2o = 5.5*percent
        f3o = 5.0*percent
        f4o = 2.0*percent
        f5o = 1.0*percent
        tddL = 15.0*percent
        f1e = f1o * 25*percent
        f2e = f2o * 25*percent
        f3e = f3o * 25*percent
        f4e = f4o * 25*percent
        f5e = f5o * 25*percent
        
    else:
        # Generate Harmonic Factors
        f1o = 15.0*percent
        f2o = 7.0*percent
        f3o = 6.0*percent
        f4o = 2.5*percent
        f5o = 1.4*percent
        tddL = 20.0*percent
        f1e = f1o * 25*percent
        f2e = f2o * 25*percent
        f3e = f3o * 25*percent
        f4e = f4o * 25*percent
        f5e = f5o * 25*percent
    
    # Create empty array to return
    retArr = np.zeros(51)
    
    # Print out values
    if(printout):
        print("IEEE 519-2014 Distorted Current Limits:\n"+
                "---------------------------------------")
        for i in range(3, 50 + 1):
            if(3<=i and i<11):
                if(i%2): # Odd term
                    retArr[i] = f1o*IL
                else: # Even term
                    retArr[i] = f1e*IL
            elif(11<=i and i<17):
                if(i%2): # Odd term
                    retArr[i] = f2o*IL
                else: # Even term
                    retArr[i] = f2e*IL
            elif(17<=i and i<23):
                if(i%2): # Odd term
                    retArr[i] = f3o*IL
                else: # Even term
                    retArr[i] = f3e*IL
            elif(23<=i and i<35):
                if(i%2): # Odd term
                    retArr[i] = f4o*IL
                else: # Even term
                    retArr[i] = f4e*IL
            elif(35<=i and i<=50):
                if(i%2): # Odd term
                    retArr[i] = f5o*IL
                else: # Even term
                    retArr[i] = f5e*IL
            else:
                print("Internal Error Encountered!")
            print("Limit of "+str(i)+"th harmonic:", retArr[i],"A")
            if(N!=0 and N<=i):
                break
        print("---------------------------------------")
    
    # Comparison requested
    if(str(type(Ih)) == ndarr):
        maxr = min(len(retArr), len(Ih)+1)
        for k in range(3, maxr):
            if(retArr[k] < Ih[k-1]):
                print("Limit surpassed for "+str(k)+"th Harmonic term.")
    
    # Return values
    if(ret):
        return(retArr)

# Define Short-Circuit RL Current Calculator
def iscrl(V,Z,t=None,f=None,mxcurrent=True,alpha=None):
    """
    iscrl Function
    
    The Isc-RL function (Short Circuit Current for RL Circuit)
    is designed to calculate the short-circuit current for an
    RL circuit.
    
    Parameters
    ----------
    V:          float
                The absolute magnitude of the voltage.
    Z:          float
                The complex value of the impedance. (R + jX)
    t:          float, optional
                The time at which the value should be calculated,
                should be specified in seconds, default=None
    f:          float, optional
                The system frequency, specified in Hz, default=None
    mxcurrent:  bool, optional
                Control variable to enable calculating the value at
                maximum current, default=True
    alpha:      float, optional
                Angle specification, default=None
    
    Returns
    -------
    Opt 1 - (Irms, IAC, K):     The RMS current with maximum DC
                                offset, the AC current magnitude,
                                and the asymmetry factor.
    Opt 2 - (i, iAC, iDC, T):   The Instantaneous current with
                                maximum DC offset, the instantaneous
                                AC current, the instantaneous DC
                                current, and the time-constant T.
    Opt 3 - (Iac):              The RMS current without DC offset.
    """
    # Calculate omega, theta, R, and X
    if(f!=None): omega = 2*np.pi*f
    else: omega = None
    R = abs(Z.real)
    X = abs(Z.imag)
    theta = np.arctan( X/R )
    
    # If Maximum Current is Desired and No alpha provided
    if(mxcurrent and alpha==None):
        alpha = theta - np.pi/2
    elif(mxcurrent and alpha!=None):
        raise ValueError("ERROR: Inappropriate Arguments Provided.\n"+
                         "Not both mxcurrent and alpha can be provided.")
    
    # Calculate Asymmetrical (total) Current if t != None
    if(t!=None and f!=None):
        # Calculate RMS if none of the angular values are provided
        if(alpha==None and omega==None):
            # Calculate tau
            tau = t/(1/60)
            K = np.sqrt(1 + 2*np.exp(-4*np.pi*tau/(X/R)) )
            IAC = abs(V/Z)
            Irms = K*IAC
            # Return Values
            return(Irms,IAC,K)
        elif(alpha==None or omega==None):
            raise ValueError("ERROR: Inappropriate Arguments Provided.")
        # Calculate Instantaneous if all angular values provided
        else:
            # Convert Degrees to Radians
            omega = np.radians(omega)
            alpha = np.radians(alpha)
            theta = np.radians(theta)
            # Calculate T
            T = X/(2*np.pi*f*R) # seconds
            # Calculate iAC and iDC
            iAC = np.sqrt(2)*V/Z*np.sin(omega*t+alpha-theta)
            iDC = -np.sqrt(2)*V/Z*np.sin(alpha-theta)*np.exp(-t/T)
            i = iAC + iDC
            # Return Values
            return(i,iAC,iDC,T)
    elif( (t!=None and f==None) or (t==None and f!=None) ):
        raise ValueError("ERROR: Inappropriate Arguments Provided.\n"+
                         "Must provide both t and f or neither.")
    else:
        IAC = abs(V/Z)
        return(Iac)

# Define Voltage Divider Calculator
def voltdiv(Vin,R1,R2,Rload=None):
    """
    voltdiv Function
    
    This function is designed to calculate the output
    voltage of a voltage divider given the input voltage,
    the resistances (or impedances) and the load resistance
    (or impedance) if present.
    
    Parameters
    ----------
    Vin:    float
            The Input Voltage, may be real or complex
    R1:     float
            The top resistor of the divider (real or complex)
    R2:     float
            The bottom resistor of the divider, the one which
            the output voltage is measured across, may be
            either real or complex
    Rload:  float, optional
            The Load Resistor (or impedance), default=None
    
    Returns
    -------
    Vout:   float
            The Output voltage as measured across R2 and/or Rload
    """
    # Determine whether Rload is given
    if(Rload==None): # No Load Given
        Vout = Vin * R2 / (R1+R2)
    else:   # Load was given
        Rp = parallelz((R2,Rload))
        Vout = Vin * Rp / (R1+Rp)
    return(Vout)

# Define Current Divider Calculator
def curdiv(Ri,Rset,Vin=None,Iin=None,Vout=False):
    """
    curdiv Function
    
    This function is disigned to accept the input current, or input
    voltage to a resistor (or impedance) network of parallel resistors
    (impedances) and calculate the current through a particular element.
    
    Parameters
    ----------
    Ri:     float
            The Particular Resistor of Interest, should not be included in
            the tuple passed to Rset.
    Rset:   float
            Tuple of remaining resistances (impedances) in network.
    Vin:    float, optional
            The input voltage for the system, default=None
    Iin:    float, optional
            The input current for the system, default=None
    Vout:   bool, optional
            Control Argument to enable return of the voltage across the
            resistor (impecance) of interest (Ri)
    
    Returns
    -------
    Opt1 - Ii:          The Current through the resistor (impedance) of interest
    Opt2 - (Ii,Vi):     The afore mentioned current, and voltage across the
                        resistor (impedance) of interest
    """
    # Validate Tuple
    if not isinstance(Rset,tuple):
        Rset = (Rset,) # Set as Tuple
    # Calculate The total impedance
    Rtot = parallelz( Rset + (Ri,) ) # Combine tuples, then calculate total resistance
    # Determine Whether Input was given as Voltage or Current
    if(Vin!=None and Iin==None): # Vin Provided
        Iin = Vin / Rtot # Calculate total current
        Ii = Iin * Rtot/Ri # Calculate the current of interest
    elif(Vin==None and Iin!=None): # Iin provided
        Ii = Iin * Rtot/Ri # Calculate the current of interest
    else:
        raise ValueError("ERROR: Too many or too few constraints provided.")
    if(Vout): # Asked for voltage across resistor of interest
        Vi = Ii * Ri
        return(Ii, Vi)
    else:
        return(Ii)

# Define Instantaneous Power Calculator
def instpower(P,Q,t,f=60):
    """
    instpower Function
    
    This function is designed to calculate the instantaneous power at a
    specified time t given the magnitudes of P and Q.
    
    Parameters
    ----------
    P:  float
        Magnitude of Real Power
    Q:  float
        Magnitude of Reactive Power
    t:  float
        Time at which to evaluate
    f:  float, optional
        System frequency (in Hz), default=60
    
    Returns
    -------
    Pinst:  float
            Instantaneous Power at time t
    """
    # Evaluate omega
    w = 2*np.pi*f
    # Calculate
    Pinst = P + P*np.cos(2*w*t) - Q*np.sin(2*w*t)
    return(Pinst)

# Define Delta-Wye Impedance Network Calculator
def dynetz(delta=None,wye=None,round=None):
    """
    dynetz Function
    
    This function is designed to act as the conversion utility
    to transform delta-connected impedance values to wye-
    connected and vice-versa.
    
    Parameters
    ----------
    delta:  tuple of float, exclusive
            Tuple of the delta-connected impedance values as:
            { Z12, Z23, Z31 }, default=None
    wye:    tuple of float, exclusive
            Tuple of the wye-connected impedance valuse as:
            { Z1, Z2, Z3 }, default=None
    
    Returns
    -------
    delta-set:  tuple of float
                Delta-Connected impedance values { Z12, Z23, Z31 }
    wye-set:    tuple of float
                Wye-Connected impedance values { Z1, Z2, Z3 }
    """
    # Determine which set of impedances was provided
    if(delta!=None and wye==None):
        Z12, Z23, Z31 = delta # Gather particular impedances
        Zsum = Z12 + Z23 + Z31 # Find Sum
        # Calculate Wye Impedances
        Z1 = Z12*Z31 / Zsum
        Z2 = Z12*Z23 / Zsum
        Z3 = Z23*Z31 / Zsum
        Zset = ( Z1, Z2, Z3 )
        if round!=None: Zset = np.around(Zset,round)
        return(Zset) # Return Wye Impedances
    elif(delta==None and wye!=None):
        Z1, Z2, Z3 = wye # Gather particular impedances
        Zmultsum = Z1*Z2 + Z2*Z3 + Z3*Z1
        Z23 = Zmultsum / Z1
        Z31 = Zmultsum / Z2
        Z12 = Zmultsum / Z3
        Zset = ( Z12, Z23, Z31 )
        if round!=None: Zset = np.around(Zset,round)
        return(Zset) # Return Delta Impedances
    else:
        raise ValueError("ERROR: Either delta or wye impedances must be specified.")

# Define Single Line Power Flow Calculator
def powerflow( Vsend, Vrec, Zline ):
    """
    powerflow Function:
    
    This function is designed to calculate the ammount of real
    power transferred from the sending end to the recieving end
    of an electrical line given the sending voltage (complex),
    the receiving voltage (complex) and the line impedance.
    
    Parameters
    ----------
    Vsend:      complex
                The sending-end voltage, should be complex
    Vrec:       complex
                The receiving-end voltage, should be complex
    Zline:      complex
                The line impedance, should be complex
    
    Returns
    -------
    pflow:      complex
                The power transferred from sending-end to
                receiving-end, positive values denote power
                flow from send to receive, negative values
                denote vice-versa.
    """
    # Evaluate the Input Terms
    Vs = abs( Vsend )
    ds = c.phase( Vsend )
    Vr = abs( Vrec )
    dr = c.phase( Vrec )
    # Calculate Power Flow
    pflow = (Vs * Vr)/(Zline) * np.sin( ds-dr )
    return( pflow )
        
# Define Heat-Sink "Resistance" Calculator
def heatsink(P=None,Tjunct=None,Tamb=None,Rjc=None,
             Rcs=None,Rsa=None,Rca=None):
    """
    HEATSINK Function
    
    Purpose:
    --------
    This function is designed to find the missing variable
    from the set of provided variables for a heat-sink
    electro-mechanical analog system.
    
    Tjunct *---[Rjc]----[Rcs]--[Rsa]--------* Tamb
           |          |               |     |
           |          |-----[Rca]-----|     |
          (^)P                             (+)Tamb
           |                                |
           |--------------------------------|
    
    Required Arguments:
    -------------------
    None - This function needs n-1 inputs to solve for nth input
    
    Optional Arguments:
    -------------------
    P:      The Power to be disipated (in watts), default=None
    Tjunct: The temperature at the junction (°C), default=None
    Tamb:   The ambient temperature (°C), default=None
    Rjc:    The thermal-resistance between the junction and the
            case, default=None
    Rcs:    The thermal-resistance between the case and the sink,
            default=None
    Rsa:    The thermal-resistance between the sink and ambient,
            default=None
    Rca:    The thermal-resistance between the case and ambient,
            default=None
            
    Returns:
    --------
    P:      Provided set: { Tjunct, Tamb, Rjc, Rca }    or
                          { Tjunct, Tamb, Rjc, Rca, Rcs, Rsa }
    Tjunct: Provided set: { P, Tamb, Rjc, Rca }         or
                          { P, Tamb, Rjc, Rca, Rcs, Rsa }
    Tamb:   Provided set: { P, Tjunct, Rjc, Rca }       or
                          { P, Tjunct, Rjc, Rca, Rcs, Rsa }
    Rjc:    Provided set: { P, Tjunct, Tamb, Rca }      or
                          { P, Tjunct, Tamb, Rca, Rcs, Rsa }
    Rca:    Provided set: { P, Tjunct, Tamb, Rjc }      or
                          { P, Tjunct, Tamb, Rjc, Rcs, Rsa }
    Rcs:    Provided set: { P, Tjunct, Tamb, Rjc, Rca, Rsa }
    Rsa:    Provided set: { P, Tjunct, Tamb, Rjc, Rca, Rcs }
    """
    # Function needs some serious development
    return(0)
            
    
# Define Impedance From Power and X/R
def zsource(S,V,XoverR,Sbase=None,Vbase=None):
    """
    zsource Function
    
    Used to calculate the source impedance given the apparent power
    magnitude and the X/R ratio.
    
    Parameters
    ----------
    S:          float
                The (rated) apparent power magnitude of the source.
                This may also be refferred to as the "Short-Circuit MVA"
    V:          float
                The (rated) voltage of the source terminals.
    XoverR:     float
                The X/R ratio rated for the source.
    Sbase:      float, optional
                The per-unit base for the apparent power. If set to
                None, will automatically force Sbase to equal S.
                If set to True will treat S as the per-unit value.
    Vbase:      float, optional
                The per-unit base for the terminal voltage. If set to
                None, will automaticlaly force Vbase to equal V. If
                set to True, will treat V as the per-unit value.
    
    Returns
    -------
    Zsource_pu: complex
                The per-unit evaluation of the source impedance.
    """
    # Force Sbase and Vbase if needed
    if Vbase == None:
        Vbase = V
    if Sbase == None:
        Sbase = S
    # Prevent scaling if per-unit already applied
    if Vbase == True:
        Vbase = 1
    if Sbase == True:
        Sbase = 1
    # Set to per-unit
    Spu = S/Sbase
    Vpu = V/Vbase
    # Evaluate Zsource Magnitude
    Zsource_pu = Vpu**2/Spu
    # Evaluate the angle
    nu = np.degrees(np.arctan(XoverR))
    Zsource_pu = eep.phasor(Zsource_pu, nu)
    return(Zsource_pu)

# Define Impedance Decomposer
def zdecompose(Zmag,XoverR):
    """
    zdecompose Function
    
    A function to decompose the impedance magnitude into its
    corresponding resistance and reactance using the X/R ratio.
    
    It is possible to "neglect" R, or make it a very small number;
    this is done by setting the X/R ratio to a very large number
    (X being much larger than R).
    
    Parameters
    ----------
    Zmag:       float
                The magnitude of the impedance.
    XoverR:     float
                The X/R ratio.
    
    Returns
    -------
    R:          float
                The resistance (in ohms)
    X:          float
                The reactance (in ohms)
    """
    # Evaluate Resistance
    R = Zmag/np.sqrt(XoverR**2+1)
    # Evaluate Reactance
    X = R * XoverR
    # Return
    return(R,X)

# Define HP to Watts Calculation
def watts(hp):
    """
    watts Formula
    
    Calculates the power (in watts) given the
    horsepower.
    
    Parameters
    ----------
    hp:         float
                The horspower to compute.
    
    Returns
    -------
    watts:      float
                The power in watts.
    """
    return(hp * 745.699872)

# Define Watts to HP Calculation
def horsepower(watts):
    """
    horsepower Function
    
    Calculates the power (in horsepower) given
    the power in watts.
    
    Parameters
    ----------
    watts:      float
                The wattage to compute.
    
    Returns
    -------
    hp:         float
                The power in horsepower.
    """
    return(watts / 745.699872)
    
# Define Power Reactance Calculator
def powerimpedance(S,V,parallel=False):
    """
    powerimpedance Function
    
    Function to determine the ohmic resistance/reactance
    (impedance) represented by the apparent power (S).
    
    Formula:    Z = V^2 / S           (series components)
                Z = V^2 / (3*S)       (parallel components)
    
    Parameters
    ----------
    S:          complex
                The apparent power of the passive element,
                may be purely resistive or purely reactive.
    V:          float
                The operating voltage of the passive element.
    parallel:   bool, optional
                Control point to specify whether the ohmic
                impedance should be returned as series components
                (False opt.) or parallel components (True opt.).
    
    Returns
    -------
    R:          float
                The ohmic resistance required to consume the
                specified apparent power (S) at the rated
                voltage (V).
    X:          float
                The ohmic reactance required to consume the
                specified apparent power (S) at the rated
                voltage (V).
    """
    # Condition Inputs
    V = abs(V)
    # Test for Parallel Component Option and Evaluate
    if isinstance(S,complex):
        # Complex Power (both R and X)
        if parallel:
            R = V**2 / (3*S.real)
            X = V**2 / (3*S.imag)
        else:
            R = V**2 / (S.real)
            X = V**2 / (S.imag)
        return( R, X )
    # Not Complex (just R)
    R = V**2 / S
    return( R )

# Define Cold-Junction-Voltage Calculator
def coldjunction(Tcj,coupletype="K",To=None,Vo=None,P1=None,P2=None,
                 P3=None,P4=None,Q1=None,Q2=None,round=None):
    """
    coldjunction Function
    
    Function to calculate the expected cold-junction-voltage given
    the temperature at the cold-junction.
    
    Parameters
    ----------
    Tcj:        float
                The temperature (in degrees C) that the junction is
                currently subjected to.
    coupletype: string, optional
                Thermocouple Type, may be one of (B,E,J,K,N,R,S,T), default="K"
    To:         float, optional
                Temperature Constant used in Polynomial.
    Vo:         float, optional
                Voltage Constant used in Polynomial.
    P1:         float, optional
                Polynomial constant.
    P2:         float, optional
                Polynomial constant.
    P3:         float, optional
                Polynomial constant.
    P4:         float, optional
                Polynomial constant.
    Q1:         float, optional
                Polynomial constant.
    Q2:         float, optional
                Polynomial constant.
    Q3:         float, optional
                Polynomial constant.
    round:      int, optional
                Control input to specify how many decimal places the result
                should be rounded to, default=1.
    
    Returns
    -------
    Vcj:        float
                The calculated cold-junction-voltage in volts.
    """
    # Condition Inputs
    coupletype = coupletype.upper()
    # Validate Temperature Range
    if coupletype == "B":
        if not (0 < Tcj and Tcj < 70):
            raise ValueError("Temperature out of range.")
    else:
        if not (-20 < Tcj and Tcj < 70):
            raise ValueError("Temperature out of range.")
    # Define Constant Lookup System
    lookup = ["B","E","J","K","N","R","S","T"]
    if not (coupletype in lookup):
        raise ValueError("Invalid Thermocouple Type")
    index = lookup.index(coupletype)
    # Define Constant Dictionary
    constants = {   "To" : [4.2000000E+01,2.5000000E+01,2.5000000E+01,2.5000000E+01,7.0000000E+00,2.5000000E+01,2.5000000E+01,2.5000000E+01],
                    "Vo" : [3.3933898E-04,1.4950582E+00,1.2773432E+00,1.0003453E+00,1.8210024E-01,1.4067016E-01,1.4269163E-01,9.9198279E-01],
                    "P1" : [2.1196684E-04,6.0958443E-02,5.1744084E-02,4.0514854E-02,2.6228256E-02,5.9330356E-03,5.9829057E-03,4.0716564E-02],
                    "P2" : [3.3801250E-06,-2.7351789E-04,-5.4138663E-05,-3.8789638E-05,-1.5485539E-04,2.7736904E-05,4.5292259E-06,7.1170297E-04],
                    "P3" : [-1.4793289E-07,-1.9130146E-05,-2.2895769E-06,-2.8608478E-06,2.1366031E-06,-1.0819644E-06,-1.3380281E-06,6.8782631E-07],
                    "P4" : [-3.3571424E-09,-1.3948840E-08,-7.7947143E-10,-9.5367041E-10,9.2047105E-10,-2.3098349E-09,-2.3742577E-09,4.3295061E-11],
                    "Q1" : [-1.0920410E-02,-5.2382378E-03,-1.5173342E-03,-1.3948675E-03,-6.4070932E-03,2.6146871E-03,-1.0650446E-03,1.6458102E-02],
                    "Q2" : [-4.9782932E-04,-3.0970168E-04,-4.2314514E-05,-6.7976627E-05,8.2161781E-05,-1.8621487E-04,-2.2042420E-04,0.0000000E+00]
                }
    # Load Data Into Terms
    if To == None:
        To = constants["To"][index]
    if Vo == None:
        Vo = constants["Vo"][index]
    if P1 == None:
        P1 = constants["P1"][index]
    if P2 == None:
        P2 = constants["P2"][index]
    if P3 == None:
        P3 = constants["P3"][index]
    if P4 == None:
        P4 = constants["P4"][index]
    if Q1 == None:
        Q1 = constants["Q1"][index]
    if Q2 == None:
        Q2 = constants["Q2"][index]
    # Define Formula Terms
    tx = (Tcj-To)
    num = tx*(P1+tx*(P2+tx*(P3+P4*tx)))
    den = 1+tx*(Q1+Q2*tx)
    Vcj = Vo + num/den
    # Round Value if Allowed
    if round != None:
        Vcj = np.around(Vcj, round)
    # Return in milivolts
    return(Vcj*m)
    
# Define Thermocouple Temperature Calculation
def thermocouple(V,coupletype="K",fahrenheit=False,cjt=None,To=None,Vo=None,P1=None,
                 P2=None,P3=None,P4=None,Q1=None,Q2=None,Q3=None,round=1):
    """
    thermocouple Function
    
    Utilizes polynomial formula to calculate the temperature being monitored
    by a thermocouple. Allows for various thermocouple types (B,E,J,K,N,R,S,T)
    and various cold-junction-temperatures.
    
    Parameters
    ----------
    V:          float
                Measured voltage (in Volts)
    coupletype: string, optional
                Thermocouple Type, may be one of (B,E,J,K,N,R,S,T), default="K"
    fahrenheit: bool, optional
                Control to enable return value as Fahrenheit instead of Celsius,
                default=False
    cjt:        float, optional
                Cold-Junction-Temperature
    To:         float, optional
                Temperature Constant used in Polynomial.
    Vo:         float, optional
                Voltage Constant used in Polynomial.
    P1:         float, optional
                Polynomial constant.
    P2:         float, optional
                Polynomial constant.
    P3:         float, optional
                Polynomial constant.
    P4:         float, optional
                Polynomial constant.
    Q1:         float, optional
                Polynomial constant.
    Q2:         float, optional
                Polynomial constant.
    Q3:         float, optional
                Polynomial constant.
    round:      int, optional
                Control input to specify how many decimal places the result
                should be rounded to, default=1.
    
    Returns
    -------
    T:          float
                The temperature (by default in degrees C, but optionally in
                degrees F) as computed by the function.
    """
    # Condition Inputs
    coupletype = coupletype.upper()
    V = V/m # Scale volts to milivolts
    # Determine Cold-Junction-Voltage
    if cjt != None:
        Vcj = coldjunction(cjt,coupletype,To,Vo,P1,P2,P3,P4,Q1,Q2,round)
        V += Vcj/m
    # Define Constant Lookup System
    lookup = ["B","E","J","K","N","R","S","T"]
    if not (coupletype in lookup):
        raise ValueError("Invalid Thermocouple Type")
    # Define Voltage Ranges
    voltages = {"J" : [-8.095, 0,      21.840, 45.494, 57.953, 69.553],
                "K" : [-6.404, -3.554, 4.096,  16.397, 33.275, 69.553],
                "B" : [0.291,  2.431,  13.820, None,   None,   None],
                "E" : [-9.835, -5.237, 0.591,  24.964, 53.112, 76.373],
                "N" : [-4.313, 0,      20.613, 47.513, None,   None],
                "R" : [-0.226, 1.469,  7.461,  14.277, 21.101, None],
                "S" : [-0.236, 1.441,  6.913,  12.856, 18.693, None],
                "T" : [-6.18,  -4.648, 0,      9.288,  20.872, None]}
    # Determine Array Selection
    vset = voltages[coupletype]
    if V < vset[0]*m:
        raise ValueError("Voltage Below Lower Bound")
    elif vset[0] <= V < vset[1]:
        select = 0
    elif vset[1] <= V < vset[2]:
        select = 1
    elif vset[2] <= V < vset[3]:
        select = 2
    elif vset[3] <= V < vset[4]:
        select = 3
    elif vset[4] <= V <= vset[5]:
        select = 4
    elif vset[5] < V:
        raise ValueError("Voltage Above Upper Bound")
    else:
        raise ValueError("Internal Error!")
    # Define Dictionary of Arrays
    data = {"J" : [[-6.4936529E+01,2.5066947E+02,6.4950262E+02,9.2510550E+02,1.0511294E+03],
                   [-3.1169773E+00,1.3592329E+01,3.6040848E+01,5.3433832E+01,6.0956091E+01],
                   [2.2133797E+01,1.8014787E+01,1.6593395E+01,1.6243326E+01,1.7156001E+01],
                   [2.0476437E+00,-6.5218881E-02,7.3009590E-01,9.2793267E-01,-2.5931041E+00],
                   [-4.6867532E-01,-1.2179108E-02,2.4157343E-02,6.4644193E-03,-5.8339803E-02],
                   [-3.6673992E-02,2.0061707E-04,1.2787077E-03,2.0464414E-03,1.9954137E-02],
                   [1.1746348E-01,-3.9494552E-03,4.9172861E-02,5.2541788E-02,-1.5305581E-01],
                   [-2.0903413E-02,-7.3728206E-04,1.6813810E-03,1.3682959E-04,-2.9523967E-03],
                   [-2.1823704E-03,1.6679731E-05,7.6067922E-05,1.3454746E-04,1.1340164E-03]],
            "K" : [[-1.2147164E+02,-8.7935962E+00,3.1018976E+02,6.0572562E+02,1.0184705E+03],
                   [-4.1790858E+00,-3.4489914E-01,1.2631386E+01,2.5148718E+01,4.1993851E+01],
                   [3.6069513E+01,2.5678719E+01,2.4061949E+01,2.3539401E+01,2.5783239E+01],
                   [3.0722076E+01,-4.9887904E-01,4.0158622E+00,4.6547228E-02,-1.8363403E+00],
                   [7.7913860E+00,-4.4705222E-01,2.6853917E-01,1.3444400E-02,5.6176662E-02],
                   [5.2593991E-01,-4.4869203E-02,-9.7188544E-03,5.9236853E-04,1.8532400E-04],
                   [9.3939547E-01,2.3893439E-04,1.6995872E-01,8.3445513E-04,-7.4803355E-02],
                   [2.7791285E-01,-2.0397750E-02,1.1413069E-02,4.6121445E-04,2.3841860E-03],
                   [2.5163349E-02,-1.8424107E-03,-3.9275155E-04,2.5488122E-05,0.0]],
            "B" : [[5.0000000E+02,1.2461474E+03],
                   [1.2417900E+00,7.2701221E+00],
                   [1.9858097E+02,9.4321033E+01],
                   [2.4284248E+01,7.3899296E+00],
                   [-9.7271640E+01,-1.5880987E-01],
                   [-1.5701178E+01,1.2681877E-02],
                   [3.1009445E-01,1.0113834E-01],
                   [-5.0880251E-01,-1.6145962E-03],
                   [-1.6163342E-01,-4.1086314E-06]],
            "E" : [[-1.1721668E+02,-5.0000000E+01,2.5014600E+02,6.0139890E+02,8.0435911E+02],
                   [-5.9901698E+00,-2.7871777E+00,1.7191713E+01,4.5206167E+01,6.1359178E+01],
                   [2.3647275E+01,1.9022736E+01,1.3115522E+01,1.2399357E+01,1.2759508E+01],
                   [1.2807377E+01,-1.7042725E+00,1.1780364E+00,4.3399963E-01,-1.1116072E+00],
                   [2.0665069E+00,-3.5195189E-01,3.6422433E-02,9.1967085E-03,3.5332536E-02],
                   [8.6513472E-02,4.7766102E-03,3.9584261E-04,1.6901585E-04,3.3080380E-05],
                   [5.8995860E-01,-6.5379760E-02,9.3112756E-02,3.4424680E-02,-8.8196889E-02],
                   [1.0960713E-01,-2.1732833E-02,2.9804232E-03,6.9741215E-04,2.8497415E-03],
                   [6.1769588E-03,0.0,3.3263032E-05,1.2946992E-05,0.0]],
            "N" : [[-5.9610511E+01,3.1534505E+02,1.0340172E+03],
                   [-1.5000000E+00,9.8870997E+00,3.7565475E+01],
                   [4.2021322E+01,2.7988676E+01,2.6029492E+01],
                   [4.7244037E+00,1.5417343E+00,-6.0783095E-01],
                   [-6.1153213E+00,-1.4689457E-01,-9.7742562E-03],
                   [-9.9980337E-01,-6.8322712E-03,-3.3148813E-06],
                   [1.6385664E-01,6.2600036E-02,-2.5351881E-02],
                   [-1.4994026E-01,-5.1489572E-03,-3.8746827E-04],
                   [-3.0810372E-02,-2.8835863E-04,1.7088177E-06]],
            "R" : [[1.3054315E+02,5.4188181E+02,1.0382132E+03,1.5676133E+03],
                   [8.8333090E-01,4.9312886E+00,1.1014763E+01,1.8397910E+01],
                   [1.2557377E+02,9.0208190E+01,7.4669343E+01,7.1646299E+01],
                   [1.3900275E+02,6.1762254E+00,3.4090711E+00,-1.0866763E+00],
                   [3.3035469E+01,-1.2279323E+00,-1.4511205E-01,-2.0968371E+00],
                   [-8.5195924E-01,1.4873153E-02,6.3077387E-03,-7.6741168E-01],
                   [1.2232896E+00,8.7670455E-02,5.6880253E-02,-1.9712341E-02],
                   [3.5603023E-01,-1.2906694E-02,-2.0512736E-03,-2.9903595E-02],
                   [0.0,0.0,0.0,-1.0766878E-02]],
            "S" : [[1.3792630E+02,4.7673468E+02,9.7946589E+02,1.6010461E+03],
                   [9.3395024E-01,4.0037367E+00,9.3508283E+00,1.6789315E+01],
                   [1.2761836E+02,1.0174512E+02,8.7126730E+01,8.4315871E+01],
                   [1.1089050E+02,-8.9306371E+00,-2.3139202E+00,-1.0185043E+01],
                   [1.9898457E+01,-4.2942435E+00,-3.2682118E-02,-4.6283954E+00],
                   [9.6152996E-02,2.0453847E-01,4.6090022E-03,-1.0158749E+00],
                   [9.6545918E-01,-7.1227776E-02,-1.4299790E-02,-1.2877783E-01],
                   [2.0813850E-01,-4.4618306E-02,-1.2289882E-03,-5.5802216E-02],
                   [0.0,1.6822887E-03,0.0,-1.2146518E-02]],
            "T" : [[-1.9243000E+02,-6.0000000E+01,1.3500000E+02,3.0000000E+02],
                   [-5.4798963E+00,-2.1528350E+00,5.9588600E+00,1.4861780E+01],
                   [5.9572141E+01,3.0449332E+01,2.0325591E+01,1.7214707E+01],
                   [1.9675733E+00,-1.2946560E+00,3.3013079E+00,-9.3862713E-01],
                   [-7.8176011E+01,-3.0500735E+00,1.2638462E-01,-7.3509066E-02],
                   [-1.0963280E+01,-1.9226856E-01,-8.2883695E-04,2.9576140E-04],
                   [2.7498092E-01,6.9877863E-03,1.7595577E-01,-4.8095795E-02],
                   [-1.3768944E+00,-1.0596207E-01,7.9740521E-03,-4.7352054E-03],
                   [-4.5209805E-01,-1.0774995E-02,0.0,0.0]]}
    # Load Data Into Terms
    if To == None:
        To = data[coupletype][0][select]
    if Vo == None:
        Vo = data[coupletype][1][select]
    if P1 == None:
        P1 = data[coupletype][2][select]
    if P2 == None:
        P2 = data[coupletype][3][select]
    if P3 == None:
        P3 = data[coupletype][4][select]
    if P4 == None:
        P4 = data[coupletype][5][select]
    if Q1 == None:
        Q1 = data[coupletype][6][select]
    if Q2 == None:
        Q2 = data[coupletype][7][select]
    if Q3 == None:
        Q3 = data[coupletype][8][select]
    # Calculate Temperature in Degrees C
    num = (V-Vo)*(P1+(V-Vo)*(P2+(V-Vo)*(P3+P4*(V-Vo))))
    den = 1 + (V-Vo)*(Q1+(V-Vo)*(Q2+Q3*(V-Vo)))
    temp = To + num/den
    # Return Temperature
    if fahrenheit:
        temp = (temp*9/5)+32
    temp = np.around(temp,round)
    return(temp)

# END OF FILE