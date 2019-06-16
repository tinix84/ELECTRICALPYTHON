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
tuple = "<class 'tuple'>"
ndarr = "<class 'numpy.ndarray'>"
tint = "<class 'int'>"
tfloat = "<class 'float'>"
tfun = "<class 'function'>"
tstr = "<class 'str'>"

# Define Phasor Generator
def phasor( mag, ang ):
    """
    PHASOR Function:
    
    Purpose:
    --------
    Generates the standard Pythonic complex representation
    of a phasor voltage or current when given the magnitude
    and angle of the specific voltage or current.
    
    Required Arguments:
    -------------------
    mag:        The Magnitude of the Voltage/Current
    ang:        The Angle (in degrees) of the Voltage/Current
    
    Returns:
    --------
    phasor:     Standard Pythonic Complex Representation of
                the specified voltage or current.
    """
    return( c.rect( mag, np.radians( ang ) ) )

# Define Reactance Calculator
def reactance(z,f=60,sensetivity=1e-12):
    """
    REACTANCE Function:
    
    Purpose:
    --------
    Calculates the Capacitance or Inductance in Farads or Henreys
    (respectively) provided the impedance of an element.
    Will return capacitance (in Farads) if ohmic impedance is
    negative, or inductance (in Henrys) if ohmic impedance is
    positive. If imaginary: calculate with j factor (imaginary number).
    
    Required Arguments:
    -------------------
    z:      The Impedance Provided, may be complex (R+jI)
    
    Optional Arguments:
    -------------------
    f:              The Frequency Base for Provided Impedance, default=60
    sensetivity:    The sensetivity used to check if a resistance was
                    provided, default=1e-12
    
    Returns:
    --------
    out:    Capacitance or Inductance of Impedance
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
    CPRINT Function
    
    Purpose:
    --------
    This function is designed to accept a complex value (val) and print
    the value in the standard electrical engineering notation:
    
        magnitude ∠ angle °
    
    This function will print the magnitude in degrees, and can print
    a unit and label in addition to the value itself.
    
    Required Arguments:
    -------------------
    val:        The Complex Value to be Printed, may be singular value,
                tuple of values, or list/array.
    
    Optional Arguments:
    -------------------
    unit:       The string to be printed corresponding to the unit mark.
                default=""
    label:      The pre-pended string used as a descriptive labeling string.
                default=""
    printval:   Control argument enabling/disabling printing of the string.
                default=True
    ret:        Control argument allowing the evaluated value to be returned.
                default=False
    round:      Control argument specifying how many decimals of the complex
                value to be printed. May be negative to round to spaces
                to the left of the decimal place (follows standard round()
                functionality). default=3
    
    Returns:
    --------
    numarr:     The array of values corresponding to the magnitude and angle,
                values are returned in the form:
                
                    [[ mag, ang ],
                     [ mag, ang ],
                          ...    ,
                     [ mag, ang ]]
                
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
    PHASORZ Function:
    
    Purpose:
    --------
    This function's purpose is to generate the phasor-based
    impedance of the specified input given as either the
    capacitance (in Farads) or the inductance (in Henreys).
    The function will return the phasor value (in Ohms).
    
    Required Arguments:
    -------------------
    Either C or L must be specified.
    
    Optional Arguments:
    -------------------
    C:       The capacitance value (specified in Farads),
             default=None
    L:       The inductance value (specified in Henreys),
             default=None
    f:       The system frequency to be calculated upon, default=60
    complex: Control argument to specify whether the returned
             value should be returned as a complex value.
             default=True
    
    Returns:
    --------
    Z:      The ohmic impedance of either C or L (respectively).
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
    PARALLELZ Function:
    
    Purpose:
    --------
    This function is designed to generate the total parallel
    impedance of a set (tuple) of impedances specified as real
    or complex values.
    
    Required Arguments:
    -------------------
    Z:      The tupled input set of impedances, may be a tuple
            of any size greater than 2. May be real, complex, or
            a combination of the two.
    
    Returns:
    --------
    Zp:     The calculated parallel impedance of the input tuple.
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
    PHASELINE Function
    
    Purpose:
    --------
    This function is designed to return the phase- or line-equivalent
    of the voltage/current provided. It is designed to be used when
    converting delta- to wye-connections and vice-versa.
    Given a voltage of one type, this function will return the
    voltage of the opposite type. The same is true for current.
    
    Required Arguments:
    -------------------
    NONE - All optional, one optional argument MUST be specified.
    
    Optional Arguments:
    -------------------
    VLL:        The Line-to-Line Voltage; default=None
    VLN:        The Line-to-Neutral Voltage; default=None
    Iline:      The Line-Current; default=None
    Iphase:     The Phase-Current; default=None
    complex:    Control to return value in complex form; default=False
    
    Returns:
    --------
    out:        The opposite type of the input; i.e.:
                    GIVEN       RETURNED
                    VLL         VLN
                    VLN         VLL
                    Iline       Iphase
                    Iphase      Iline
                Output may be returned as complex if desired.
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
    POWERSET Function
    
    Purpose:
    --------
    This function is designed to calculate all values
    in the set { P, Q, S, PF } when two (2) of the
    values are provided. The equations in this
    function are prepared for AC values, that is:
    real and reactive power, apparent power, and power
    factor.
    
    Required Arguments:
    -------------------
    NONE; a minimum of two of the optional arguments
          must be entered for proper execution.
    
    Optional Arguments:
    -------------------
    P:      Real Power, unitless; default=None
    Q:      Reactive Power, unitless; default=None
    S:      Apparent Power, unitless; default=None
    PF:     Power Factor, unitless, provided as a
            decimal value, lagging is positive,
            leading is negative; default=None
    
    Returns:
    --------
    ( P, Q, S, PF ):    Completely calculated set,
                        all terms are as described
                        above.
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
    POWERTRIANGLE Function
    
    Purpose:
    --------
    This function is designed to draw a power triangle given
    values for the complex power system.
    
    Required Arguments:
    -------------------
    NONE; a minimum of two of the optional arguments
          must be entered for proper execution.
    
    Optional Arguments:
    -------------------
    P:          Real Power, unitless; default=None
    Q:          Reactive Power, unitless; default=None
    S:          Apparent Power, unitless; default=None
    PF:         Power Factor, unitless, provided as a
                decimal value, lagging is positive,
                leading is negative; default=None
    color:      The color of the power triangle lines;
                default="red"
    text:       The title of the power triangle plot,
                default="Power Triangle"
    printval:   Control argument to allow the numeric
                values to be printed on the plot,
                default="False"
    
    Returns:
    --------
    NONE;   plots generated.
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
#
#   Calculates Req and Xeq, or Rc and Xm, or both sets given three
#   values from a specific set of inputs { Poc, Voc, Ioc,  Psc,
#   Vsc, Isc }.
#
#   Requires one or both of two sets: { Poc, Voc, Ioc }, or
#   { Psc, Vsc, Isc }.
#   Will return: { Rc, Xm } given first set, { Req, Xeq } given
#   second set.
#   All values given must be given as absolute value, not complex.
#   All values returned are given with respect to high-side/primary
###################################################################
def transformertest(Poc=False,Voc=False,Ioc=False,Psc=False,Vsc=False,
               Isc=False):
    """
    TRANSFORMERTEST Function
    
    Purpose:
    --------
    This function will determine the non-ideal circuit components of
    a transformer (Req and Xeq, or Rc and Xm) given the test-case
    parameters for the open-circuit test and/or the closed-circuit
    test. Requires one or both of two sets: { Poc, Voc, Ioc }, or
    { Psc, Vsc, Isc }.
    All values given must be given as absolute value, not complex.
    All values returned are given with respect to primary.
    
    Required Arguments:
    -------------------
    NONE,   A minimum of one complete set of optional arguments must
            be provided for function to complete successfully.
            Optional Arg. Sets are: { Poc, Voc, Ioc }, or
            { Psc, Vsc, Isc }.
    
    Optional Arguments:
    -------------------
    Poc:    The open-circuit measured power (real power), default=None
    Voc:    The open-circuit measured voltage (measured on X),
            default=None
    Ioc:    The open-circuit measured current (measured on primary),
            default=None
    Psc:    The short-circuit measured power (real power), default=None
    Vsc:    The short-circuit measured voltage (measured on X),
            default=None
    Isc:    The short-circuit measured current (measured on X),
            default=None
    
    Returns:
    --------
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
def phasorplot(phasor,title="Phasor Diagram",legend=False,bg="#d5de9c",radius=1.2):
    """
    PHASORPLOT Function
    
    Purpose:
    --------
    This function is designed to plot a phasor-diagram with angles in degrees
    for up to 12 phasor sets. Phasors must be passed as a complex number set,
    (e.g. [ m+ja, m+ja, m+ja, ... , m+ja ] ).
    
    Required Arguments:
    -------------------
    phasor:     The set of phasors to be plotted.
    
    Optional Arguments:
    -------------------
    title:      The Plot Title, default="Phasor Diagram"
    legend:     Control argument to enable displaying the legend, must be passed
                as an array or list of strings, default=False
    bg:         Background-Color control, default="#d5de9c"
    radius:     The diagram radius, default=1.2
    """
    numphs = len(phasor)
    
    colors = ["#FF0000","#800000","#FFFF00","#808000","#00ff00","#008000",
            "#00ffff","#008080","#0000ff","#000080","#ff00ff","#800080"]
    
    if numphs > 12:
        raise ValueError("ERROR: No more than 12 phasors allowed.")
    
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

###################################################################
#   Define Non-Linear Power Factor Calculator
#
#   Returns the the unknown variable of the set { PFtrue, PFdist,
#   PFdisp }.
#
#   Requires any two of the three inputs
###################################################################
def nlinpf(PFtrue=False,PFdist=False,PFdisp=False):
    """
    NLINPF Function
    
    Purpose:
    --------
    This function is designed to evaluate one of three unknowns
    given the other two. These particular unknowns are the arguments
    and as such, they are described in the representative sections
    below.
    
    Required Arguments:
    -------------------
    None.
    
    Optional Arguments:
    -------------------
    PFtrue:     The "True" power-factor, default=None
    PFdist:     The "Distorted" power-factor, default=None
    PFdisp:     The "Displacement" power-factor, default=None
    
    Returns:
    --------
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
    ISCRL Function
    
    Purpose:
    --------
    The Isc-RL function (Short Circuit Current for RL Circuit)
    is designed to calculate the short-circuit current for an
    RL circuit.
    
    Required Arguments:
    -------------------
    V:          The absolute magnitude of the voltage.
    Z:          The complex value of the impedance. (R + jX)
    
    Optional Arguments:
    -------------------
    t:          The time at which the value should be calculated,
                should be specified in seconds, default=None
    f:          The system frequency, specified in Hz, default=None
    mxcurrent:  Control variable to enable calculating the value at
                maximum current, default=True
    alpha:      Angle specification, default=None
    
    Returns:
    --------
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
    VOLTDIV Function
    
    Purpose:
    --------
    This function is designed to calculate the output
    voltage of a voltage divider given the input voltage,
    the resistances (or impedances) and the load resistance
    (or impedance) if present.
    
    Required Arguments:
    -------------------
    Vin:    The Input Voltage, may be real or complex
    R1:     The top resistor of the divider (real or complex)
    R2:     The bottom resistor of the divider, the one which
            the output voltage is measured across, may be
            either real or complex
    
    Optional Arguments:
    -------------------
    Rload:  The Load Resistor (or impedance), default=None
    
    Returns:
    --------
    Vout:   The Output voltage as measured across R2 and/or Rload
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
    CURDIV Function
    
    Purpose:
    --------
    This function is disigned to accept the input current, or input
    voltage to a resistor (or impedance) network of parallel resistors
    (impedances) and calculate the current through a particular element.
    
    Required Arguments:
    -------------------
    Ri:     The Particular Resistor of Interest, should not be included in
            the tuple passed to Rset.
    Rset:   Tuple of remaining resistances (impedances) in network.
    
    Optional Arguments:
    -------------------
    Vin:    The input voltage for the system, default=None
    Iin:    The input current for the system, default=None
    Vout:   Control Argument to enable return of the voltage across the
            resistor (impecance) of interest (Ri)
    
    Returns:
    --------
    Opt1 - Ii:          The Current through the resistor (impedance) of interest
    Opt2 - (Ii,Vi):     The afore mentioned current, and voltage across the
                        resistor (impedance) of interest
    """
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
    INSTPOWER Function
    
    Purpose:
    --------
    This function is designed to calculate the instantaneous power at a
    specified time t given the magnitudes of P and Q.
    
    Required Arguments:
    -------------------
    P:  Magnitude of Real Power
    Q:  Magnitude of Reactive Power
    t:  Time at which to evaluate
    
    Optional Arguments:
    -------------------
    f:  System frequency (in Hz), default=60
    
    Returns:
    --------
    Pinst:  Instantaneous Power at time t
    """
    # Evaluate omega
    w = 2*np.pi*f
    # Calculate
    Pinst = P + P*np.cos(2*w*t) - Q*np.sin(2*w*t)
    return(Pinst)

# Define Delta-Wye Impedance Network Calculator
def dynetz(delta=None,wye=None,round=None):
    """
    DYNETZ Function
    
    Purpose:
    --------
    This function is designed to act as the conversion utility
    to transform delta-connected impedance values to wye-
    connected and vice-versa.
    
    Required Arguments:
    -------------------
    *Requires Either delta or wye*
    
    Optional Arguments:
    -------------------
    delta:  Tuple of the delta-connected impedance values as:
            { Z12, Z23, Z31 }, default=None
    wye:    Tuple of the wye-connected impedance valuse as:
            { Z1, Z2, Z3 }, default=None
    
    Returns:
    delta-set: Delta-Connected impedance values { Z12, Z23, Z31 }
    wye-set:   Wye-Connected impedance values { Z1, Z2, Z3 }
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
    POWERFLOW Function:
    
    Purpose:
    --------
    This function is designed to calculate the ammount of real
    power transferred from the sending end to the recieving end
    of an electrical line given the sending voltage (complex),
    the receiving voltage (complex) and the line impedance.
    
    Required Arguments:
    -------------------
    Vsend:      The sending-end voltage, should be complex
    Vrec:       The receiving-end voltage, should be complex
    Zline:      The line impedance, should be complex
    
    Optional Arguments:
    -------------------
    None.
    
    Returns:
    --------
    pflow:      The power transferred from sending-end to
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
            
    
    

# END OF FILE