###################################################################
#   EEPOWER.PY
#
#   A library of functions, constants and more
#   that are related to Power in Electrical Engineering.
#
#   February 13, 2019
#   September 3, 2018
#   August 30, 2018
#   February 20, 2019
#
#   Written by Joe Stanley
#
#   Special Thanks To:
#   Paul Ortmann - Idaho Power
#
#   Included Constants:
#   - 'A' Operator for Symmetrical Components: a
#   - Not a Number value (NaN): NAN
#
#   Symmetrical Components Matricies:
#   - ABC to 012 Conversion:        abc012
#   - 012 to ABC Conversion:        i012abc
#
#   Included Functions
#   - Phasor V/I Generator:         phasor
#   - Phasor Impedance Generator:   phasorz
#   - Complex Display Function:     cprint
#   - Parallel Impedance Adder:     parallelz
#   - V/I Line/Phase Converter:     phaseline
#   - Power Set Values:             powerset
#   - Power Triangle Function:      powertriangle
#   - Transformer SC OC Tests:      trans_scoc
#   - Phasor Plot Generator:        phasorplot
#   - Total Harmonic Distortion:    thd
#   - Total Demand Distortion:      tdd
#   - Reactance Calculator:         reactance
#   - Non-Linear PF Calc:           pf_nonlin
#   - Harmonic Limit Calculator:    harmoniclimit
#
#   Additional functions available in sub-modules:
#   - capacitor.py
#   - perunit.py
#   - systemsolution.py
###################################################################
name = "eepower"
ver = "1.6.11"

# Import Submodules
from .capacitor import *
from .perunit import *
from .systemsolution import *

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
abc012 = 1/3 * np.array([[ 1, 1, 1    ],
                         [ 1, a, a**2 ],
                         [ 1, a**2, a ]])
i012abc = np.array([[ 1, 1, 1    ],
                    [ 1, a**2, a ],
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
def reactance(z,f):
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
    z:      The Impedance Provided
    f:      The Frequency Base for Provided Impedance
    
    Returns:
    --------
    out:    Capacitance or Inductance of Impedance
    """
    w = 2*np.pi*f
    if isinstance(z, complex):
        if (z.imag > 0):
            out = z/(w*1j)
        else:
            out = 1/(w*1j*z)
        out = abs(out)
    else:
        if (z > 0):
            out = z/(w)
        else:
            out = 1/(w*z)
        out = abs(out)
    return(out)

# Define display function
def cprint(val,unit="",label="",printval=True,ret=False,decimals=3):
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
    val:        The Complex Value to be Printed
    
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
    decimals:   Control argument specifying how many decimals of the complex
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
        row, col = val.shape
        sz = val.size
        mult = True
        if label=="":
            label = np.array([])
            for _ in range(sz):
                label = np.append(label,[""])
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
        mag = round( mag, decimals ) #Round
        ang = round( ang, decimals ) #Round
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
def phasorz(f,C=None,L=None,complex=True):
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
    f:      The system frequency to be calculated upon
    
    Optional Arguments:
    -------------------
    C:       The capacitance value (specified in Farads),
             default=None
    L:       The inductance value (specified in Henreys),
             default=None
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
def parallelz(Z):
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
    L = len(Z)
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

    #Plot
    if plot:
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

###################################################################
#   Define Transformer Short-Circuit/Open-Circuit Function
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
def trans_scoc(Poc=False,Voc=False,Ioc=False,Psc=False,Vsc=False,
               Isc=False):
    SC = False
    OC = False
    # Given Open-Circuit Values
    if (Poc!=False) and (Voc!=False) and (Ioc!=False):
        PF = Poc/(Voc*Ioc)
        Y = c.rect(Ioc/Voc,-np.arccos(PF))
        Rc = 1/Y.real
        Xm = -1/Y.imag
        OC = True
    # Given Short-Circuit Values
    if (Psc!=False) and (Vsc!=False) and (Isc!=False):
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
def pf_dist(I1=False,IRMS=False,Ih=False):
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
def pf_nonlin(PFtrue=False,PFdist=False,PFdisp=False):
    if(PFtrue!=False and PFdist!=False and PFdisp!=False):
        raise ValueError("ERROR: Too many constraints, no solution.") 
    elif ( PFdist!=False and PFdisp!=False ):
        return( PFdist * PFdisp )
    elif ( PFtrue!=False and PFdisp!=False ):
        return( PFtrue / PFdisp )
    elif ( PFtrue!=False and PFdist!=False ):
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
    