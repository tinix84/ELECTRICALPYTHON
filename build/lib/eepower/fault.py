###################################################################
#   FAULT.PY
#
#   A library of functions, constants and more that are 
#   related to Power-System-Faults in Electrical Engineering.
#
#   March 18, 2019
#
#   Written by Joe Stanley
#
#   Special Thanks To:
#   Dr. Yacine Chakhchoukh - Faculty University of Idaho
#
#   Included Functions:
#   - Single Line to Ground             phs1g
#   - Double Line to Ground             phs2g
#   - Line to Line                      phs2
#   - Three-Phase Fault                 phs3
#   - Faulted Bus Voltage               busvolt
#   - CT Saturation Function            ct_saturation
#   - CT C-Class Calculator             ct_cclass
#   - CT Sat. V at rated Burden         ct_satratburden
#   - CT Voltage Peak Formula           ct_vpeak
#   - CT Time to Saturation             ct_timetosat
#   - Transient Recovery Voltage Calc.  pktransrecvolt
#   - TRV Reduction Resistor            trvresistor
#   - Natural Frequency Calculator      natfreq
#   - TOC Trip Time                     toctriptime
#   - TOC Reset Time                    tocreset
#   - Pickup Setting Assistant          pickup
#   - Radial TOC Coordination Tool      tdcoordradial
####################################################################

# Import Necessary Libraries
import numpy as np
from scipy.optimize import fsolve

# Import Local Dependencies
from .__init__ import Aabc, A012

# Define Single Line to Ground Fault Function
def phs1g(Vsrc,Xseq,Rf=0,load=None,sequence=True):
    """
    PHS1G Function
    
    Purpose:
    --------
    This function will evaluate the 0, Positive, and Negative
    sequence currents for a single-line-to-ground fault with the
    option of calculating with or without a load.
    
    Required Arguments:
    -------------------
    Vsrc:       The Source Voltage
    Xseq:       Tupple of sequence reactances as (X0, X1, X2)
    
    Optional Arguments:
    -------------------
    Rf:         The fault resistance, default=0
    load:       The load conditions, default=None
    sequence:   Control argument to force return into ABC-Domain Currents
    
    Returns:
    --------
    Ifault: The Fault Current, equal for 0, pos., and neg. seq.
    """
    # Decompose Reactance Tuple
    X0, X1, X2 = Xseq
    # Load Condition
    if(load!=None):
        print("nothing yet")
    else:
        # Ensure that X-components are imaginary
        if(not isinstance(X0, complex)): X0 *= 1j
        if(not isinstance(X1, complex)): X1 *= 1j
        if(not isinstance(X2, complex)): X2 *= 1j
        # Calculate Fault Current
        Ifault = Vsrc / (X0 + X1 + X2 + 3*Rf)
        Ifault = np.array([ Ifault, Ifault, Ifault ])
    # Prepare Value for return
    if not sequence:
        Ifault = A012.dot( Ifault ) # Convert to ABC-Domain
    # Return Value
    return(Ifault)
    
# Define Double Line to Ground Fault Current Calculator
def phs2g(Vsrc,Xseq,Rf=0,load=None,sequence=True):
    """
    PHS2G Function
    
    Purpose:
    --------
    This function will evaluate the 0, Positive, and Negative
    sequence currents for a double-line-to-ground fault with the
    option of calculating with or without a load.
    
    Required Arguments:
    -------------------
    Vsrc:       The Source Voltage
    Xseq:       Tupple of sequence reactances as (X0, X1, X2)
    
    Optional Arguments:
    -------------------
    Rf:         The fault resistance, default=0
    load:       The load conditions, default=None
    sequence:   Control argument to force return into ABC-Domain Currents
    
    Returns:
    --------
    Ifault: The Array of Fault Currents as (If0, If1, If2)
    """
    # Decompose Reactance Tuple
    X0, X1, X2 = Xseq
    # Load Condition
    if(load!=None):
        print("nothing yet")
    else:
        # Ensure that X-components are imaginary
        if(not isinstance(X0, complex)): X0 *= 1j
        if(not isinstance(X1, complex)): X1 *= 1j
        if(not isinstance(X2, complex)): X2 *= 1j
        # Calculate Fault Currents
        If1 = Vsrc / (X1 + (X2*(X0+3*Rf))/(X0+X2+3*Rf))
        If2 = -(Vsrc - X1*If1)/X2
        If0 = -(Vsrc - X1*If1)/(X0+3*Rf)
        faults = np.array([If0, If1, If2])
    # Return Currents
    if not sequence:
        faults = A012.dot(faults.T)
    return(faults)

# Define Phase-to-Phase Fault Current Calculator
def phs2(Vsrc,Xseq,Rf=0,load=None,sequence=True):
    """
    PHS2 Function
    
    Purpose:
    --------
    This function will evaluate the 0, Positive, and Negative
    sequence currents for a phase-to-phase fault with the
    option of calculating with or without a load.
    
    Required Arguments:
    -------------------
    Vsrc:       The Source Voltage
    Xseq:       Tupple of sequence reactances as (X0, X1, X2)
    
    Optional Arguments:
    -------------------
    Rf:         The fault resistance, default=0
    load:       The load conditions, default=None
    sequence:   Control argument to force return into ABC-Domain Currents
    
    Returns:
    --------
    Ifault: The Array of Fault Currents as (If0, If1, If2)
    """
    # Decompose Reactance Tuple
    X0, X1, X2 = Xseq
    # Load Condition
    if(load!=None):
        print("nothing yet")
    else:
        # Ensure that X-components are imaginary
        if(not isinstance(X0, complex)): X0 *= 1j
        if(not isinstance(X1, complex)): X1 *= 1j
        if(not isinstance(X2, complex)): X2 *= 1j
        # Calculate Fault Currents
        If0 = 0
        If1 = Vsrc / (X1 + X2 + Rf)
        If2 = -If1
        faults = np.array([If0, If1, If2])
    # Return Currents
    if not sequence:
        faults = A012.dot(faults.T)
    return(faults)

# Define Three-Phase Fault Current Calculator
def phs3(Vsrc,Xseq,Rf=0,load=None,sequence=True):
    """
    PHS3 Function
    
    Purpose:
    --------
    This function will evaluate the 0, Positive, and Negative
    sequence currents for a three-phase fault with the
    option of calculating with or without a load.
    
    Required Arguments:
    -------------------
    Vsrc:       The Source Voltage
    Xseq:       Tupple of sequence reactances as (X0, X1, X2)
    
    Optional Arguments:
    -------------------
    Rf:         The fault resistance, default=0
    load:       The load conditions, default=None
    sequence:   Control argument to force return into ABC-Domain Currents
    
    Returns:
    --------
    Ifault: The Fault Current, equal for 0, pos., and neg. seq.
    """
    # Decompose Reactance Tuple
    X0, X1, X2 = Xseq
    # Load Condition
    if(load!=None):
        print("nothing yet")
    else:
        # Ensure that X-components are imaginary
        if(not isinstance(X1, complex)): X1 *= 1j
        # Calculate Fault Currents
        Ifault = Vsrc/(X1 + Rf)
        Ifault = np.array([ 0, Ifault, 0 ])
    # Prepare to Return Value
    if not sequence:
        Ifault = A012.dot( Ifault ) # Convert to ABC-Domain
    return(Ifault)


# Define Faulted Bus Voltage Calculator
def busvolt(k,n,Vpf,Z0,Z1,Z2,If,sequence=True):
    """
    BUSVOLT Function
    
    Purpose:
    --------
    This function is designed to calculate the bus voltage(s)
    given a specific set of fault characteristics.
    
    Required Arguments:
    -------------------
    k:          Bus at which to calculate faulted voltage
    n:          Bus at which fault occurred
    Vpf:        Voltage Pre-Fault, Singular Number
    Z0:         Zero-Sequence Impedance Matrix
    Z1:         Positive-Sequence Impedance Matrix
    Z2:         Negative-Sequence Impedance Matrix
    If:         Sequence Fault Current Evaluated at Bus *n*
    
    Optional Arguments:
    -------------------
    sequence:   Control argument to force return into ABC-Domain Currents
    
    Returns:
    --------
    Vf:         The Fault Voltage, set of sequence or phase voltages as
                specified by *sequence*
    """
    # Condition Inputs
    k = k-1
    n = n-1
    Z0 = np.asarray(Z0)
    Z1 = np.asarray(Z1)
    Z2 = np.asarray(Z2)
    If = np.asarray(If)
    # Generate Arrays For Calculation
    Vfmat = np.array([0, Vpf, 0]).T
    Zmat = np.array([[Z0[k,n], 0, 0],
                     [0, Z1[k,n], 0],
                     [0, 0, Z2[k,n]]])
    # Perform Calculation
    Vf = Vfmat - Zmat.dot(If)
    if not sequence:
        Vf = A012.dot( Vf ) # Convert to ABC-Domain
    return(Vf)


# Define CT Saturation Function
def ct_saturation(XR,Imag,Vrated,Irated,CTR,Rb,Xb,remnance=0,freq=60,ALF=20):
    """
    ct_saturation Function
    
    A function to determine the saturation value and a boolean indicator
    showing whether or not CT is -in fact- saturated.
    
    Parameters
    ----------
    XR:         float
                The X-over-R ratio of the system.
    Imag:       float
                The (maximum) current magnitude to use for calculation,
                typically the fault current.
    Vrated:     float
                The rated voltage (accompanying the C-Class value) of
                the CT.
    Irated:     float
                The rated secondary current for the CT.
    CTR:        float
                The CT Ratio (primary/secondary, N) to be used.
    Rb:         float
                The total burden resistance in ohms.
    Xb:         float
                The total burden reactance in ohms.
    remnance:   float, optional
                The system flux remnance, default=0.
    freq:       float, optional
                The system frequency in Hz, default=60.
    ALF:        float, optional
                The Saturation Constant which must be satisfied,
                default=20.
    
    Returns
    -------
    result:     float
                The calculated Saturation value.
    saturation: bool
                Boolean indicator to mark presence of saturation.
    """
    # Define omega
    w = 2*np.pi*freq
    # Find Lb
    Lb = Xb/w
    # Re-evaluate Vrated
    Vrated = Vrated*(1-remnance)
    # Calculate each "term" (multiple)
    t1 = (1+XR)
    t2 = (Imag/(Irated*CTR))
    t3 = abs(Rb+1j*w*Lb)*100/Vrated
    # Evaluate
    result = t1*t2*t3
    # Test for saturation
    saturation = result >= ALF
    # Return Results
    return(result,saturation)


# Define C-Class Calculator
def ct_cclass(XR,Imag,Irated,CTR,Rb,Xb,remnance=0,freq=60,ALF=20):
    """
    ct_cclass Function
    
    A function to determine the C-Class rated voltage for a CT.
    
    Parameters
    ----------
    XR:         float
                The X-over-R ratio of the system.
    Imag:       float
                The (maximum) current magnitude to use for calculation,
                typically the fault current.
    Irated:     float
                The rated secondary current for the CT.
    CTR:        float
                The CT Ratio (primary/secondary, N) to be used.
    Rb:         float
                The total burden resistance in ohms.
    Xb:         float
                The total burden reactance in ohms.
    remnance:   float, optional
                The system flux remnance, default=0.
    freq:       float, optional
                The system frequency in Hz, default=60.
    ALF:        float, optional
                The Saturation Constant which must be satisfied,
                default=20.
    
    Returns
    -------
    c_class:    float
                The calculated C-Class rated voltage.
    """
    # Define omega
    w = 2*np.pi*freq
    # Find Lb
    Lb = Xb/w
    # Calculate each "term" (multiple)
    t1 = (1+XR)
    t2 = (Imag/(Irated*CTR))
    t3 = abs(Rb+1j*w*Lb)*100/ALF
    # Evaluate
    Vr_w_rem = t1*t2*t3
    c_class = Vr_w_rem/(1-remnance)
    # Return Result
    return(c_class)


# Define Saturation Voltage at Rated Burden
def ct_satratburden(Inom,VArat=None,ANSIv=None,ALF=20,):
    """
    ct_satratburden Function
    
    A function to determine the Saturation at rated burden.
    
    Parameters
    ----------
    Inom:       float
                Nominal Current
    VArat:      float, optional, exclusive
                The apparent power (VA) rating of the CT.
    ANSIv:      float, optional, exclusive
                The ANSI voltage requirement to meet.
    ALF:        float, optional
                Accuracy Limit Factor, default=20.
    
    Returns
    -------
    Vsat:       float
                The saturated voltage.
    """
    # Validate Inputs
    if VArat == None and ANSIv == None:
        raise ValueError("VArat or ANSIv must be specified.")
    elif VArat==None:
        # Calculate VArat from ANSIv
        Zrat = ANSIv/(20*Inom)
        VArat = Inom**2 * Zrat
    # Determine Vsaturation
    Vsat = ALF * VArat/Inom
    return(Vsat)


# Define CT Vpeak Formula
def ct_vpeak(Zb,Ip,N):
    """
    ct_vpeak Function
    
    Simple formula to calculate the Peak Voltage of a CT.
    
    Parameters
    ----------
    Zb:         float
                The burden impedance magnitude (in ohms).
    Ip:         float
                The peak current for the CT.
    N:          float
                The CTR turns ratio of the CT.
    
    Returns
    -------
    Vpeak:      float
                The peak voltage.
    """
    return(np.sqrt(3.5*Zb*Ip*N))


# Define Saturation Time Calculator
def ct_timetosat(Vknee,XR,Rb,CTR,Imax,ts=None,npts=100,freq=60,plot=False):
    """
    ct_timetosat Function
    
    Function to determine the "time to saturate" for an underrated C-Class
    CT using three standard curves described by Juergen Holbach.
    
    Parameters
    ----------
    Vknee:      float
                The knee-voltage for the CT.
    XR:         float
                The X-over-R ratio of the system.
    Rb:         float
                The total burden resistance in ohms.
    CTR:        float
                The CT Ratio (primary/secondary, N) to be used.
    Imax:       float
                The (maximum) current magnitude to use for calculation,
                typically the fault current.
    ts:         numpy.ndarray or float, optional
                The time-array or particular (floatint point) time at which
                to calculate the values. default=np.linspace(0,0.1,freq*npts)
    npts:       float, optional
                The number of points (per cycle) to calculate if ts is not
                specified, default=100.
    freq:       float, optional
                The system frequency in Hz, default=60.
    plot:       bool, optional
                Control argument to enable plotting of calculated curves,
                default=False.
    """
    # Calculate omega
    w = 2*np.pi*freq
    # Calculate Tp
    Tp = XR/w
    # If ts isn't specified, generate it
    if ts==None:
        ts = np.linspace(0,0.1,freq*npts)
    # Calculate inner term
    term = -XR*(np.exp(-ts/Tp)-1)
    # Calculate Vsaturation terms
    Vsat1 = Imax*Rb*(term+1)
    Vsat2 = Imax*Rb*(term-np.sin(w*ts))
    Vsat3 = Imax*Rb*(1-np.cos(w*ts))
    # If plotting requested
    if plot and isinstance(ts,np.ndarray):
        plt.plot(ts,Vsat1,label="Vsat1")
        plt.plot(ts,Vsat2,label="Vsat2")
        plt.plot(ts,Vsat3,label="Vsat3")
        plt.axhline(Vknee,label="V-knee",linestyle='--')
        plt.title("Saturation Curves")
        plt.xlabel("Time (ts)")
        plt.legend()
        plt.show()
    elif plot:
        print("Unable to plot a single point, *ts* must be a numpy-array.")
    # Determine the crossover points for each saturation curve
    Vsat1c = Vsat2c = Vsat3c = 0
    if isinstance(ts,np.ndarray):
        for i in range(len(ts)):
            if Vsat1[i]>Vknee and Vsat1c==0:
                Vsat1c = ts[i-1]
            if Vsat2[i]>Vknee and Vsat2c==0:
                Vsat2c = ts[i-1]
            if Vsat3[i]>Vknee and Vsat3c==0:
                Vsat3c = ts[i-1]
        results = (Vsat1c,Vsat2c,Vsat3c)
    else:
        results = (Vsat1,Vsat2,Vsat3)
    return(results)

# Define Function to Calculate TRV
def pktransrecvolt(C,L,R=0,VLL=None,VLN=None,freq=60):
    """
    pktransrecvolt Function
    
    Peak Transient Recovery Voltage calculation
    function, evaluates the peak transient
    recovery voltage (restriking voltage) and
    the Rate-of-Rise-Recovery Voltage.
    
    Parameters
    ----------
    C:          float
                Capacitance Value in Farads.
    L:          float
                Inductance in Henries.
    R:          float, optional
                The resistance of the system used for
                calculation, default=0.
    VLL:        float, exclusive
                Line-to-Line voltage, exclusive
                optional argument.
    VLN:        float, exclusive
                Line-to-Neutral voltage, exclusive
                optional argument.
    freq:       float, optional
                System frequency in Hz.
    
    Returns
    -------
    Vcpk:       float
                Peak Transient Recovery Voltage in volts.
    RRRV:       float
                The RRRV (Rate-of-Rise-Recovery Voltage)
                calculated given the parameters in volts
                per second.
    """
    # Evaluate alpha, omega-n, and fn
    alpha = R/(2*L)
    wn = 1/np.sqrt(L*C) - alpha
    fn = wn/(2*np.pi)
    # Evaluate Vm
    if VLL!=None:
        Vm = np.sqrt(2/3)*VLL
    elif VLN!=None:
        Vm = np.sqrt(2)*VLN
    else:
        raise ValueError("One voltage must be specified.")
    # Evaluate Vcpk (worst case)
    Vcpk = wn**2/(wn**2-2*np.pi*freq)*Vm*2
    # Evaluate RRRV
    RRRV = 2*Vm*fn/0.5
    return(Vcpk,RRRV)

# Define TRV Reduction Resistor Function
def trvresistor(C,L,reduction,Rd0=500,wd0=260e3,tpk0=10e-6):
    """
    trvresistor Function
    
    Function to find the resistor value that
    will reduce the TRV by a specified
    percentage.
    
    Parameters
    ----------
    C:          float
                Capacitance Value in Farads.
    L:          float
                Inductance in Henries.
    reduction:  float
                The percentage that the TRV
                should be reduced by.
    Rd0:        float, optional
                Damping Resistor Evaluation Starting Point
                default=500
    wd0:        float, optional
                Omega-d evaluation starting point, default=260*k
    tpk0:       float, optional
                Time of peak voltage evaluation starting point,
                default=10*u
    
    Returns
    -------
    Rd:         float
                Damping resistor value, in ohms.
    wd:         float
                Omega-d
    tpk:        float
                Time of peak voltage.
    """
    # Evaluate omega-n
    wn = 1/np.sqrt(L*C)
    # Generate Constant Factor
    fctr = (1-reduction)*2 - 1
    # Define Function Set
    def equations(data):
        Rd, wd, tpk = data
        X = np.sqrt(wn**2-(1/(2*Rd*C))**2) - wd
        Y = np.exp(-tpk/(2*Rd*C))-fctr
        Z = wd*tpk - np.pi
        return(X,Y,Z)
    Rd, wd, tpk = fsolve(equations, (Rd0,wd0,tpk0))
    return(Rd, wd, tpk)

# Define Natural Frequency/Resonant Frequency Calculator
def natfreq(C,L,Hz=True):
    """
    natfreq Function
    
    Evaluates the natural frequency (resonant frequency)
    of a circuit given the circuit's C and L values. Defaults
    to returning values in Hz, but may also return in rad/sec.
    
    Parameters
    ----------
    C:          float
                Capacitance Value in Farads.
    L:          float
                Inductance in Henries.
    Hz:         bool, optional
                Control argument to set return value in either
                Hz or rad/sec; default=True.
    
    Returns
    -------
    freq:       float
                Natural (Resonant) frequency, will be in Hz if
                argument *Hz* is set True (default), or rad/sec
                if argument is set False.
    """
    # Evaluate Natural Frequency in rad/sec
    freq = 1/np.sqrt(L*C)
    # Convert to Hz as requested
    if Hz:
        freq = freq / (2*np.pi)
    return(freq)

# Define Time-Overcurrent Trip Time Function
def toctriptime(I,Ipickup,TD,curve="U1"):
    """
    toctriptime Function
    
    Time-OverCurrent Trip Time Calculator, evaluates the time
    to trip for a specific TOC (51) element given the curve
    type, current characteristics and time-dial setting.
    
    Parameters
    ----------
    I:          float
                Measured Current in Amps
    Ipickup:    float
                Fault Current Pickup Setting (in Amps)
    TD:         float
                Time Dial Setting
    curve:      string, optional
                Name of specified TOC curve, may be entry from set:
                {U1,U2,U3,U4,U5,C1,C2,C3,C4,C5}, default=U1
    
    Returns
    -------
    tt:         float
                Time-to-Trip for characterized element.
    """
    # Condition Inputs
    curve = curve.upper()
    # Define Dictionary of Constants
    const = {   "U1" : {"A": 0.0104, "B": 0.2256, "P": 0.02},
                "U2" : {"A": 5.95, "B": 0.180, "P": 2.00},
                "U3" : {"A": 3.88, "B": 0.0963, "P": 2.00},
                "U4" : {"A": 5.67, "B": 0.352, "P": 2.00},
                "U5" : {"A": 0.00342, "B": 0.00262, "P": 0.02},
                "C1" : {"A": 0.14, "B":0, "P": 0.02},
                "C2" : {"A": 13.5, "B":0, "P": 2.00},
                "C3" : {"A": 80.0, "B":0, "P": 2.00},
                "C4" : {"A": 120.0, "B":0, "P": 2.00},
                "C5" : {"A": 0.05, "B":0, "P": 0.04}}
    # Load Constants
    A = const[curve]["A"]
    B = const[curve]["B"]
    P = const[curve]["P"]
    # Evaluate M
    M = I / Ipickup
    # Evaluate Trip Time
    tt = TD * (A/(M**P-1)+B)
    return(tt)

# Define Time Overcurrent Reset Time Function
def tocreset(I,Ipickup,TD,curve="U1"):
    """
    tocreset Function
    
    Function to calculate the time to reset for a TOC
    (Time-OverCurrent, 51) element.
    
    Parameters
    ----------
    I:          float
                Measured Current in Amps
    Ipickup:    float
                Fault Current Pickup Setting (in Amps)
    TD:         float
                Time Dial Setting
    curve:      string, optional
                Name of specified TOC curve, may be entry from set:
                {U1,U2,U3,U4,U5,C1,C2,C3,C4,C5}, default=U1
    
    Returns
    -------
    tr:         float
                Time-to-Reset for characterized element.
    """
    # Condition Inputs
    curve = curve.upper()
    # Define Dictionary of Constants
    C = {   "U1" : 1.08,"U2" : 5.95,"U3" : 3.88,
            "U4" : 5.67,"U5" : 0.323,"C1" : 13.5,
            "C2" : 47.3,"C3" : 80.0,"C4" : 120.0,
            "C5" : 4.85}
    # Evaluate M
    M = I / Ipickup
    # Evaluate Reset Time
    tr = TD * (C[curve]/(1-M**2))
    return(tr)

# Define Pickup Current Calculation
def pickup(Iloadmax,Ifaultmin,scale=0,printout=False,units="A"):
    """
    pickup Function
    
    Used to assist in evaluating an optimal phase-over-current pickup
    setting. Uses maximum load and minimum fault current to provide
    user assistance.
    
    Parameters
    ----------
    Iloadmax:   float
                The maximum load current in amps.
    Ifaultmin:  float
                The minimum fault current in amps.
    scale:      int, optional
                Control scaling to set number of significant figures.
                default=0
    printout:   boolean, optional
                Control argument to enable printing of intermediate
                stages, default=False.
    units:      string, optional
                String to be appended to any printed output denoting
                the units of which are being printed, default="A"
    
    Returns
    -------
    setpoint:   float
                The evaluated setpoint at which the function suggests
                the phase-over-current pickup setting be placed.
    """
    IL2 = 2*Iloadmax
    IF2 = Ifaultmin/2
    exponent = len(str(IL2).split('.')[0])
    setpoint = np.ceil(IL2*10**(-exponent+1+scale))*10**(exponent-1-scale)
    if printout:
        print("Range Min:",IL2,units,"\t\tRange Max:",IF2,units)
    if IF2 < setpoint:
        setpoint = IL2
        if IL2 > IF2:
            raise ValueError("Invalid Parameters.")
    if printout:
        print("Current Pickup:",setpoint,units)
    return(setpoint)

# Define Time-Dial Coordination Function
def tdcoordradial(I,CTI,Ipu_up,Ipu_dn,TDdn,curve="U1",scale=1,freq=60):
    """
    tdcoordradial Function
    
    Function to evaluate the Time-Dial (TD) setting in radial schemes
    where the Coordinating Time Interval (CTI) and the up/downstream
    pickup settings are known along with the TD setting for the
    downstream protection.
    
    Parameters
    ----------
    I:          float
                Measured fault current in Amps, typically set using the
                maximum fault current available.
    CTI:        float
                Coordinating Time Interval in cycles.
    Ipu_up:     float
                Pickup setting for upstream protection,
                specified in amps
    Ipu_dn:     float
                Pickup setting for downstream protection,
                specified in amps
    TDdn:       float
                Time-Dial setting for downstream protection,
                specified in seconds
    curve:      string, optional
                Name of specified TOC curve, may be entry from set:
                {U1,U2,U3,U4,U5,C1,C2,C3,C4,C5}, default=U1
    scale:      int, optional
                Scaling value used to evaluate a practical TD
                setting, default=1
    freq:       float, optional
                System operating frequency, default=60
    
    Returns
    -------
    TD:         float
                Calculated Time-Dial setting according to radial
                scheme logical analysis.
    """
    # Condition Inputs
    curve = curve.upper()
    CTI = CTI/freq # Evaluate in seconds from cycles
    # Define Dictionary of Constants
    const = {   "U1" : {"A": 0.0104, "B": 0.2256, "P": 0.02},
                "U2" : {"A": 5.95, "B": 0.180, "P": 2.00},
                "U3" : {"A": 3.88, "B": 0.0963, "P": 2.00},
                "U4" : {"A": 5.67, "B": 0.352, "P": 2.00},
                "U5" : {"A": 0.00342, "B": 0.00262, "P": 0.02},
                "C1" : {"A": 0.14, "B":0, "P": 0.02},
                "C2" : {"A": 13.5, "B":0, "P": 2.00},
                "C3" : {"A": 80.0, "B":0, "P": 2.00},
                "C4" : {"A": 120.0, "B":0, "P": 2.00},
                "C5" : {"A": 0.05, "B":0, "P": 0.04}}
    # Load Constants
    A = const[curve]["A"]
    B = const[curve]["B"]
    P = const[curve]["P"]
    # Evaluate M
    M = I / Ipu_dn
    # Evaluate Trip Time
    tpu_desired = TDdn * (A/(M**P-1)+B) + CTI
    # Re-Evaluate M
    M = I / Ipu_up
    # Calculate TD setting
    TD = tpu_desired / (A/(M**2-1)+B)
    # Scale and Round
    TD = np.ceil(TD*10**scale)/10**scale
    return(TD)


# END OF FILE