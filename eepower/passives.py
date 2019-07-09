###################################################################
#   PASSIVES.PY
#
#   A library of functions, constants and more that are related to
#   Inductors and Capacitors in Electrical Engineering.
#
#   Written by Joe Stanley
#
#   Special Thanks To:
#   Paul Ortmann - Idaho Power
#   Brian Johnson - University of Idaho
#
#   Included Functions
#   - Inductor Charge:              inductorcharge
#   - Inductor Discharge:           inductordischarge
#   - Inductor Stored Energy:       inductorenergy      
#   - Back-to-Back Cap. Surge:      capbacktoback  
#   - Capacitor Stored Energy:      energy
#   - Cap. Voltage after Time:      VafterT
#   - Cap. Voltage Discharge:       vcapdischarge
#   - Cap. Voltage Charge:          vcapcharge
#   - Rectifier Cap. Calculation:   rectifier
#   - Cap. VAR to FARAD Conversion: farads
#   DOCUMENTATION NOT UP TO DATE!!!!!
####################################################################

# Import libraries as needed:
import numpy as np

# Define Capacitor Voltage Discharge Function
def vcapdischarge(t,Vs,R,C):
    """
    vcapdischarge Function
    
    Function to calculate the voltage of a
    capacitor that is discharging given the time.
    
    Parameters
    ----------
    t:          float
                The time at which to calculate the voltage.
    Vs:         float
                The starting voltage for the capacitor.
    R:          float
                The ohmic value of the resistor being used
                to discharge.
    C:          float
                Capacitive value (in Farads).
    
    Returns
    -------
    Vc:         float
                The calculated voltage of the capacitor.
    """
    Vc = Vs*(np.exp(-t/(R*C)))
    return(Vc)

# Define Capacitor Voltage Charge Function
def vcapcharge(t,Vs,R,C):
    """
    vcapcharge Function
    
    Function to calculate the voltage of a
    capacitor that is charging given the time.
    
    Parameters
    ----------
    t:          float
                The time at which to calculate the voltage.
    Vs:         float
                The charging voltage for the capacitor.
    R:          float
                The ohmic value of the resistor being used
                to discharge.
    C:          float
                Capacitive value (in Farads).
    
    Returns
    -------
    Vc:         float
                The calculated voltage of the capacitor.
    """
    Vc = Vs*(1-np.exp(-t/(R*C)))
    return(Vc)
    
# Define Capacitive Energy Transfer Function
def captransfer(t,Vs,R,Cs,Cd):
    """
    captransfer Function
    
    Calculate the voltage across a joining
    resistor (R) that connects Cs and Cd, the
    energy-source and -destination capacitors,
    respectively. Calculate the final voltage
    across both capacitors.
    
    Parameters
    ----------
    t:          float
                Time at which to calculate resistor voltage.
    Vs:         float
                Initial voltage across source-capacitor (Cs).
    R:          float
                Value of resistor that connects capacitors.
    Cs:         float
                Source capacitance value in Farads.
    Cd:         float
                Destination capacitance value in Farads.
    
    Returns
    -------
    rvolt:      float
                Voltage across the resistor at time t.
    vfinal:     float
                Final voltage that both capacitors settle to.
    """
    tau = (R*Cs*Cd) / (Cs+Cd)
    rvolt = Vs*np.exp(-t/tau)
    vfinal = Vs*Cs/(Cs+Cd)
    return(rvolt,vfinal)
    
# Define Inductor Energy Formula
def inductorenergy(L,I):
    """
    inductorenergy Function
    
    Function to calculate the energy stored in an inductor
    given the inductance (in Henries) and the current.
    
    Parameters
    ----------
    L:          float
                Inductance Value (in Henries)
    I:          float
                Current traveling through inductor.
    
    Returns
    -------
    E:          float
                The energy stored in the inductor (in Joules).
    """
    return(1/2 * L * I**2)

# Define Inductor Charge Function
def inductorcharge(t,Vs,R,L):
    """
    inductorcharge Function
    
    Calculates the Voltage and Current of an inductor
    that is charging/storing energy.
    
    Parameters
    ----------
    t:          float
                Time at which to calculate voltage and current.
    Vs:         float
                Charging voltage across inductor and resistor.
    R:          float
                Resistance related to inductor.
    L:          float
                Inductance value in Henries.
    
    Returns
    -------
    V1:         float
                Voltage across inductor at time t.
    I1:         float
                Current through inductor at time t.
    """
    Vl = Vs*np.exp(-R*t/L)
    Il = Vs/R*(1-np.exp(-R*t/L))
    return(Vl,Il)

# Define Capacitive Back-to-Back Switching Formula
def capbacktoback(C1,C2,Lm,VLN=None,VLL=None):
    """
    capbacktoback Function
    
    Function to calculate the maximum current and the 
    frequency of the inrush current of two capacitors
    connected in parallel when one (energized) capacitor
    is switched into another (non-engergized) capacitor.
    
    Note: This formula is only valid for three-phase systems.
    
    Parameters
    ----------
    C1:         float
                The capacitance of the
    VLN:        float, exclusive
                The line-to-neutral voltage experienced by
                any one of the (three) capacitors in the
                three-phase capacitor bank.
    VLL:        float, exclusive
                The line-to-line voltage experienced by the
                three-phase capacitor bank.
    """
    # Evaluate Max Current
    imax = np.sqrt(2/3)*VLL*np.sqrt((C1*C2)/((C1+C2)*Lm))
    # Evaluate Inrush Current Frequency
    ifreq = 1/(2*np.pi*np.sqrt(Lm*(C1*C2)/(C1+C2)))
    return(imax,ifreq)

# Define Inductor Discharge Function
def inductordischarge(t,Io,R,L):
    """
    inductordischarge Function
    
    Calculates the Voltage and Current of an inductor
    that is discharging its stored energy.
    
    Parameters
    ----------
    t:          float
                Time at which to calculate voltage and current.
    Io:         float
                Initial current traveling through inductor.
    R:          float
                Resistance being discharged to.
    L:          float
                Inductance value in Henries.
    
    Returns
    -------
    V1:         float
                Voltage across inductor at time t.
    I1:         float
                Current through inductor at time t.
    """
    Il = Io*np.exp(-R*t/L)
    Vl = Io*R*(1-np.exp(-R*t/L))
    return(Vl,Il)
    
# Define Apparent Power to Farad Conversion
def farads(VAR,V,freq=60):
    """
    farads Formula
    
    Function to calculate the required capacitance
    in Farads to provide the desired power rating
    (VARs).
    
    Parameters
    ----------
    VAR:        float
                The rated power to meet.
    V:          float
                The voltage across the capacitor;
                not described as VLL or VLN, merely
                the capacitor voltage.
    freq:       float, optional
                The System frequency
    
    Returns
    -------
    C:          float
                The evaluated capacitance (in Farads).
    """
    return(VAR / (2*np.pi*freq*V**2))

# Define Capacitor Energy Calculation
def capenergy(C,v):
    """
    capenergy Function
    
    A simple function to calculate the stored voltage (in Joules)
    in a capacitor with a charged voltage.
    
    Parameters
    ----------
    C:          float
                Capacitance in Farads.
    v:          float
                Voltage across capacitor.
    
    Returns
    -------
    energy:     float
                Energy stored in capacitor (Joules).
    """
	energy = 1/2 * C * v**2
	return(energy)

# Define Capacitor Voltage Discharge Function
def loadedvcapdischarge(t,vo,C,P):
    """
    loadedvcapdischarge Function
    
    Returns the voltage of a discharging capacitor after time (t - 
    seconds) given initial voltage (vo - volts), capacitor size
    (cap - Farads), and load (P - Watts).
    
    Parameters
    ----------
    t:          float
                Time at which to calculate voltage.
    vo:         float
                Initial capacitor voltage.
    C:          float
                Capacitance (in Farads)
    P:          float
                Load power consumption (in Watts).
    
    Returns
    -------
    Vt:         float
                Voltage of capacitor at time t.
    """
	Vt = np.sqrt(vo**2 - 2*P*t/C)
	return(Vt)
	
# Define Capacitor Discharge Function
def timedischarge(Vinit,Vmin,C,P,dt=1e-3,RMS=True,Eremain=False):
    """
    timedischarge Function
    
    Returns the time to discharge a capacitor to a specified
    voltage given set of inputs.
    
    Parameters
    ----------
    Vinit:      float
                Initial Voltage (in volts)
    Vmin:       float
                Final Voltage (the minimum allowable voltage) (in volts)
    C:          float
                Capacitance (in Farads)
    P:          float
                Load Power being consumed (in Watts)
    dt:         float, optional
                Time step-size (in seconds) (defaults to 1e-3 | 1ms)
    RMS:        bool, optional
                if true converts RMS Vin to peak
    Eremain:    bool, optional
                if true: also returns the energy remaining in cap
    
    Returns
    -------
    Returns time to discharge from Vinit to Vmin in seconds.
    May also return remaining energy in capacitor if Eremain=True
    """
    t = 0 # start at time t=0
    if RMS:
        vo = Vinit*np.sqrt(2) # convert RMS to peak
    else:
        vo = Vinit
    vc = loadedvcapdischarge(t,vo,C,P) # set initial cap voltage
    while(vc >= Vmin):
        t = t+dt # increment the time
        vcp = vc # save previous voltage
        vc = loadedvcapdischarge(t,vo,C,P) # calc. new voltage
    if(Eremain):
        E = energy(C,vcp) # calc. energy
        return(t-dt,E)
    else:
        return(t-dt)


# Define Rectifier Capacitor Calculator
def rectifier(Iload, fswitch, dVout):
    """
    rectifier Function
    
    Returns the capacitance (in Farads) for a needed capacitor in
    a rectifier configuration given the system frequency (in Hz),
    the load (in amps) and the desired voltage ripple.
    
    Parameters
    ----------
    Iload:      float
                The load current that must be met.
    fswitch:    float
                The switching frequency of the system.
    dVout:      float
                Desired delta-V on the output.
    
    Returns
    -------
    C:          float
                Required capacitance (in Farads) to meet arguments.
    """
	C = Iload / (fswitch * dVout)
	return(C)

    
# END OF FILE