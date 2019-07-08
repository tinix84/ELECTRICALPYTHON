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
    
def captransfer(t,Vs,R,Cs,Cd):
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
                The energy stored in the inductor.
    """
    return(1/2 * L * I**2)

def inductorcharge(t,Vs,R,L):
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

def inductordischarge(t,Io,R,L):
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

###################################################################
#   Define Capacitor Energy Calculation
#
#   Returns energy (in Joules) of a capacitor given capacitor size
#   (in Farads) and voltage (in Volts).
###################################################################
def capenergy(cap,v):
	energy = 1/2 * cap * v**2
	return(energy)

###################################################################
#   Define Capacitor Voltage Discharge Function
#
#   Returns the voltage of a discharging capacitor after time (t - 
#   seconds) given initial voltage (vo - volts), capacitor size
#   (cap - Farads), and load (P - Watts).
###################################################################
def VafterT(t,vo,cap,P):
	Vt = np.sqrt(vo**2 - 2*P*t/cap)
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
    vc = VafterT(t,vo,C,P) # set initial cap voltage
    while(vc >= Vmin):
        t = t+dt # increment the time
        vcp = vc # save previous voltage
        vc = VafterT(t,vo,C,P) # calc. new voltage
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