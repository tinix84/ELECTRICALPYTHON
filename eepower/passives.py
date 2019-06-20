###################################################################
#   PASSIVES.PY
#
#   A library of functions, constants and more
#   that are related to Capacitors in Electrical Engineering.
#
#   February 13, 2019
#
#   Written by Joe Stanley
#
#   Special Thanks To:
#   Paul Ortmann - Idaho Power
#
#   Included Functions
#   - Capacitor Stored Energy:      energy
#   - Cap. Voltage after Time:      VafterT
#   - Cap. Voltage Discharge:       discharge
#   - Rectifier Cap. Calculation:   rectifier
#   DOCUMENTATION NOT UP TO DATE!!!!!
####################################################################

# Import libraries as needed:
import numpy as np

def vcapdischarge(t,Vs,R,C):
    Vc = Vs*(np.exp(-t/(R*C)))
    return(Vc)


def vcapcharge(t,Vs,R,C):
    Vc = Vs*(1-np.exp(-t/(R*C)))
    return(Vc)
    
def captransfer(t,Vs,R,Cs,Cd):
    tau = (R*Cs*Cd) / (Cs+Cd)
    rvolt = Vs*np.exp(-t/tau)
    vfinal = Vs*Cs/(Cs+Cd)
    return(rvolt,vfinal)
    
# Define Inductor Energy Formula
def inductorenergy(L,I):
    return(1/2 * L * I**2)

def inductorcharge(t,Vs,R,L):
    Vl = Vs*np.exp(-R*t/L)
    Il = Vs/R*(1-np.exp(-R*t/L))
    return(Vl,Il)

# Define Capacitive Back-to-Back Switching Formula
def capbacktoback(VLL,C1,C2,Lm):
    # Evaluate Max Current
    imax = np.sqrt(2/3)*VLL*np.sqrt((C1*C2)/((C1+C2)*Lm))
    # Evaluate Inrush Current Frequency
    ifreq = 1/(2*np.pi*np.sqrt(Lm*(C1*C2)/(C1+C2)))
    return(imax,ifreq)

def inductordischarge(t,Io,R,L):
    Il = Io*np.exp(-R*t/L)
    Vl = Io*R*(1-np.exp(-R*t/L))
    return(Vl,Il)

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
	
###################################################################
#   Define Capacitor Discharge Function
#
#   Returns the time to discharge a capacitor to a specified
#   voltage given set of inputs:
#
#   Vinit: Initial Voltage (in volts)
#   Vmin: Final Voltage (the minimum allowable voltage) (in volts)
#   cap: Capacitance (in Farads)
#   P: Load Power being consumed (in Watts)
#   dt: Time step-size (in seconds) (defaults to 1e-3 | 1ms)
#   RMS: if true converts RMS Vin to peak
#   Eremain: if true: also returns the energy remaining in cap
#
#   Returns time to discharge from Vinit to Vmin in seconds.
#   May also return remaining energy in capacitor if Eremain=True
###################################################################
def discharge(Vinit,Vmin,cap,P,dt=1e-3,RMS=True,Eremain=False):
    t = 0 # start at time t=0
    if RMS:
        vo = Vinit*np.sqrt(2) # convert RMS to peak
    else:
        vo = Vinit
    vc = VafterT(t,vo,cap,P) # set initial cap voltage
    while(vc >= Vmin):
        t = t+dt # increment the time
        vcp = vc # save previous voltage
        vc = VafterT(t,vo,cap,P) # calc. new voltage
    if(Eremain):
        E = energy(cap,vcp) # calc. energy
        return(t-dt,E)
    else:
        return(t-dt)

###################################################################
#   Define Rectifier Capacitor Calculator
#
#   Returns the capacitance (in Farads) for a needed capacitor in
#   a rectifier configuration given the system frequency (in Hz),
#   the load (in amps) and the desired voltage ripple.
###################################################################
def rectifier(Iload, fsys, dVout):
	C = Iload / (fsys * dVout)
	return(C)