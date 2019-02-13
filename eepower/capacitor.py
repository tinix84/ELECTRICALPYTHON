###################################################################
#   CAPACITOR.PY
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
####################################################################

# Import libraries as needed:
import numpy as np

###################################################################
#   Define Capacitor Energy Calculation
#
#   Returns energy (in Joules) of a capacitor given capacitor size
#   (in Farads) and voltage (in Volts).
###################################################################
def energy(cap,v):
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
    vc = C_VafterT(t,vo,cap,P) # set initial cap voltage
    while(vc >= Vmin):
        t = t+dt # increment the time
        vcp = vc # save previous voltage
        vc = C_VafterT(t,vo,cap,P) # calc. new voltage
    if(Eremain):
        E = C_energy(cap,vcp) # calc. energy
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