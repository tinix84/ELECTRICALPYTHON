###################################################################
#   ELECTRONICS.PY
#
#   A library of functions, constants and more that are related
#   to Power-Electronics and converters in Electrical Engineering.
#
#   April 10th, 2019
#
#   Written by Joe Stanley
#
#   Special Thanks To:
#   Dr. Brian Johnson - Faculty University of Idaho
#
#   Included Functions:
#   - VSC DC Bus Voltage Calculator:        vscdcbus
#   - PLL-VSC Gains Calculator:             vscgains
####################################################################

# Import the Necessary Library
import numpy as np

# Define function to find VDC setpoint
def vscdcbus(VLL,Zs,P,Q=0,mmax=0.8,debug=False):
    """
    VSCDCBUS Function:
    
    Purpose:
    --------
    The purpose of this function is to calculate the
    required DC-bus voltage for a Voltage-Sourced-
    Converter (VSC) given the desired P/Q parameters
    and the known source impedance (Vs) of the VSC.
    
    Required Arguments:
    -------------------
    VLL:    Line-to-Line voltage on the line-side of
            the source impedance.
    Zs:     The source impedance of the VSC
    P:      The desired real-power output
    
    Optional Arguments:
    Q:      The desired reactive-power output, default=0
    mmax:   The maximum of the m value for the converter
            default=0.8
    debug:  Control value to enable printing stages of
            the calculation, default=False
            
    Return:
    -------
    VDC:    The DC bus voltage.
    """
    # Determine the Load Current
    Iload = np.conj((P+1j*Q) / (VLL*np.sqrt(3)))
    # Evaluate the Terminal Voltage
    Vtln = abs(VLL/np.sqrt(3) + Iload*Zs)
    # Find the Peak Terminal Voltage
    Vtpk = np.sqrt(2)*Vtln
    # Calculate the VDC value
    VDC = 2*Vtpk / mmax
    if debug:
        print("Iload", Iload)
        print("Vtln", Vtln)
        print("Vtpk", Vtpk)
        print("VDC", VDC)
    return(VDC)

# Define kp/ki/w0L calculating function
def vscgains(Rs,Ls,tau=0.005,f=60):
    """
    VSCGAINS Function:
    
    Purpose:
    --------
    This function is designed to calculate the kp, ki,
    and omega-not-L values for a Phase-Lock-Loop based VSC.
    
    Required Arguments:
    -------------------
    Rs:      The equiv-resistance (in ohms) of the VSC
    Ls:      The equiv-inductance (in Henrys) of the VSC
    
    Optional Arguments:
    -------------------
    tau:     The desired time-constant, default=0.005
    f:       The system frequency (in Hz), default=60
    
    Returns:
    --------
    kp:      The Kp-Gain Value
    ki:      The Ki-Gain Value
    w0L:     The omega-not-L gain value
    """
    # Calculate kp
    kp = Ls / tau
    # Calculate ki
    ki = kp*Rs/Ls
    # Calculate w0L
    w0L = 2*np.pi*f*Ls
    return(kp,ki,w0L)