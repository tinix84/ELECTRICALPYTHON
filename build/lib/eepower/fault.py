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
#   - 
####################################################################

# Import Necessary Libraries
import numpy as np

# Import Local Dependencies
from .__init__ import Aabc, A012

# Define Single Line to Ground Fault Function
def phs1g(Vsrc,Xseq,Rf=0,load=None):
    """
    PHS1G Function
    
    Purpose:
    --------
    This function will evaluate the 0, Positive, and Negative
    sequence currents for a single-line-to-ground fault with the
    option of calculating with or without a load.
    
    Required Arguments:
    -------------------
    Vsrc:   The Source Voltage
    Xseq:   Tupple of sequence reactances as (X0, X1, X2)
    
    Optional Arguments:
    -------------------
    Rf:     The fault resistance, default=0
    load:   The load conditions, default=None
    
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
    Vsrc:   The Source Voltage
    Xseq:   Tupple of sequence reactances as (X0, X1, X2)
    
    Optional Arguments:
    -------------------
    Rf:     The fault resistance, default=0
    load:   The load conditions, default=None
    
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
    Vsrc:   The Source Voltage
    Xseq:   Tupple of sequence reactances as (X0, X1, X2)
    
    Optional Arguments:
    -------------------
    Rf:     The fault resistance, default=0
    load:   The load conditions, default=None
    
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
def phs3(Vsrc,Xseq,Rf=0,load=None):
    """
    PHS3 Function
    
    Purpose:
    --------
    This function will evaluate the 0, Positive, and Negative
    sequence currents for a three-phase fault with the
    option of calculating with or without a load.
    
    Required Arguments:
    -------------------
    Vsrc:   The Source Voltage
    Xseq:   Tupple of sequence reactances as (X0, X1, X2)
    
    Optional Arguments:
    -------------------
    Rf:     The fault resistance, default=0
    load:   The load conditions, default=None
    
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
        # Calculate Fault Currents
        Ifault = Vsrc/X1
    return(Ifault)


# END OF FILE