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
def phs2g(Vsrc,Xseq,Rf=0,load=None):
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
    Ifault: The Tupled Fault Currents as (If0, If1, If2)
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
        If1 = Vsrc / (X1 + (1/X2 + 1/(X0 + 3*Rf))**(-1))
        If2 = -If1 * (X0 + 3*Rf)/(X2 + X0 + 3*Rf)
        If0 = -If1 * X2/(X2 + X0 + 3*Rf)
    # Return Currents
    return(If0, If1, If2)

# Define Phase-to-Phase Fault Current Calculator
def phs2(Vsrc,Xseq,Rf=0,load=None):
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
    Ifault: The Tupled Fault Currents as (If0, If1, If2)
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
    return(If0, If1, If2)

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