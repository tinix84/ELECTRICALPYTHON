#################################################################################
#   BODE.PY
#
#   This file contains a set of bode plot functions related to signals.
#   These items will commonly be used in Electrical Engineering Applications.
#
#   February 21, 2019
#
#   Written by Joe Stanley
#   Special Thanks To and Code Support From:
#   Dr. Dennis Sullivan
#
#   Included Functions:
#   - Transfer Function Bode Plot Generator         bode
#   - S-Domain Bode Plot Generator                  sbode
#   - Z-Domain Bode Plot Generator                  zbode
#################################################################################

# Import External Dependencies
import matplotlib.pyplot as plt
import numpy as np
import scipy.signal as sig
from numpy import pi
from cmath import exp

# Import Local Dependencies
from .filter import sys_condition

# Define System Bode Plotting Function
def bode(system,mn=-2,mx=3,npts=100,gtitle="",xlim=False,ylim=False,sv=False,
         disp3db=False,lowcut=None,magnitude=True,angle=True):
    """ System Bode Plotting Function
    
    A simple function to generate the Bode Plot for magnitude
    and frequency given a transfer function system.
    
    Required Arguments
    ------------------
    system:        The Transfer Function; can be provided as the following:
                - 1 (instance of lti)
                - 2 (num, den)
                - 3 (zeros, poles, gain)
                - 4 (A, B, C, D)
                
    Optional Arguments
    ------------------
    mn:            The minimum frequency (as an exponent to 10, e.g. 10^mn)
                   to be calculated for. default=-2.
    mx:            The maximum frequency (as an exponent to 10, e.g. 10^mx)
                   to be calculated for. default=3.
    npts:          The number of points over which to calculate the system.
                   default=100.
    gtitle:        Additional string to be added to plot titles;
                   default="".
    xlim:          Limit in x-axis for graph plot. Accepts tuple of: (xmin, xmax).
                   Default is False.
    ylim:          Limit in y-axis for graph plot. Accepts tuple of: (ymin, ymax).
                   Default is False.
    sv:            Save the plots as PNG files. Default is False.
    disp3db:       Control argument to enable the display of the 3dB line,
                   default=False.
    lowcut:        An additional marking line that can be plotted, default=None
    magnitude:     Control argument to enable plotting of magnitude, default=True
    angle:         Control argument to enable plotting of angle, default=True
    
    Returns
    -------
    NONE:    Generates plot of magnitude and phase, does not return
            any numerical values.
    """
    # Condition system input to ensure proper execution
    system = sys_condition(system,False)
    
    # Generate the frequency range to calculate over
    wover = np.logspace(mn,mx,npts)
    
    # Calculate the bode system
    w, mag, ang = sig.bode(system, wover)
    
    # Plot Magnitude
    if(magnitude):
        magTitle = "Magnitude "+gtitle
        plt.title(magTitle)
        plt.plot(w, mag)
        plt.xscale("log")
        plt.grid(which="both")
        plt.ylabel("Magnitude (dB)")
        plt.xlabel("Frequency (rad/sec)")
        if disp3db:
            plt.axhline(-3)
        if lowcut!=None:
            plt.axhline(lowcut)
        if xlim!=False:
            plt.xlim(xlim)
        if ylim!=False:
            plt.ylim(ylim)
        if sv:
            plt.savefig(magTitle+".png")
        plt.show()

    # Plot Angle
    if(angle):
        angTitle = "Angle "+gtitle
        plt.title(angTitle)
        plt.plot(w, ang)
        plt.xscale("log")
        plt.grid(which="both")
        plt.ylabel("Angle (degrees)")
        plt.xlabel("Frequency (rad/sec)")
        if xlim!=False:
            plt.xlim(xlim)
        if ylim!=False:
            plt.ylim(ylim)
        if sv:
            plt.savefig(angTitle+".png")
        plt.show()

def sbode(f,NN=1000,gtitle="",xlim=False,ylim=False,mn=0,mx=1000,
          disp3db=False,lowcut=None,magnitude=True,angle=True):
    """
    SBODE Function
    
    Required Arguments:
    -------------------
    f:    The Input Function, must be callable
    
    Optional Arguments:
    -------------------
    NN:            The Interval over which to be generated, default=1000
    gtitle:        Additional string to be added to plot titles;
                   default="".
    xlim:          Limit in x-axis for graph plot. Accepts tuple of: (xmin, xmax).
                   Default is False.
    ylim:          Limit in y-axis for graph plot. Accepts tuple of: (ymin, ymax).
                   Default is False.
    mn:            The minimum W value to be generated, default=0
    mx:            The maximum W value to be generated, default=1000
    sv:            Save the plots as PNG files. Default is False.
    disp3db:       Control argument to enable the display of the 3dB line,
                   default=False.
    lowcut:        An additional marking line that can be plotted, default=None
    magnitude:     Control argument to enable plotting of magnitude, default=True
    angle:         Control argument to enable plotting of angle, default=True
    
    Returns:
    --------
    NONE - Plots Generated
    """
    W = np.linspace(mn,mx,NN)
    H = np.zeros(NN, dtype = np.complex)

    for n in range(0,NN):
        s = 1j*W[n]
        H[n] = f(s)
    if(magnitude):
        plt.semilogx(W,20*np.log10(abs(H)),'k')
        plt.ylabel('|H| dB')
        plt.xlabel('Frequency (rad/sec)')
        plt.title(gtitle+" Magnitude")
        plt.grid(which='both')
        if disp3db:
            plt.axhline(-3)
        if lowcut!=None:
            plt.axhline(lowcut)
        if xlim!=False:
            plt.xlim(xlim)
        if ylim!=False:
            plt.ylim(ylim)
        plt.show()

    aaa = np.angle(H)
    for n in range(NN):
        if aaa[n] > pi:
            aaa[n] = aaa[n] - 2*pi

    if(angle):
        plt.title(gtitle+" Phase")
        plt.semilogx(W,(180/pi)*aaa,'k')
        plt.ylabel('H phase (degrees)')
        plt.xlabel('Frequency (rad/sec)')
        plt.grid(which='both')
        if xlim!=False:
            plt.xlim(xlim)
        if ylim!=False:
            plt.ylim(ylim)
        plt.show()


def zbode(f,dt=0.01,NN=1000,gtitle="",mn=0,mx=2*pi,xlim=False,ylim=False,
          approx=False,disp3db=False,lowcut=None,magnitude=True,angle=True):
    """
    ZBODE Function
    
    Required Arguments:
    -------------------
    f:      The Input Function, must be callable.
            Must be specified as transfer function of type:
                - S-Domain (when approx=False, default)
                - Z-Domain (when approx=True)
    
    Optional Arguments:
    -------------------
    dt:            The time-step used, default=0.01
    NN:            The Interval over which to be generated, default=1000
    mn:            The minimum phi value to be generated, default=0
    mx:            The maximum phi value to be generated, default=2*pi
    approx:        Control argument to specify whether input funciton
                   should be treated as Z-Domain function or approximated
                   Z-Domain function. default=False
    gtitle:        Additional string to be added to plot titles;
                   default="".
    xlim:          Limit in x-axis for graph plot. Accepts tuple of: (xmin, xmax).
                   Default is False.
    ylim:          Limit in y-axis for graph plot. Accepts tuple of: (ymin, ymax).
                   Default is False.
    sv:            Save the plots as PNG files. Default is False.
    disp3db:       Control argument to enable the display of the 3dB line,
                   default=False.
    lowcut:        An additional marking line that can be plotted, default=None
    magnitude:     Control argument to enable plotting of magnitude, default=True
    angle:         Control argument to enable plotting of angle, default=True
    
    Returns:
    --------
    NONE - Plots Generated
    """
    phi = np.linspace(mn,mx,NN)

    H = np.zeros(NN, dtype = np.complex)
    for n in range(0,NN):
        z = exp(1j*phi[n])
        if(approx!=False): # Approximated Z-Domain
            s = approx(z,dt) # Pass current z-value and dt
            H[n] = f(s)
        else: # Z-Domain Transfer Function Provided
            H[n] = dt*f(z)
            
    if(magnitude):
        plt.semilogx((180/pi)*phi,20*np.log10(abs(H)),'k')
        plt.ylabel('|H| dB')
        plt.xlabel('Frequency (degrees)')
        plt.title(gtitle+" Magnitude")
        plt.grid(which='both')
        if disp3db:
            plt.axhline(-3)
        if lowcut!=None:
            plt.axhline(lowcut)
        if xlim!=False:
            plt.xlim(xlim)
        if ylim!=False:
            plt.ylim(ylim)
        plt.show()

    aaa = np.angle(H)
    for n in range(NN):
        if aaa[n] > pi:
            aaa[n] = aaa[n] - 2*pi

    if(angle):
        plt.semilogx((180/pi)*phi,(180/pi)*aaa,'k')
        plt.ylabel('H (degrees)')
        plt.grid(which='both')
        plt.xlabel('Frequency (degrees)')
        plt.title(gtitle+" Phase")
        if xlim!=False:
            plt.xlim(xlim)
        if ylim!=False:
            plt.ylim(ylim)
        plt.show()


# End of BODE.PY