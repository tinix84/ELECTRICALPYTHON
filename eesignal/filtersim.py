#################################################################################
#   FILTERSIM.PY
#
#   This file contains a variety of functions and constants related to signals.
#   These items will commonly  be used in Electrical Engineering Applications.
#
#   February 21, 2019
#
#   Written by Joe Stanley
#   Special Thanks To and Code Support From:
#   Dr. Dennis Sullivan
#
#   Included Functions:
#   - Z-Domain Filter Simulator       zfiltersim
#   - System Response Plotter:        sysresponse
#   - State Space Simulator:          statespace
#
#   Local Dependencies:
#   - Tupple to Matrix Converter:       tuple_to_matrix
#   - Numpy Array to Matrix Converter:  nparr_to_matrix
#   - Function Concatinator:            c_func_concat
#################################################################################

# Import Required External Libraries
import numpy as np
import matplotlib.pyplot as plt

# Import Local Dependencies
from .__init__ import tuple_to_matrix, nparr_to_matrix, c_func_concat

# Define constants
NAN = float('nan')
matrix = "<class 'numpy.matrixlib.defmatrix.matrix'>"
tuple = "<class 'tuple'>"
ndarr = "<class 'numpy.ndarray'>"
tint = "<class 'int'>"
tfloat = "<class 'float'>"
tfun = "<class 'function'>"
tnfloat = "<class 'numpy.float64'>"

def zfiltersim( fin, filter, NN=1000, title="",plotinput=True,legend=True):
    """
    ZFILTERSIM Function
    
    Purpose:
    --------
    Given an input function and filter parameters (specified in
    the z-domain) this function will plot the input function over
    NN time-steps of an unspecified size (the step-size must be
    specified outside the scope of this function). This function
    will also plot the resultant (output) function over the NN
    time-steps after having been filtered by that filter which
    is specified.
    
    Required Arguments:
    -------------------
    fin:    The input function, must be callable with specified
            step-size.
    filter: The filter parameter set as described at the end of
            this help/comment section.
    
    Optional Arguments:
    NN:         The number of time-steps to be plotted; default=1000
    title:      The title presented on each plot; default=""
    plotinput:  An argument to control whether the input is plotted
                separately, default=True.
    legend:     An argument to control whether the legend is shown,
                default=True.
    
    Returns:
    --------
    NONE, this function generates plots.
    
    # NOTICE: ------------------------------------
    # The *filter* argument should be provided as:
    # [[ a11, a12, b10, b11, b12],
    #  [ a21, a22, b20, b21, b22],
    #  [           ...          ],
    #  [ an1, an2, bn0, bn1, bn2]]
    # Where each row corresponds to a 1- or 2-pole
    # filter with the structure as follows:
    #          b0 + b1z^-1 + b2z^-2
    # H(z) = -----------------------------
    #           1 - a1z^-1 - a2z^-2
    # --------------------------------------------
    """
    # Start with arrays set to zero
    x = np.zeros(NN)
    y = np.zeros(NN)
    
    # ----- The input  -----
    for k in range(NN):
        x[k] = fin(k)
    if plotinput:
        plt.figure()
        plt.title(title)
        plt.plot(x)
        plt.show()
    
    # Identify how many rows were provided
    rows, cols = filter.shape
    # Operate with each individual filter set
    x_tmp = np.copy( x )
    nsteps = NN - 4
    for row_n in range(rows):
        row = filter[row_n] # Capture individual row
        A1 = row[0]
        A2 = row[1]
        B0 = row[2]
        B1 = row[3]
        B2 = row[4]
        T = 3
        for _ in range(nsteps):
            T = T + 1
            # Apply Filtering Specified by Individual Row
            y[T] = A1*y[T-1] + A2*y[T-2] + B0*x_tmp[T] + B1*x_tmp[T-1] +  B2*x_tmp[T-2]
        # Copy New output into temporary input
        x_tmp = np.copy( y )
    # Copy finalized output into *ytime* for plotting
    ytime = np.copy( x_tmp )
    # Plot Filtered Output
    plt.figure(1)
    plt.plot(x,'k--',label="Input")
    plt.plot(ytime,'k',label="Output")
    plt.title(title)
    plt.grid(which='both')
    if legend: plt.legend()
    plt.show()
    
# Define System Response Plotter function
def sys_response(system,npts=1000,dt=0.01,combine=True,gtitle="",xlim=False,
                stepResponse=True,rampResponse=False,parabolicResponse=False,sv=False):
    """ System Response Plotter Function
    
    Given a transfer function, plots the response against step, ramp, and
    parabolic inputs and plots the error for each of these functions.
    
    Required Arguments
    ------------------
    system:        The Transfer Function; can be provided as the following:
                - 1 (instance of lti)
                - 2 (num, den)
                - 3 (zeros, poles, gain)
                - 4 (A, B, C, D)
    
    Optional Arguments
    ------------------
    npts:                Number of steps to calculate over; default is 1000.
    dt:                    Difference between each data point, default is 0.01.
    combine:            If combination of numerator and denominator is needed.
                        This value should be set to "True" if the parts should be
                        combined to show the complete system with feedback.
                        Default is True.
    gtitle:                Additional string to be added to plot titles;
                        default is "".
    stepResponse:        Plot the step-response and corresponding error;
                        default is True.
    rampResponse:        Plot the ramp-response and corresponding error;
                        default is False.
    parabolicResponse:    Plot the parabolic-response and corresponding error;
                        default is False.
    xlim:                Limit in x-axis for graph plot. Accepts tuple of: (xmin, xmax).
                        Default is False.
    sv:                    Save the figures plotted. Default is False.
    
    Returns
    -------
    NONE:        Generates the plot of the desired responses,
                does not return numerical values.
    """
    # Define Time Axis
    TT = np.arange(0,npts*dt,dt)
    
    # Condition system input to ensure proper execution
    system = sys_condition(system,combine)    
    
    # Allocate space for all outputs
    step = np.zeros(npts)
    ramp = np.zeros(npts)
    parabola = np.zeros(npts)
    errS = np.zeros(npts)
    errR = np.zeros(npts)
    errP = np.zeros(npts)
    
    # Generate Inputs
    for i in range(npts):
        step[i] = 1.0
        ramp[i] = (dt*i)
        parabola[i] = (dt*i)**(2)
    
    # Simulate Response for each input (step, ramp, parabola)
    # All 'x' values are variables that are considered don't-care
    x, y1, x = sig.lsim((system),step,TT)
    x, y2, x = sig.lsim((system),ramp,TT)
    x, y3, x = sig.lsim((system),parabola,TT)
    
    # Calculate error over all points
    for k in range(npts):
        errS[k] = step[k] - y1[k]
        errR[k] = ramp[k] - y2[k]
        errP[k] = parabola[k] - y3[k]
    
    # Plot responses if allowed
    if (stepResponse):
        plt.figure()
        plt.subplot(121)
        plt.title("Step Response "+gtitle)
        plt.plot(TT,y1,'k--', label="Step Response")
        plt.plot(TT,step,'k', label="Step Function")
        plt.grid()
        plt.legend()
        plt.xlabel("Time (seconds)")
        if xlim != False:
            plt.xlim(xlim)
        plt.subplot(122)
        plt.title("Step Response Error "+gtitle)
        plt.plot(TT,errS,'k', label="Error")
        plt.grid()
        plt.legend()
        plt.xlabel("Time (seconds)")
        if xlim != False:
            plt.xlim(xlim)
        plt.subplots_adjust(wspace=0.3)
        if sv:
            plt.savefig("Step Response ("+gtitle+").png")
        plt.show()
    if (rampResponse):
        plt.figure()
        plt.subplot(121)
        plt.title("Ramp Response "+gtitle)
        plt.plot(TT,y2,'k--', label="Ramp Response")
        plt.plot(TT,ramp,'k', label="Ramp Function")
        plt.grid()
        plt.legend()
        plt.xlabel("Time (seconds)")
        if xlim != False:
            plt.xlim(xlim)
        plt.subplot(122)
        plt.title("Ramp Response Error "+gtitle)
        plt.plot(TT,errR,'k', label="Error")
        plt.grid()
        plt.legend()
        plt.xlabel("Time (seconds)")
        if xlim != False:
            plt.xlim(xlim)
        plt.subplots_adjust(wspace=0.3)
        if sv:
            plt.savefig("Ramp Response ("+gtitle+").png")
        plt.show()
    if (parabolicResponse):
        plt.figure()
        plt.subplot(121)
        plt.title("Parabolic Response "+gtitle)
        plt.plot(TT,y3,'k--', label="Parabolic Response")
        plt.plot(TT,parabola,'k', label="Parabolic Function")
        plt.grid()
        plt.legend()
        plt.xlabel("Time (seconds)")
        if xlim != False:
            plt.xlim(xlim)
        plt.subplot(122)
        plt.title("Parabolic Response Error "+gtitle)
        plt.plot(TT,errP,'k', label="Error")
        plt.grid()
        plt.legend()
        plt.xlabel("Time (seconds)")
        if xlim != False:
            plt.xlim(xlim)
        plt.subplots_adjust(wspace=0.3)
        if sv:
            plt.savefig("Parabolic Response ("+gtitle+").png")
        plt.show()
        
def statespace(A,B,x=0,f=0,solution=2,C=False,D=False,npts=9999,NN=10000,dt=0.01,
        xlim=False,ylim=False,gtitle="",ret=False,plot=True,pltfn=False,sv=False):
    """ Plots the state-space simulation of an arbitrary set of matricies.

    Required Arguments:
    --------------------
    A :            Matrix A; if not type=numpy.matrix, converts to numpy.matrix
    B :            Matrix B; if not type=numpy.matrix, converts to numpy.matrix

    Optional Arguments:
    --------------------
    x :            Matrix x; if not type=numpy.matrix, converts to numpy.matrix
    f :            Forcing Function; must be provided as callable function that
                will return any/all forcing function Arguments needed as
                numpy matrix (preferred), numpy array, or tuple.
                Forcing function(s) can be provided as tuple of function
                handles, system will automatically concatenate their output
                to a matrix that can be handled.
    solution:    Determines What Type of Solution Simulation is Required;
                Default of solution is 2
                0=zero-input    ( No Forcing Function )
                1=zero-state    ( No Initial Conditions )
                2=total            ( Both Initial Conditions and Forcing Function )
                3=total, output ( Both ICs and FFs, also plot combined output )
    npts:         Changes the range of simulation; defualt=9999
    NN:            Number of descrete points; default=10,000
    dt:            Delta-t, step-size; default=0.01
    xlim:        Limit in x-axis for graph plot.
    ylim:        Limit in y-axis for graph plot.
    gtitle:        Additional String for Plot Title
    ret:        If true: returns state space terms
    plot:        If true: Plots individual state space terms
    pltfn:        If true: Plots original Forcing Functions

    Figures:
    --------
    Forcing Functions:        The plot of forcing functions, only provided if pltfn is true.
    State Variables:        The plot of state variables, always provided if plot is true.
    Combined Output:        The plot of the combined terms in the output, provided if C and D are not False.

    """

    # Test for NN and npts
    if (npts >= NN):
        print("WARNING: NN must be greater than npts; NN="+str(NN)+"npts="+str(npts))
        print("Autocorrecting npts to be NN-1.")
        npts = NN-1

    # Test for C and D matricies
    mC = str(type(C))
    mD = str(type(D))
    if (((mC==matrix) or (mC==ndarr) or (mC==tuple)) and
        ((mD==matrix) or (mD==ndarr) or (mD==tuple))):
        if (solution!=3):
            print("WARNING: C and D matricies provided, but solution requested "+
                    "does not include combined output.")
    elif (((mC==matrix) or (mC==ndarr) or (mC==tuple)) and
        (mD!=matrix) and (mD!=ndarr) and (mD!=tuple)):
        if (D==False):
            print("WARNING: D matrix not provided; D now assumed to be 0.")
            D = np.matrix('0')
    else:
        C = np.matrix('0')
        D = np.matrix('0')
        if (solution==3):
            print("WARNING: Combined output requested, but no matricies C and D given.")
            print("         Solution being set to: 2 - Complete Simulation")
            solution = 2


    # Create values for input testing
    mA = str(type(A))
    mB = str(type(B))
    mx = str(type(x))
    mC = str(type(C))
    mD = str(type(D))
    if (str(type(f)) == tfun): # if f is a function, test as one
        mF = str(type(f(1))) # f should return: int, float, tuple, np.arr, np.matrix
    elif (str(type(f)) == tuple): # if f is tupple of arguments
        if (str(type(f[0])) == tfun): #if first argument is a function
            c_funcs = c_func_concat(f) # concatinate functions into one
            mF = "MultiFunctions" # label as multiple concatenated functions
        else:
            mF = "NA" # Can't handle function type
    else:
        mF = "NA" # Can't handle function type

    # Test for x input
    if (mx!=matrix) and (mx!=ndarr) and (mx!=tuple):
        if x==0: # No specified initial conditions
            if (mA==matrix): # Use A matrix as reference
                rA = A.shape[0]
                x = np.asmatrix(np.zeros(rA)).T
                mx = str(type(x))
                print("WARNING: No input x (Initial Condition) given.")
                if (solution!=1) and (solution!=3):
                    solution = 1
                    print("\n         Solution type changed to 1: Zero-State.")
            elif (mB==matrix): # Use B matrix as reference
                rB = B.shape[0]
                x = np.asmatrix(np.zeros(rB)).T
                mx = str(type(x))
                print("WARNING: No input x (Initial Condition) given.")
                if (solution!=1) and (solution!=3):
                    solution = 1
                    print("\n         Solution type changed to 1: Zero-State.")
            else:
                raise ValueError("ERROR: No x matrix (Initial Conditions) given,"+
                                "\nNot enough additional information to infer x matrix.")

    # Test for matrix inputs
    if (mA!=matrix) or (mB!=matrix) or (mx!=matrix) or (mC!=matrix) or (mD!=matrix):
        # One or more arguments are not of type numpy.matrix
        # Convert to type matrix
        print("WARNING: Converting one or more input matricies to type: numpy matrix")
        A = np.asmatrix(A)
        B = np.asmatrix(B)
        x = np.asmatrix(x)

    # Gather dimensions of inputs
    rA, cA = A.shape
    rB, cB = B.shape
    rx, cx = x.shape
    rC, cC = C.shape
    rD, cD = D.shape
    rF, cF = 1, 1 # Defualt for a function returning one value

    if (mF==tuple): # If function returns tuple
        print("WARNING: Converting Forcing Function from type: tuple")
        fn = lambda x: tuple_to_matrix(x, f) # Use conversion function
        rF, cF = fn(1).shape # Prepare for further testing
    elif (mF==ndarr): # If function returns numpy array
        print("WARNING: Converting Forcing Function from type: numpy array")
        fn = lambda x: nparr_to_matrix(x, f) # Use conversion function
        rF, cF = fn(1).shape # Prepare for further testing
    elif (mF==tint) or (mF==tfloat) or (mF==tnfloat): # If function returns int or float or numpy float
        fn = f # Pass function handle
    elif (mF==matrix): # If function returns matrix
        fn = f # Pass function handle
        rF, cF = fn(1).shape # Prepare for further testing
    elif (mF=="MultiFunctions"): # There are multiple functions in one argument
        print("WARNING: Casting multiple forcing functions to output type: numpy matrix")
        fn = c_funcs.func_c # Gather function handle from function concatenation class
        rF, cF = fn(1).shape # Prepare for further testing
    elif (mF=="NA"): # Function doesn't meet requirements
        raise ValueError("ERROR: Forcing function does not meet requirements."+
                        "\nFunction doesn't return data type: int, float, numpy.ndarray"+
                        "\n or numpy.matrixlib.defmatrix.matrix. Nor does function "+
                        "\ncontain tuple of function handles. Please review function.")

    # Test for size correlation between matricies
    if (cA != rA): # A isn't nxn matrix
        raise ValueError("ERROR: Matrix 'A' is not NxN matrix.")
    elif (rA != rB): # A and B matricies don't have same number of rows
        if (B.size % rA) == 0: # Elements in B divisible by rows in A
            print("WARNING: Reshaping 'B' matrix to match 'A' matrix.")
            B = np.matrix.reshape(B,(rA,int(B.size/rA))) # Reshape Matrix
        else:
            raise ValueError("ERROR: 'A' matrix dimensions don't match 'B' matrix dimensions.")
    elif (rA != rx): # A and x matricies don't have same number of rows
        if (x.size % rA) == 0: # Elements in x divisible by rows in A
            print("WARNING: Reshaping 'x' matrix to match 'A' matrix.")
            x = np.matrix.reshape(x,(rA,1)) # Reshape Matrix
        else:
            raise ValueError("ERROR: 'A' matrix dimensions don't match 'B' matrix dimensions.")
    elif (cB != rF) or (cF != 1): # Forcing Function matrix doesn't match B matrix
        raise ValueError("ERROR: 'B' matrix dimensions don't match forcing function dimensions.")
    elif (solution==3) and (cC != cA) or (rC != 1): # Number of elements in C don't meet requirements
        raise ValueError("ERROR: 'C' matrix dimensions don't match state-space variable dimensions.")
    elif (solution==3) and ((cD != rF) or (rD != 1)): # Number of elements in D don't meet requirements
        if (cD == rD) and (cD == 1) and (D[0] == 0): # D matrix is set to [0]
            D = np.asmatrix(np.zeros(rF)) # Re-create D to meet requirements
            print("WARNING: Autogenerating 'D' matrix of zeros to match forcing functions.")
        else:
            raise ValueError("ERROR: 'D' matrix dimensions don't match forcing function dimensions.")

    # Test for forcing function
    if (f==0) and (solution!=0):
        print("WARNING: No forcing function provided.\n         "+
                "Solution type changed to 0: Zero-Input")
        solution = 0 # Change to Zero-Input calculation

    # Start by defining Constants
    T = 0
    TT = np.arange(0,(dt*(NN)),dt)
    yout = 0

    # Define list of strings for plot output
    soltype = ["(Zero-Input)","(Zero-State)","(Complete Simulation)","(Complete Sim., Combined Output)"]

    # Create a keyed list of state-space variables
    xtim = {}
    xtim_len = rA # Number of Rows in A matrix
    for n in range(xtim_len):
        key = n #Each key should be the iterative variable
        xtim_init = np.zeros(NN) #Define the initial array
        xtim[key] = xtim_init #Create each xtim

    # Create a keyed list of function outputs
    if (mF!=tint) and (mF!=tfloat):
        fn_arr = {}
        for n in range(rF):
            key = n #Each key should be the iterative variable
            fn_init = np.zeros(NN) #Define the initial array
            fn_arr[key] = fn_init #Create each fn_arr
            fnc = rF
    else:
        fn_arr = np.zeros(NN) #Create the fn_arr
        fnc = 1

    # When asked to find zero-state, set all ICs to zero
    if solution == 1:
        for n in range(xtim_len):
            x[n] = 0 #Set each value to zero

    # Finite-Difference Simulation
    for i in range(0,npts):
        for n in range(xtim_len):
            xtim[n][i] = x[n] #xtim[state-variable][domain] = x[state-variable]
        # Create Forcing Function output

        if fnc > 1: # More than one forcing function
            for n in range(fnc):
                fn_arr[n][i] = np.asarray(fn(T))[n][0]
        else: # only one forcing function
            fn_arr[i] = fn(T)

        if solution == 0: #Zero-input, no added function input
            x = x + dt*A*x
        else: #Zero-state or Total, add function input
            x = x + dt*A*x + dt*B*fn(T)
            if solution==3:
                yout = yout + dt*D*fn(T)

        T = T+dt #Add discrete increment to T

    # Plot Forcing Functions
    if (pltfn):
        fffig = plt.figure("Forcing Functions")
        if fnc > 1:
            for x in range(fnc):
                plt.plot(TT,fn_arr[x],label="f"+str(x+1))
        else:
            plt.plot(TT,fn_arr,label="f1")
        if xlim!=False:
            plt.xlim(xlim)
        if ylim!=False:
            plt.ylim(ylim)
        plt.title("Forcing Functions "+gtitle)
        plt.xlabel("Time (seconds)")
        plt.legend(title="Forcing Functions")
        plt.grid()
        if sv:
            plt.savefig('Simulation Forcing Functions.png')
        if plot:
            plt.show()

    # Plot each state-variable over time
    stvfig = plt.figure("State Variables")
    for x in range(xtim_len):
        plt.plot(TT,xtim[x],label="x"+str(x+1))
    if xlim!=False:
            plt.xlim(xlim)
    if ylim!=False:
        plt.ylim(ylim)
    plt.title("Simulated Output Terms "+soltype[solution]+gtitle)
    plt.xlabel("Time (seconds)")
    plt.legend(title="State Variable")
    plt.grid()
    if sv:
        plt.savefig('Simulation Terms.png')
    if plot:
        plt.show()

    # Plot combined output
    if (solution==3):
        cofig = plt.figure("Combined Output")
        C = np.asarray(C) # convert back to array for operation
        for i in range(cC):
            yout = yout + xtim[i]*C[0][i] # Sum all st-space var mult. by their coeff
        yout = np.asarray(yout) # convert output to array for plotting purposes
        plt.plot(TT,yout[0])
        if xlim!=False:
            plt.xlim(xlim)
        if ylim!=False:
            plt.ylim(ylim)
        plt.title("Combined Output "+gtitle)
        plt.xlabel("Time (seconds)")
        plt.grid()
        if sv:
            plt.savefig('Simulation Combined Output.png')
        if plot:
            plt.show()

    # Return Variables if asked to
    if ret:
        return(TT, xtim)
    
# End of FILTERSIM.PY