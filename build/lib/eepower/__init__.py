###################################################################
#   EEPOWER.PY
#
#   A library of functions, constants and more
#   that are related to Power in Electrical Engineering.
#
#   February 13, 2019
#   September 3, 2018
#   August 30, 2018
#   February 20, 2019
#
#   Written by Joe Stanley
#
#   Special Thanks To:
#   Paul Ortmann - Idaho Power
#
#   Included Constants:
#   - 'A' Operator for Symmetrical Components: a
#   - Not a Number value (NaN): NAN
#
#   Included Functions
#   - Complex Display Function:     cprint
#   - Impedance Conversion:         phasorz
#   - Parallel Impedance Adder:     parallelz
#   - V/I Line/Phase Converter:     phaseline
#   - Power Triangle Function:      powertriangle
#   - Transformer SC OC Tests:      trans_scoc
#   - Phasor Plot Generator:        phasorplot
#   - Total Harmonic Distortion:    thd
#   - Total Demand Distortion:      tdd
#   - Reactance Calculator:         reactance
#   - Non-Linear PF Calc:           pf_nonlin
#   - Harmonic Limit Calculator:    harmoniclimit
#
#   Additional functions available in sub-modules:
#   - capacitor.py
#   - perunit.py
#   - systemsolution.py
###################################################################
name = "eepower"
ver = "1.1.2"

# Import Submodules
from .capacitor import *
from .perunit import *
from .systemsolution import *

# Import libraries as needed:
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import cmath as c

# Define constants
a = c.rect(1,np.radians(120)) # A Operator for Sym. Components
NAN = float('nan')
VLLcVLN = c.rect(np.sqrt(3),np.radians(30)) # Conversion Operator
ILcIP = c.rect(np.sqrt(3),np.radians(-30)) # Conversion Operator

# Define type constants
matrix = "<class 'numpy.matrixlib.defmatrix.matrix'>"
tuple = "<class 'tuple'>"
ndarr = "<class 'numpy.ndarray'>"
tint = "<class 'int'>"
tfloat = "<class 'float'>"
tfun = "<class 'function'>"

###################################################################
#   Define Reactance Calculator
#
#   Accepts reactance (in ohms) and frequency (in Hertz).
#
#   If imaginary: calculate with j factor (imaginary number)
#
#   Returns capacitance (in Farads) if ohmic value is negative, or
#   inductance (in Henrys) if ohmic value is positive.
###################################################################
def reactance(z,f):
	w = 2*np.pi*f
	if isinstance(z, complex):
		if (z.imag > 0):
			out = z/(w*1j)
		else:
			out = 1/(w*1j*z)
		out = abs(out)
	else:
		if (z > 0):
			out = z/(w)
		else:
			out = 1/(w*z)
		out = abs(out)
	return(out)

###################################################################
#   Define display function
#
#   Print a complex voltage or current in polar form,
#   show angle in degrees.
#
#   Requires voltage or current be provided as complex value.
###################################################################
def cprint(val,unit=False,label=False,printval=True,ret=False,decimals=3):
	mag, ang_r = c.polar(val) #Convert to polar form
	ang = np.degrees(ang_r) #Convert to degrees
    mag = round( mag, decimals ) #Round
    ang = round( ang, decimals ) #Round
  	# Print values (by default)
	if printval and not unit and not label:
		print(mag,"∠",ang,"°")
	elif printval and unit and not label:
		print(mag,"∠",ang,"°",unit)
	elif printval and unit and label:
		print(label,mag,"∠",ang,"°",unit)
	# Return values when requested
	if ret:
		return(mag,ang)

###################################################################
#   Define Impedance Conversion function
#
#   Converts either Capacitance (in Farads) or Inductance (in
#   Henrys) to Impedance value (in ohms).
#
#   Returns C or L as value in Ohms.
###################################################################
def phasorz(f,C=False,L=False,imaginary=True):
	w = 2*np.pi*f
	#C Given in ohms, return as Z
	if (C!=False):
		Z = -1/(w*C)
	#L Given in ohms, return as Z
	if (L!=False):
		Z = w*L
	#If asked for imaginary number
	if (imaginary):
		Z *= 1j
	return(Z)

###################################################################
#   Define Parallel Impedance Adder
#
#   Calculate parallel impedance given a tuple of impedances.
#
#   Requires all impedance inputs given in phasor form.
###################################################################
def parallelz(Z):
	L = len(Z)
	Zp = (1/Z[0]+1/Z[1])**(-1)
	if(L > 2):
		for i in range(2,L):
			Zp = (1/Zp+1/Z[i])**(-1)
	return(Zp)

###################################################################
#   Define Phase/Line Converter
#
#   Convert Line-Line voltage to Line-Neutral, or vice-versa.
#   Converts Line current to Phase current, or vice-versa.
#   Can only convert one voltage or current at a time.
#
#   Input may be provided as absolute value or complex.
#   Output may be specified as complex, but defaults to abs. val.
###################################################################
def phaseline(VLL=False,VLN=False,Iline=False,Iphase=False,complex=False):
output = 0
	#Given VLL, convert to VLN
	if (VLL!=False):
		VLN = VLL/(VLLcVLN)
        output = VLN
	#Given VLN, convert to VLL
	elif (VLN!=False):
		VLL = VLN*VLLcVLN
        output = VLL
	#Given Iphase, convert to Iline
	elif (Iphase!=False):
		Iline = Iphase*ILcIP
        output = Iline
	#Given Iline, convert to Iphase
	elif (Iline!=False):
		Iphase = Iline/ILcIP
        output = Iphase
	#Neither given, error encountered
	else:
		print("ERROR: No value given"+
				"or innapropriate value"+
				"given.")
        return(0)
    #Return as complex only when requested
    if complex:
        return( output )
    return(abs( output ))
###################################################################
#   Define Power Triangle Function
#
#   Draws a Power Triangle given a set of two values from set of:
#   { P, Q, S, PF }. Can also return all values from said set.
#   Plot title can be added to, and plot may be turned off.
#   Color of plot may be changed, and multiple figures may be shown.
#
#   Requires two values from set: { P, Q, S, PF }
#   All values given must be given as absolute value, not complex.
###################################################################
def powertriangle(P=False,Q=False,S=False,PF=False,color="red",
			   text="Power Triangle",figure=1,printval=False,ret=False,plot=True):
	#Given P and Q
	if (P!=False) and (Q!=False):
		S = np.sqrt(P**2+Q**2)
		PF = P/S
		if Q<0:
			PF=-PF
	#Given S and PF
	elif (S!=False) and (PF!=False):
		P = abs(S*PF)
		Q = np.sqrt(S**2-P**2)
		if PF<0:
			Q=-Q
	#Given P and PF
	elif (P!=False) and (PF!=False):
		S = P/PF
		Q = Q = np.sqrt(S**2-P**2)
		if PF<0:
			Q=-Q
	else:
		print("ERROR: Invalid Parameters or too few"+
			 " parameters given to calculate.")
		return(0)

	#Generate Lines
	Plnx = [0,P]
	Plny = [0,0]
	Qlnx = [P,P]
	Qlny = [0,Q]
	Slnx = [0,P]
	Slny = [0,Q]

	#Plot
	if plot:
		plt.figure(figure)
		plt.title(text)
		plt.plot(Plnx,Plny,color=color)
		plt.plot(Qlnx,Qlny,color=color)
		plt.plot(Slnx,Slny,color=color)
		plt.xlabel("Real Power (W)")
		plt.ylabel("Reactive Power (VAR)")
		mx = max(abs(P),abs(Q))

		if P>0:
			plt.xlim(0,mx*1.1)
			x=mx
		else:
			plt.xlim(-mx*1.1,0)
			x=-mx
		if Q>0:
			plt.ylim(0,mx*1.1)
			y=mx
		else:
			plt.ylim(-mx*1.1,0)
			y=-mx
		if PF > 0:
			PFtext = " Lagging"
		else:
			PFtext = " Leading"
		text = "P:   "+str(P)+" W\n"
		text = text+"Q:   "+str(Q)+" VAR\n"
		text = text+"S:   "+str(S)+" VA\n"
		text = text+"PF:  "+str(abs(PF))+PFtext+"\n"
		text = text+"ΘPF: "+str(np.degrees(np.arccos(PF)))+"°"+PFtext
		# Print all values if asked to
		if printval:
			 plt.text(x/20,y*4/5,text,color=color)
		plt.show()

	# Return when requested
	if ret:
		return(P,Q,S,PF)

###################################################################
#   Define Transformer Short-Circuit/Open-Circuit Function
#
#   Calculates Req and Xeq, or Rc and Xm, or both sets given three
#   values from a specific set of inputs { Poc, Voc, Ioc,  Psc,
#   Vsc, Isc }.
#
#   Requires one or both of two sets: { Poc, Voc, Ioc }, or
#   { Psc, Vsc, Isc }.
#   Will return: { Rc, Xm } given first set, { Req, Xeq } given
#   second set.
#   All values given must be given as absolute value, not complex.
#   All values returned are given with respect to high-side/primary
###################################################################
def trans_scoc(Poc=False,Voc=False,Ioc=False,Psc=False,Vsc=False,
				Isc=False):
	SC = False
	OC = False
	# Given Open-Circuit Values
	if (Poc!=False) and (Voc!=False) and (Ioc!=False):
		PF = Poc/(Voc*Ioc)
		Y = c.rect(Ioc/Voc,-np.arccos(PF))
		Rc = 1/Y.real
		Xm = -1/Y.imag
		OC = True
	# Given Short-Circuit Values
	if (Psc!=False) and (Vsc!=False) and (Isc!=False):
		PF = Psc/(Vsc*Isc)
		Zeq = c.rect(Vsc/Isc,np.arccos(PF))
		Req = Zeq.real
		Xeq = Zeq.imag
		SC = True
	# Return All if Found
	if OC and SC:
		return(Req,Xeq,Rc,Xm)
	elif OC:
		return(Rc,Xm)
	elif SC:
		return(Req,Xeq)
	else:
		print("An Error Was Encountered.\n"+
				"Not enough arguments were provided.")

###################################################################
#   Define Phasor Plot Generator
#
#   Plots a phasor-diagram with angles in degrees for a number of
#   Phasors. Phasors must be passed as a set of complex numbers.
#   (e.g. [ m+ja, m+ja, m+ja, ... , m+ja ] ) No more than 12
#   Phasors are allowed to be plotted at one time.
###################################################################
def phasorplot(phasor,title="Phasor Diagram",bg="#d5de9c",radius=1.2):
	numphs = len(phasor)
	
	colors = ["#FF0000","#800000","#FFFF00","#808000","#00ff00","#008000",
			"#00ffff","#008080","#0000ff","#000080","#ff00ff","#800080"]
	
	if numphs > 12:
		raise ValueError("ERROR: No more than 12 phasors allowed.")
	
	# Force square figure and square axes looks better for polar, IMO
	width, height = matplotlib.rcParams['figure.figsize']
	size = min(width, height)
	# Make a square figure
	fig = plt.figure(figsize=(size, size))
	ax = fig.add_axes([0.1, 0.1, 0.8, 0.8], polar=True, facecolor=bg)
	ax.set_rmax(radius)
	plt.grid(True)
	
	# Plot the diagram
	plt.title(title+"\n")
	for i in range(numphs):
		mag, ang_r = c.polar(phasor[i])
		plt.arrow(0,0,ang_r,mag,color=colors[i])
	plt.show()

###################################################################
#   Define Total Demand Distortion function
#
#   Returns the Total demand distortion given an array of the
#   distortion factors.
#
#   Harmonic array should contain the fundamental frequency h1
#   IL: Peak Demand load current (RMS) at fundamental frequency
###################################################################
def tdd(harmonics,IL):
	# Sum all terms of 2*fundamental and up
	sum = 0
	for h in range(1,len(harmonics)):
		sum += harmonics[h]**2
	
	# Take square-root of sum and divide by IL
	TDD = np.sqrt(sum)/IL
	return(TDD)

###################################################################
#   Define Total Harmonic Distortion function
#
#   Returns the Total harmonic distortion given an array of the
#   distortion factors.
#
#   Harmonic array should contain the fundamental frequency h1
#   PFdist: the distorted power factor, can be used to find thd
###################################################################
def thd(harmonics=False,PFdist=False):
	if(PFdist != False):
		# Use PFdistortion equation to find THD
		THD = np.sqrt(1/(PFdist**2)-1)
	else:
		# Sum all terms of 2*fundamental and up
		sum = 0
		for h in range(1,len(harmonics)):
			sum += harmonics[h]**2
		# Take Square Root of Sum
		sum = np.sqrt(sum)
		# Divide by magnitude of fundamental frequency
		THD = sum/harmonics[0]
	return(THD)

###################################################################
#   Define Distorted Power Factor function
#
#   Returns the distorted power factor value given I1 and IRMS or
#   the set of harmonic current magnitudes.
#
#   Ih array should contain the fundamental frequency h1
###################################################################
def pf_dist(I1=False,IRMS=False,Ih=False):
	if (I1 != False and IRMS != False):
		# Find PFdist by using fundamental and RMS current
		PFdist = I1/IRMS
	else:
		# Find PFdist by using THD
		THD = thd(Ih) # Calculate THD
		PFdist = 1/np.sqrt(1+THD**2)
	
	return(PFdist)

###################################################################
#   Define Non-Linear Power Factor Calculator
#
#   Returns the the unknown variable of the set { PFtrue, PFdist,
#   PFdisp }.
#
#   Requires any two of the three inputs
###################################################################
def pf_nonlin(PFtrue=False,PFdist=False,PFdisp=False):
	if(PFtrue!=False and PFdist!=False and PFdisp!=False):
		raise ValueError("ERROR: Too many constraints, no solution.") 
	elif ( PFdist!=False and PFdisp!=False ):
		return( PFdist * PFdisp )
	elif ( PFtrue!=False and PFdisp!=False ):
		return( PFtrue / PFdisp )
	elif ( PFtrue!=False and PFdist!=False ):
		return( PFtrue / PFdist )
	else:
		raise ValueError("ERROR: Function requires at least two arguments.") 

###################################################################
#   Define Harmonic Current Limit function
#
#   Returns the limits of harmonic currents given the load
#   characteristics (Short-Circuit Current, Peak Demand Current).
#   Prints results when printout=True.
#
#   By default, prints the results, does not return values as an
#   numpy array. Returns this array when ret=True.
#
#   Compares to measured harmonic currents if Ih is provided.
#
#   N is the maximum harmonic term.
###################################################################
def harmoniclimit(Isc,IL,N=0,Ih=0,printout=True,ret=False):
	percent = 1/100 # Use for scaling
	Ir = Isc/IL # compute ratio
	if(Ir < 20):
		# Generate Harmonic Factors
		f1o = 4.0*percent
		f2o = 2.0*percent
		f3o = 1.5*percent
		f4o = 0.6*percent
		f5o = 0.3*percent
		tddL = 5.0*percent
		f1e = f1o * 25*percent
		f2e = f2o * 25*percent
		f3e = f3o * 25*percent
		f4e = f4o * 25*percent
		f5e = f5o * 25*percent
		
	elif(20 <= Ir and Ir < 50):
		# Generate Harmonic Factors
		f1o = 7.0*percent
		f2o = 3.5*percent
		f3o = 2.5*percent
		f4o = 1.0*percent
		f5o = 0.5*percent
		tddL = 8.0*percent
		f1e = f1o * 25*percent
		f2e = f2o * 25*percent
		f3e = f3o * 25*percent
		f4e = f4o * 25*percent
		f5e = f5o * 25*percent
		
	elif(50 <= Ir and Ir < 100):
		# Generate Harmonic Factors
		f1o = 10.0*percent
		f2o = 4.5*percent
		f3o = 4.0*percent
		f4o = 1.5*percent
		f5o = 0.7*percent
		tddL = 12.0*percent
		f1e = f1o * 25*percent
		f2e = f2o * 25*percent
		f3e = f3o * 25*percent
		f4e = f4o * 25*percent
		f5e = f5o * 25*percent
		
	elif(100 <= Ir and Ir < 1000):
		# Generate Harmonic Factors
		f1o = 12.0*percent
		f2o = 5.5*percent
		f3o = 5.0*percent
		f4o = 2.0*percent
		f5o = 1.0*percent
		tddL = 15.0*percent
		f1e = f1o * 25*percent
		f2e = f2o * 25*percent
		f3e = f3o * 25*percent
		f4e = f4o * 25*percent
		f5e = f5o * 25*percent
		
	else:
		# Generate Harmonic Factors
		f1o = 15.0*percent
		f2o = 7.0*percent
		f3o = 6.0*percent
		f4o = 2.5*percent
		f5o = 1.4*percent
		tddL = 20.0*percent
		f1e = f1o * 25*percent
		f2e = f2o * 25*percent
		f3e = f3o * 25*percent
		f4e = f4o * 25*percent
		f5e = f5o * 25*percent
	
	# Create empty array to return
	retArr = np.zeros(51)
	
	# Print out values
	if(printout):
		print("IEEE 519-2014 Distorted Current Limits:\n"+
				"---------------------------------------")
		for i in range(3, 50 + 1):
			if(3<=i and i<11):
				if(i%2): # Odd term
					retArr[i] = f1o*IL
				else: # Even term
					retArr[i] = f1e*IL
			elif(11<=i and i<17):
				if(i%2): # Odd term
					retArr[i] = f2o*IL
				else: # Even term
					retArr[i] = f2e*IL
			elif(17<=i and i<23):
				if(i%2): # Odd term
					retArr[i] = f3o*IL
				else: # Even term
					retArr[i] = f3e*IL
			elif(23<=i and i<35):
				if(i%2): # Odd term
					retArr[i] = f4o*IL
				else: # Even term
					retArr[i] = f4e*IL
			elif(35<=i and i<=50):
				if(i%2): # Odd term
					retArr[i] = f5o*IL
				else: # Even term
					retArr[i] = f5e*IL
			else:
				print("Internal Error Encountered!")
			print("Limit of "+str(i)+"th harmonic:", retArr[i],"A")
			if(N!=0 and N<=i):
				break
		print("---------------------------------------")
	
	# Comparison requested
	if(str(type(Ih)) == ndarr):
		maxr = min(len(retArr), len(Ih)+1)
		for k in range(3, maxr):
			if(retArr[k] < Ih[k-1]):
				print("Limit surpassed for "+str(k)+"th Harmonic term.")
	
	# Return values
	if(ret):
		return(retArr)
	