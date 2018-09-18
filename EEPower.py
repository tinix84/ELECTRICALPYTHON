###################################################################
#   EEPower.py
#
#   A library of functions, constants and more
#   that are related to Power in Electrical Engineering.
#
#   September 3, 2018
#   August 30, 2018
#
#   Written by Joe Stanley
#
#   Included Constants:
#   - 'A' Operator for Symmetrical Components: a_op
#   -
#
#   Included Functions
#   - Complex Display Function:		vi_cprint
#   - Impedance Conversion:			z_mk
#   - Parallel Impedance Adder: 	P_Zadd
#   - 3-Phase Voltage Converter: 	v_convert
#   - 3-Phase Current Converter: 	i_convert
#   - Power Triangle Function: 		P_triangle
#   - Transformer SC OC Tests:		trans_scoc
#   - Per Unit Base Creator:		pu
#   - Per Unit Base Converter:		pu_conv
#   - Phasor Plot Generator:		phasor_plot
###################################################################

# Import libraries as needed:
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import cmath as c

# Define constants
a_op = c.rect(1,np.radians(120)) # A Operator for Sym. Components
VLLcVLN = c.rect(np.sqrt(3),np.radians(30)) # Conversion Operator
ILcIP = c.rect(np.sqrt(3),np.radians(-30)) # Conversion Operator

###################################################################
#   Define per unit base creator function
#
#   Calculate per unit base given voltage (V) and power (VA)
#   Inputs: v, s, phase, z
#   When calculating 3-phase per unit bases, use 3-phase-S-base
#   and line-to-line-voltage for s and v respectively.
#
#   Returns per unit base of z if z=True, pu base of i if z=False
#   Returns as value for 3-phase by default, can also provide
#   1-phase values.
###################################################################
def pu(v,s,phase=3,z=True):
	if z:
		return(v**2/s)
	elif (phase==3) and not z:
		return(s/(np.sqrt(3)*v))
	else:
		return(v/s)

###################################################################
#   Define per unit converter function
#
#   Converts a [quantity] to a new per-unit base (puB_new) given an
#   old per-unit base (puB_old).
#
#   Returns per-unit value in new base.
###################################################################
def pu_conv(quantity, puB_old, puB_new):
	pu_new = quantity*puB_old/puB_new
	return(pu_new)

###################################################################
#   Define display function
#
#   Print a complex voltage or current in polar form,
#   show angle in degrees.
#
#   Requires voltage or current be provided as complex value.
###################################################################
def vi_cprint(val,unit=False,label=False,printval=True,ret=False):
	mag, ang_r = c.polar(val) #Convert to polar form
	ang = np.degrees(ang_r) #Convert to degrees
  	# Print values (by default)
	if printval and not unit and not label:
		print(mag,"∠",ang,"°")
	elif printval and unit and not label:
		print(mag,"∠",ang,"°"+unit)
	elif printval and unit and label:
		print(label,mag,"∠",ang,"°"+unit)
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
def z_mk(w,C=False,L=False):
	#C Given in ohms, return as Z
	if (C!=False):
		Zc = 1/(1j*w*C)
		return(Zc)
	#L Given in ohms, return as Z
	if (L!=False):
		Zl = 1j*w*L
		return(Zl)

###################################################################
#   Define Parallel Impedance Adder
#
#   Calculate parallel impedance given two inputs.
#
#   Requires both impedance inputs given in phasor form.
###################################################################
def P_Zadd(Z1,Z2):
	Zp = 1/Z1+1/Z2
	return(Zp)

###################################################################
#   Define voltage conversion
#
#   Convert Line-Line voltage to Line-Neutral, or vice-versa.
#
#   Requires that voltage is provided in complex form.
###################################################################
def v_convert(VLL=False,VLN=False):
	#Given VLL, convert to VLN
	if (VLL!=False):
		VLN = VLL/(VLLcVLN)
		return(VLN)
	#Given VLN, convert to VLL
	elif (VLN!=False):
		VLL = VLN*VLLcVLN
		return(VLL)
	#Neither given, error encountered
	else:
		print("ERROR: No value given"+
				"or innapropriate value"+
				"given.")

###################################################################
#   Define current conversion
#
#   Converts Line current to Phase current, or vice-versa.
#
#   Requires that current is provided in complex form.
###################################################################
def i_convert(Iline=False,Iphase=False):
	#Given Iphase, convert to Iline
	if (Iphase!=False):
		Iline = Iphase*ILcIP
		return(Iline)
	#Given Iline, convert to Iphase
	elif (Iline!=False):
		Iphase = Iline/ILcIP
		return(Iphase)
	#Neither given, error encountered
	else:
		print("ERROR: No value given"+
				"or innapropriate value"+
				"given.")

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
def P_triangle(P=False,Q=False,S=False,PF=False,color="red",
			   text="",figure=1,printval=False,ret=False,plot=True):
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
		plt.title("Power Triangle"+str(text))
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

				
# Define Phasor Plot Generator
def phasor_plot(phasor,title="Phasor Diagram",bg="#d5de9c",radius=1.2):
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
		