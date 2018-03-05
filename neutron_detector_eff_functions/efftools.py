#!/usr/bin/env python

from pylab import exp, zeros
import numpy as np


def efficiency2particles(d,R1,R2,sigma):
	"""This script is the implementation of Equations from 3.11 to 3.16 from page 63 on Francesco's PhD thesis (download from twiki). The output is efficiency for a BS layer and a T layer,
	varibles must be scalar, d is the thickness of the coating, R1 and R2 effective ranges of the two partciles emitted back-to-back (as defined at page 60 according to an energy threshold),
	sigma the macroscopic cross section. Note that if you want ot consider an inclined layer you just need to feed into sigma the effective sigma (Equation 3.9 page 62) which is sigma/sin(theta),
	with theta the angle between neutron beam and the layer (defined in Fig. 3.1 page 59)."""

	# assert isinstance(R1,float) or isinstance(R1,int)
	# assert isinstance(R2,float) or isinstance(R2,int)

	d=float(d)
	R1=float(R1)
	R2=float(R2)
	sigma=float(sigma)
	if R1<R2:
		R1,R2=R2,R1

	assert R2<=R1

	eff=zeros((2,1))

	if d!=0:
		if d<R2:
			# d<R2<R1 (3.11)
			# Back scattering
			eff[0,0]=(1.-1./(2.*sigma*R1)-1./(2.*sigma*R2))*(1.-exp(-sigma*d))+(1./(2.*R1)+1./(2.*R2))*d*exp(-sigma*d)
			# Transmission
			eff[1,0]=(1.+1./(2.*sigma*R1)+1./(2.*sigma*R2))*(1.-exp(-sigma*d))-(1./(2.*R1)+1./(2.*R2))*d
		elif d>=R2 and d<R1:
			# R2<d<R1
			# Back scattering
			eff[0,0] = (1.-1./(2.*R1*sigma)-1./(2.*R2*sigma))+((exp(-sigma*R2))/(2.*sigma*R2))-(1.-1./(sigma*R1)-d/R1)*((exp(-sigma*d))/2.)
			# Transmission
			eff[1,0] = ((exp((R2-d)*sigma))/(2.*R2*sigma))-(1.+1./(2.*R1*sigma)+1/(2.*R2*sigma))*(exp(-sigma*d))+ (1./2.)*(1.+1./(sigma*R1)-d/R1)

		else:
			# R2<R1<d
			# Back scattering
			eff[0,0] = (1.-1./(2.*R1*sigma)-1./(2.*R2*sigma))+ ((exp(-sigma*R2))/(2.*sigma*R2))+((exp(-sigma*R1))/(2.*sigma*R1))
			# Transmission
			eff[1,0] = (exp(-sigma*d))*(-1.-1./(2.*R1*sigma)-1./(2.*R2*sigma)+((exp(sigma*R2))/(2.*sigma*R2))+((exp(sigma*R1))/(2.*sigma*R1)))
	return eff


def efficiency4boron(d, Ralpha94, RLi94, Ralpha06, RLi06, sigma):
	"""Francesco's PhD thesis (download from twiki) is the reference for these scripts. The output is efficiency for a blade (2 layers BS+T in this order), a BS layer and a T layer of B4C,
	all varibles must be scalar, d is the thickness of the coating, R are the 4 effective ranges of the two partciles emitted back-to-back (as defined at page 60 according to an energy threshold),
	sigma the macroscopic cross section. Note that if you want ot consider an inclined layer you just need to feed into sigma the effective sigma (Equation 3.9 page 62) which is sigma/sin(theta),
	with theta the angle between neutron beam and the layer (defined in Fig. 3.1 page 59)"""

	d=float(d)
	Ralpha94=float(Ralpha94)
	RLi94=float(RLi94)
	Ralpha06=float(Ralpha06)
	RLi06=float(RLi06)
	sigma=float(sigma)

	temp1=efficiency2particles(d,Ralpha94,RLi94,sigma)
	temp2=efficiency2particles(d,Ralpha06,RLi06,sigma)

	eff=zeros((3,1))

	eff[1,0]= 0.94*temp1[0,0]+0.06*temp2[0,0]          #backscatt
	eff[2,0]= 0.94*temp1[1,0]+0.06*temp2[1,0]          #transmission
	eff[0,0]= eff[1,0] + (exp(-d*sigma))*eff[2,0]      #total eff blade

	return eff


def mg_same_thick(sigma_eq, ranges, thickness, nb):
	"""calculates efficiency for configuration of multi grid with blades of same thickness, doesn't considerate aluminium.

	Args:
		sigma_eq (int):
		ranges (array[4]):  array of ranges of particles
		nb (int): number of blades
		thickness (float): thickness of b10

	Returns:
		eff: total efficiency result



	..  Original source in Matlab: https://bitbucket.org/europeanspallationsource/dg_matlabborontools/src/bcbac538ad10d074c5150a228847efc2e0269e0d/B10tools/rangesMAT.m?at=default&fileviewer=file-view-default

	"""
	eff = []
	temp = efficiency4boron(thickness, ranges[0],  ranges[1], ranges[2], ranges[3], sigma_eq)
	expi = exp(-2*sigma_eq*thickness)
	eff.append((temp[0][0]*(1-expi**nb)/(1-expi)))
	#TODO add aluminium consideration
	return eff


def data_samethick_vs_thickandnb(sigma_eq, ranges, nb, window):
	"""calculates efficiency for configuration of multi grid with blades of same thickness and plots it in a given window, doesn't considerate aluminium.

	Args:
		sigma_eq (int):
		ranges (array[4]):  array of ranges of particles
		nb (list): number of blades
		window (qt): window with figure to plot in



	..  Original source in Matlab: https://bitbucket.org/europeanspallationsource/dg_matlabborontools/src/bcbac538ad10d074c5150a228847efc2e0269e0d/MultiGrid_Optimization/MG1_Calc4monoch_sameThickBlades_VS_ThickAndNb.m?at=default&fileviewer=file-view-default

	"""
	thicklist = np.arange(0.0011, 5, 0.05)
	eff = []
	efftotal = []
	c = 0
	nb.append(5)
	nb.append(15)
	nb.append(20)
	#nb = list(set(nb) - set(nb))
	#TODO add aluminium consideration
	for b in nb:
		for n in thicklist:
			eff.append(mg_same_thick(sigma_eq, ranges, n, b))
		efftotal.append(eff)
		eff = []

	# random data
	return efftotal
	# create an axis
	ax = window.figure.add_subplot(111)
	ax.set_xlabel('d (um)')
	ax.set_ylabel('Efficiency')
	# discards the old graph
	#ax.hold(False)
	ax.grid()
	# plot data
	c=0
	ax.set_color_cycle(['red', 'blue', 'black', 'green'])
	for b in nb:
		if c==0:
			ax.plot(thicklist, efftotal[c], '-', label='Nb='+str(b), linewidth=3)
		else:
			ax.plot(thicklist, efftotal[c], '-', label='Nb=' + str(b))
		ax.plot([thicklist[efftotal[c].index(max(efftotal[c])[0])], thicklist[efftotal[c].index(max(efftotal[c])[0])]], [0, max(efftotal[c])[0]], '--',color='black')
		ax.plot([0, thicklist[efftotal[c].index(max(efftotal[c])[0])]], [max(efftotal[c]), max(efftotal[c])], '--', color='black')
		c +=1
	ax.legend()
	# refresh canvas

	window.canvas.draw()

def data_samethick_vs_thickandnb_depth(sigma_eq, ranges, blades):
	"""calculates efficiency for configuration of multi grid with blades of same thickness

	Args:
		sigma_eq (int):
		ranges (array[4]):  array of ranges of particles
		nb (int): number of blades
		blades (list): list of blades

	Returns:
		eff[0]: total efficiency result
		eff[1]: list of efficiency per blade without substrate
		eff[2]:  list of efficiency per blade with substrate

	"""
	if blades[0].substrate == 0:
		#transparent substrate
		delta = 1
	else:
		delta = blades[0].substrate
	eff1blade = []
	efftotal = []
	for b in blades:
		eff1blade.append(efficparam(b.backscatter, sigma_eq, ranges, delta))
	cumthick = 0
	c = 0
	#TOFIX
	for b in blades:
		expi = exp(-2*sigma_eq*cumthick)
		# no substrate ,  with substrate
		efftotal.append([eff1blade[c][3]*expi, eff1blade[c][4]*expi*(delta**c)])
		cumthick = cumthick+b.backscatter
		c = +1
	sinsub = 0
	sub = 0
	for e in efftotal:
		sinsub = sinsub + e[0]
		sub = sub + e[1]
	return efftotal, sinsub, sub

def metadata_samethick_vs_thickandnb(sigma_eq, ranges, nb):
	"""gets metadata for plotting effVSthick of multi grid with blades of same thickness, doesn't considerate aluminium.

	Args:
		sigma_eq (list):
		ranges (array[4]):  array of ranges of particles
		nb (int): number of blades
	returns:
		thicklist (list):  x axis
		eff  (list): y axis



	..  Original source in Matlab: https://bitbucket.org/europeanspallationsource/dg_matlabborontools/src/bcbac538ad10d074c5150a228847efc2e0269e0d/MultiGrid_Optimization/MG1_Calc4monoch_sameThickBlades_VS_ThickAndNb.m?at=default&fileviewer=file-view-default

	"""
	#TODO add aluminium consideration
	thicklist = np.arange(0.0011, 5, 0.05)
	eff = []
	for n in thicklist:
		efftemp = 0
		for s in sigma_eq:
			efftemp = efftemp + mg_same_thick(s[0], ranges, n, nb)[0] * s[1]
		eff.append(efftemp)
	return thicklist, eff,

def metadata_diffthick_vs_thickandnb(sigma_eq, ranges, nb):
	"""gets metadata for plotting effVSthick of multi grid with blades of same thickness, doesn't considerate aluminium.

	Args:
		sigma_eq (int):
		ranges (array[4]):  array of ranges of particles
		nb (int): number of blades
	returns:
		thicklist (list):  x axis
		eff  (list): y axis



	..  Original source in Matlab: https://bitbucket.org/europeanspallationsource/dg_matlabborontools/src/bcbac538ad10d074c5150a228847efc2e0269e0d/MultiGrid_Optimization/MG1_Calc4monoch_sameThickBlades_VS_ThickAndNb.m?at=default&fileviewer=file-view-default

	"""
	thicklist = np.arange(0.0011, 5, 0.05)
	eff = []
	#TODO add aluminium consideration
	for n in thicklist:
		eff.append(mg_same_thick(sigma_eq, ranges, n, nb)[0])
	return thicklist, eff,

#TODO Different thicknesses
def metadata_samethick_vs_wave(sigmaeq, blades, ranges, nb):
	eff = []
	for sigma in sigmaeq:
		# TODO add aluminium consideration
		eff.append(mg_same_thick(sigma, ranges, blades[0].backscatter, nb)[0])
	return eff

def metadata_diffthick_vs_wave(sigmaeq, blades, ranges, nb):
	eff = []
	#intermediate variable for calculating eff sumatory
	effd = 0
	c = 0
	thickness = []
	for b in blades:
		thickness.append(b.backscatter)
	for sigma in sigmaeq:
			effd = mgeff_depth_profile(thickness, ranges, sigma[0], 1)[1]
			eff.append(effd*100)
	return eff


def metadata_singleLayer_vs_wave(sigmaeq, thickness, ranges, nb):
	eff = [[],[],[]]
	for sigma in sigmaeq:
		result = efficiency4boron(thickness, ranges[0], ranges[1], ranges[2], ranges[3], sigma[0])
		print(result)
		eff[0].append(result[1][0]*100)#bs
		eff[1].append(result[2][0]*100)#trans
		eff[2].append(result[0][0]*100)#total
	return eff

def metadata_samethick_vs_thickandnb_single(sigma_eq, ranges, nb):
	"""gets metadata for plotting effVSthick of single layer, doesn't considerate aluminium.

	Args:
		sigma_eq (list):
		ranges (array[4]):  array of ranges of particles
		nb (int): number of blades
	returns:
		thicklist (list):  x axis
		eff  (list): y axis



	..  Original source in Matlab: https://bitbucket.org/europeanspallationsource/dg_matlabborontools/src/bcbac538ad10d074c5150a228847efc2e0269e0d/MultiGrid_Optimization/MG1_Calc4monoch_sameThickBlades_VS_ThickAndNb.m?at=default&fileviewer=file-view-default

	"""
	thicklist = np.arange(0.0011, 10, 0.05)
	effback = []
	efftrans = []
	efftotal = []
	for n in thicklist:
		c = 0
		efftemp = [0,0,0]
		for s in sigma_eq:
			rtemp = efficiency4boron(n, ranges[0], ranges[1], ranges[2], ranges[3], s[0])
			efftemp[0] = efftemp[0] + rtemp[1][0]*s[1]
			efftemp[1] = efftemp[1] + rtemp[2][0]*s[1]
			efftemp[2] = efftemp[2] + rtemp[0][0]*s[1]
			c += 1
		effback.append(efftemp[0])
		efftrans.append(efftemp[1])
		efftotal.append(efftemp[2])
	return thicklist, effback, efftrans, efftotal


def mgeff_depth_profile(thickness, ranges, sigma, varargin):
	if varargin is None:
		# transparent substrate
		delta = 1
	else:
		delta = varargin
	cumthick = 0
	eff1blade = []
	for i, t in enumerate(thickness):
		temp = efficparam(thickness[i], sigma, ranges, varargin)
		eff1blade.append([temp[3], temp[4]])
	eff =[]
	efftotal = 0
	for i, t in enumerate(thickness):
		expi = exp(-2*sigma*cumthick)
		eff.append([eff1blade[i][0]*expi, eff1blade[i][1]*(delta**i)*expi])
		efftotal = efftotal + eff[i][0]
		cumthick = cumthick+t
	return eff, efftotal


def efficparam(thickness,sigma_eq,ranges,varargin):
	"""calculates efficiency of a double layer blade (back and transmission) with the same thickness of bs and t

    Args:
        sigma_eq (int):
        ranges (array[4]):  array of ranges of particles
        nb (int): number of blades
    returns:
        efficiency (list):	[0]	eff back-scattering only
							[1] eff transmission only
							[2] eff trasmission after crossing back-scattering
							[3] eff double layer (bscatt first and then transmission)
							[4] eff trasmission after crossing back-scattering and substrate

..  Original source in Matlab: https://bitbucket.org/europeanspallationsource/dg_matlabborontools/src/bcbac538ad10d074c5150a228847efc2e0269e0d/MultiGrid_Optimization/MG1_Calc4monoch_sameThickBlades_VS_ThickAndNb.m?at=default&fileviewer=file-view-default

    """
	if varargin is None:
		#transparent substrate
		delta = 1
	else:
		delta = varargin
	c1 = efficiency2particles(thickness,ranges[0],ranges[1],sigma_eq)
	d1 = efficiency2particles(thickness,ranges[2],ranges[3],sigma_eq)
	efficiency = []
	# backscatt
	efficiency.append((0.94 * c1[0] + 0.06 * d1[0])[0])
	# trasmission
	efficiency.append((0.94 * c1[1] + 0.06 * d1[1])[0])
	#TODO
	# eff trasmission after crossing back-scattering
	efficiency.append((exp(-thickness * sigma_eq)) * efficiency[1])
	# BS + T
	efficiency.append(efficiency[0] + efficiency[2])
	# BS + T + substrate
	efficiency.append(efficiency[0] + delta * efficiency[2])
	return efficiency

#def sigma(wavelength,theta,massdensity,composition)
#"""Francesco's PhD thesis (download from twiki) is the reference for these scripts. The #output is the effective  macroscopic cross section  """
# sigmadqdx = wiefn
#return sigmanum

if __name__=='__main__':
	#t efficiency2particles(3,2,1,1)
    #print efficiency4boron(1, 2.8, 1.1, 3, 1.3, 3)
	#print 'ciao'
    print ()