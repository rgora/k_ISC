#!/usr/bin/env python
from sys import argv 
import numpy as np
from math import *
cmRectoHa=4.5563352812122295e-06
autos=2.4188843265052e-17
masses = {"c" : 1.0000, "o" : 1., "n" : 1., "h" : 1.}
masses = {"c" : 12.0000, "o" : 15.99491, "n" : 14.00307, "h" : 1.00782}
if len(argv)!=5:
	print """Files order: 
		1. Vib_modes of the initial state. 
		2. Vib_spectrum of the initial state (frequencies).
		3. Vib_modes of the final state.
		4. Vib_spectrum of the final state."""
else: 
	print "Extracting vibrational modes of the initial state..."
	state1=open(argv[1],'r')
	coordinates1=open('coord.ini','r')
	lines=state1.readlines()
	coordlines1=coordinates1.readlines()
	vib_modes=[]
	atomic_masses=[]
	r1=[]
	nat=0
	for c in coordlines1:
		if not c.startswith('$'):
			temp_list=c.split()
			atomic_masses.append(sqrt(masses[temp_list[-1]]))
                        atomic_masses.append(sqrt(masses[temp_list[-1]]))
                        atomic_masses.append(sqrt(masses[temp_list[-1]]))
			for j in temp_list[:-1]:
				r1.append(float(j))
			nat+=1
	print "There are",nat,"atoms."
	coord = nat*3
	modes = coord-6
	Mass=np.diag(np.array(atomic_masses))
 	print Mass
	R1=np.array(r1)
	R1.shape=(coord,1)
 	print R1
#
	list1=[]
	for line in lines[1:]:
		temp_list=line.split()
		mass_factor=sqrt(1.0/atomic_masses[int(temp_list[0])-1])
		if  len(list1)!=coord:
			for j in temp_list[2:]:
				list1.append(float(j)*mass_factor)
		elif len(list1)==coord:
#			print len(list1[6:]), "normal modes"
			for k in list1[6:]:
                        	vib_modes.append(k)
#			print len(list1)
			list1=[]
			for j in temp_list[2:]:
				list1.append(float(j)*mass_factor)
#		print line[:-1]
	for k in list1[6:]:
		vib_modes.append(k)
	list1=[]
#
	n_elements=len(vib_modes)
	print "There are", n_elements, "elements."
	Vini=np.array(vib_modes)
  	Vini.shape=(coord,modes)
 	print Vini
	state1.close()
	coordinates1.close()
#
	print "Extracting frequencies for the initial state..."
	initial_spectrum=open(argv[2],'r')
	spec=initial_spectrum.readlines()
	frequencies=[]
	vib_frequencies=[]
	for line in spec[9:-1]:
		freq_list=line.split()
		f=freq_list[2]
#		print f
		vib_frequencies.append(float(f)*cmRectoHa)
	Omega_ini=np.diag(np.array(vib_frequencies))
	initial_spectrum.close()
#	exit()
 	print Omega_ini
#
	print "Extracting vibrational modes of the final state..."
	state2=open(argv[3],'r')
	coordinates2=open('coord.fin','r')
	lines2=state2.readlines()
	coordlines2=coordinates2.readlines()
	vib_modes=[]
	r2=[]
	count=0
	for c in coordlines2:
#		print c
		if not c.startswith('$'):
			temp_list=c.split()
#			print temp_list
			for j in temp_list[:-1]:
				r2.append(float(j))
#	print r2
	R2=np.array(r2)
	R2.shape=(coord,1)
 	print Mass
	print R2
	list1=[]
	for line in lines2[1:]:
		temp_list=line.split()
		mass_factor=sqrt(1.0/atomic_masses[int(temp_list[0])-1])
		if  len(list1)!=coord:
			for j in temp_list[2:]:
				list1.append(float(j)*mass_factor)
		elif len(list1)==coord:
#			print len(list1[6:]), "normal modes"
			for k in list1[6:]:
                        	vib_modes.append(k)
#			print len(list1)
			list1=[]
			for j in temp_list[2:]:
				list1.append(float(j)*mass_factor)
#		print line[:-1]
	for k in list1[6:]:
		vib_modes.append(k)
	list1=[]
#
	print "There are", len(vib_modes), "elements."
	Vfin=np.array(vib_modes)
  	Vfin.shape=(coord,modes)
 	print Vfin
	state2.close()
	VfinT=np.transpose(Vfin)
#	print "Transpose"
#	print VfinT
#
	print "Extracting frequencies for the final state..."
	initial_spectrum=open(argv[4],'r')
	spec=initial_spectrum.readlines()
	frequencies=[]
	vib_frequencies=[]
	for line in spec[9:-1]:
		freq_list=line.split()
		f=freq_list[2]
#		print f
		vib_frequencies.append(float(f)*cmRectoHa)
	Omega_fin=np.diag(np.array(vib_frequencies))
	initial_spectrum.close()
	#exit()
 	print Omega_fin
#
#	Duschinksy Matrix calculation
#
	J=np.dot(VfinT,Vini)
        print J
#
#	Displacement vector
#
	delR=np.subtract(R2,R1)
	D=np.dot(np.dot(VfinT,Mass),delR)
        print delR
        print D
#
#	Remaining matrices
#
	M=np.subtract(np.dot(np.dot(np.transpose(J),np.dot(Omega_fin,Omega_fin)),J),np.dot(Omega_ini,Omega_ini))		
	A=np.dot(np.dot(np.transpose(J),np.dot(Omega_fin,Omega_fin)),D)
	C=0.5*(np.dot(np.dot(np.transpose(D),np.dot(Omega_fin,Omega_fin)),D))
        print M
        print A
        print C
	delE=-.0325027730
#	delE=0.061017902299965954
	soc=48.15*cmRectoHa
#
# 	Evaluating the formula for short-time approximation
#
	sum_a=0
	for i in range(modes):
		a=M[i,i]/Omega_ini[i,i]
		sum_a+=a
	sum_b=0
	for i in range(modes):
		for j in range(modes):
			b=(M[i,j])**2/(Omega_ini[i,i]*Omega_ini[j,j])
			sum_b+=b
	sum_c=0
	for i in range(modes):
		c=(A[i])**2/Omega_ini[i,i]
		sum_c+=c	
	k_short=(soc**2)*sqrt(pi/(0.0625*sum_b + 0.25*sum_c))*exp(-1*((0.25*sum_a + C[0] - delE)**2)/(0.25*sum_b+sum_c))
	print "k_short =",k_short, "[1/a.u.]"
	k_isc=k_short/autos
	print "k_short =", '%e' % k_isc,"[1/s]"
	
