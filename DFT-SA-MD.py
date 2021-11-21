# module used

import random
import numpy as np
import os

# read original Cu eam/fs file for simulated annealing

f0 = file("Cu.eam.fs","r")
fs = f0.readlines()
f0.close()

# DFT data saved for simulated annealing

f1 = file("F100DFT","r")
f100 = f1.readlines()
f1.close()

F100DFT = []

i = 0
while i < len(f100) :
	force = f100[i].split()
	fx = float(force[0])
	fy = float(force[1])
	fz = float(force[2])
	F100DFT.append([fx, fy, fz])
	i = i + 1

f2 = file("E100DFT","r")
e100 = f2.readlines()
f2.close()

E100DFT = []

i = 0
while i < len(e100) :
	energy = float(e100[i])
	E100DFT.append(energy)
	i = i + 1

f3 = file("F111DFT","r")
f111 = f3.readlines()
f3.close()

F111DFT = []

i = 0
while i < len(f111) :
	force = f111[i].split()
	fx = float(force[0])
	fy = float(force[1])
	fz = float(force[2])
	F111DFT.append([fx, fy, fz])
	i = i + 1

f4 = file("E111DFT","r")
e111 = f4.readlines()
f4.close()

E111DFT = []

i = 0
while i < len(e111) :
	energy = float(e111[i])
	E111DFT.append(energy)
	i = i + 1

f5 = file("P100111DFT","r")
p100111 = f5.readlines()
f5.close()

P100111DFT = []

i = 0
while i < len(p100111) :
        energy = float(p100111[i])
        P100111DFT.append(energy)
        i = i + 1
	
# read original LAMMPS input files for simulated annealing

r1 = file("input.100","r")
input100 = r1.readlines()
r1.close()

r2 = file("input.111","r")
input111 = r2.readlines()
r2.close()

# original parameters for Morse and Grimme potentials

s6 = 0.60

fN = 4.24
fH = 1.29

D0HC = 0.00176
alphaHC = 1.87
r0HC = 2.51

D0HA2 = 0.00157
alphaHA2 = 2.74
r0HA2 = 2.82

D0CT2 = 0.00137
alphaCT2 = 1.8
r0CT2 = 4.27

D0NH2 = 0.00418
alphaNH2 = 4.13
r0NH2 = 2.31

# potential energy for bare Cu surface

pe100 = -3620.9882289838
pe111 = -3417.4543937389

# number of molecules on Cu

N100 = 72
N111 = 56

# denominators of the cost function from DFT results

FcostF100D = 7.109291647
FcostE100D = 2.454783391

FcostF111D = 4.152302339
FcostE111D = 2.069379808

FcostP100111D = 0.023179279

# function to choose new parameters

def newparameter(x, range) :
	r01 = random.uniform(0.0, 1.0)  # choose a random floating point number between 0 and 1
	return (x - range * x + 2 * range * x * r01)

range = 0.01

# weighting factor for energy in Fcost

w = 0.041667

# simulated annealing temperature

T = 0

N = 0  # number of loops starting from 0

C = 15  # number of configurations on each surface

d = 0.1  # distance in z-direction between each configuration

Fcost = []
Nkeep = []
Fkeep = []
Fi = 10000  # the very first Fcost kept
Fkeep.append(Fi)

# ---beginning of simulated annealing loop---

n = 0  # initial number for Fcost
m = 0  # initial number for Fkeep 
while n <= N :	
    # create new Cu eam/fs file with fN and fH

	i = 50007
	j = 90008
	k = 10006
	while i < 60007 and j < 100008 and k < 20006 :
		fs[i] = float(fs[k]) * fN
		fs[j] = float(fs[k]) * fH
		fs[i] = "%.25s \n" % (fs[i])
		fs[j] = "%.25s \n" % (fs[j])
		i = i + 1
		j = j + 1
		k = k + 1

	w0 = file("Cu.fNH.eam.fs","w")
	i = 0
	while i < len(fs) :
		w0.write(str(fs[i]))
		i = i + 1
	w0.close()
	
	# create new LAMMPS input file for 100 and 111
	
	i = 6
	input100[i] = "pair_style hybrid/overlay eam/fs morse 6.0 zhou 12.0 %.4s 20.0 \n" % (s6)
	input111[i] = input100[i]
	i = 11
	input100[i] = "pair_coeff 1 7 morse %.7s %.4s %.4s # HC \n" % (D0HC, alphaHC, r0HC)
	input111[i] = input100[i]
	i = 12
	input100[i] = "pair_coeff 2 7 morse %.7s %.4s %.4s # HA2 \n" % (D0HA2, alphaHA2, r0HA2)
	input111[i] = input100[i]
	i = 14
	input100[i] = "pair_coeff 4 7 morse %.7s %.4s %.4s # CT2 \n" % (D0CT2, alphaCT2, r0CT2)
	input111[i] = input100[i]
	i = 16
	input100[i] = "pair_coeff 6 7 morse %.7s %.4s %.4s # NH2 \n" % (D0NH2, alphaNH2, r0NH2)
	input111[i] = input100[i]

	w1 = file("newinput.100","w")
	i = 0		
	while i < len(input100) : 
		w1.write(str(input100[i]))
		i = i + 1
	w1.close()
	
	w2 = file("newinput.111","w")
	i = 0		
	while i < len(input111) : 
		w2.write(str(input111[i]))
		i = i + 1
	w2.close()
	
	# ---loops executing LAMMPS for "C" HDA configurations on Cu(100) and Cu(111)---

	F100FF = []
	E100FF = []
	F111FF = []
	E111FF = []

	# read original data files for simulated annealing
	
	r3 = file("data.100-0.6N","r")
	data100 = r3.readlines()
	r3.close()

	r4 = file("data.111-0.6N","r")
	data111 = r4.readlines()
	r4.close()
	
	c = 0
	while c < C :
		# create new LAMMPS data file for 100 and 111

		i = 30
		while i < 3774 :
			data1 = data100[i].split()
			id = data1[0]
			mol = data1[1]
			type = data1[2]
			charge = data1[3]
			x1 = data1[4]
			y1 = data1[5]
			z1 = float(data1[6]) + d
			hash = data1[7]
			atom = data1[8]
			data100[i] = "%.4s %.2s %.1s %.5s %.10s %.10s %.10s %.1s %.3s \n" % (id, mol, type, charge, x1, y1, z1, hash, atom)
			i = i + 1
		
		i = 30
		while i < 2942 :
			data2 = data111[i].split()
			id = data2[0]
			mol = data2[1]
			type = data2[2]
			charge = data2[3]
			x2 = data2[4]
			y2 = data2[5]
			z2 = float(data2[6]) + d
			hash = data2[7]
			atom = data2[8]
			data111[i] = "%.4s %.2s %.1s %.5s %.10s %.10s %.10s %.1s %.3s \n" % (id, mol, type, charge, x2, y2, z2, hash, atom)
			i = i + 1

		w3 = file("data.100","w")
		i = 0		
		while i < len(data100) : 
			w3.write(str(data100[i]))
			i = i + 1
		w3.close()
		
		w4 = file("data.111","w")
		i = 0		
		while i < len(data111) : 
			w4.write(str(data111[i]))
			i = i + 1
		w4.close()
		
		# execute LAMMPS for 100 and 111
	
		os.system("lmp_openmpi -log log.100 < newinput.100")
		os.system("lmp_openmpi -log log.111 < newinput.111")		

		# read output files and create force and energy matrix for 100 and 111

		rf100 = file("dump.100","r")
		dump100 = rf100.readlines()
		rf100.close()
		
		rf111 = file("dump.111","r")
		dump111 = rf111.readlines()
		rf111.close()

		i = 70
		while i < 77 :
			force = dump100[i].split()
			fx = float(force[1])
			fy = float(force[2])
			fz = float(force[3])
			F100FF.append([fx, fy, fz])
			force = dump111[i].split()
			fx = float(force[1])
			fy = float(force[2])
			fz = float(force[3])
			F111FF.append([fx, fy, fz])
			i = i + 1
	
		i = 86
		while i < 103 :
			force = dump100[i].split()
			fx = float(force[1])
			fy = float(force[2])
			fz = float(force[3])
			F100FF.append([fx, fy, fz])
			force = dump111[i].split()
			fx = float(force[1])
			fy = float(force[2])
			fz = float(force[3])
			F111FF.append([fx, fy, fz])
			i = i + 1		
		
		re100 = file("log.100","r")
		log100 = re100.readlines()
		re100.close()

		re111 = file("log.111","r")
		log111 = re111.readlines()
		re111.close()
		
		i = 71
		energy = log100[i].split()
		Ei100 = (float(energy[1]) - pe100) / N100
		E100FF.append(Ei100)
		energy = log111[i].split()
		Ei111 = (float(energy[1]) - pe111) / N111
		E111FF.append(Ei111)

		c = c + 1

	# calculation of cost function based on Eq. 8 in J. Phys. Chem. C 2014, 118, 3366-3374 

	# force terms
	
	FcostF100N = 0.0  # N = numerator		
	FcostF111N = 0.0
	
	i = 72
	while i < 360 :
		j = 0
		while j < 3 :
			FcostF100N = FcostF100N + np.power(F100FF[i][j] - F100DFT[i][j], 2)	
			FcostF111N = FcostF111N + np.power(F111FF[i][j] - F111DFT[i][j], 2)
			j = j + 1
		i = i + 1

	# energy terms

	FcostE100N = 0.0
	FcostE111N = 0.0
	FcostP100111N = 0.0

	i = 0
	while i < 15 :
		FcostE100N = FcostE100N + np.power(E100FF[i] - E100DFT[i], 2)
		FcostE111N = FcostE111N + np.power(E111FF[i] - E111DFT[i], 2)
		FcostP100111N = FcostP100111N + np.power(E100FF[i] - E111FF[i] - P100111DFT[i], 2)
		i = i + 1

	# cost function for the nth loop
	
	FcostF = FcostF100N / FcostF100D + FcostF111N / FcostF111D
	FcostE = w * (FcostE100N / FcostE100D + FcostE111N / FcostE111D)
	FcostP = w * FcostP100111N / FcostP100111D
	
	Fcostotal = FcostF + FcostE + FcostP
	Fcost.append(Fcostotal)
	
	# keep current parameters 

	if Fcost[n] < Fkeep[m] :
		m = m + 1
		
		Nkeep.append(n)
		Fkeep.append(Fcost[n])

		s6keep = s6

		fNkeep = fN
		fHkeep = fH
		
		D0HCkeep = D0HC
		alphaHCkeep = alphaHC
		r0HCkeep = r0HC
		
		D0HA2keep = D0HA2
		alphaHA2keep = alphaHA2
		r0HA2keep = r0HA2
		
		D0CT2keep = D0CT2
		alphaCT2keep = alphaCT2
		r0CT2keep = r0CT2
		
		D0NH2keep = D0NH2
		alphaNH2keep = alphaNH2
		r0NH2keep = r0NH2

	# choose new parameters
	
	s6 = newparameter(s6keep, range * 0)
	
	fN = newparameter(fNkeep, range * 0)
	fH = newparameter(fHkeep, range * 0)
	
	D0HC = newparameter(D0HCkeep, range)
	alphaHC = newparameter(alphaHCkeep, range) 
	r0HC = newparameter(r0HCkeep, range)
	
	D0HA2 = newparameter(D0HA2keep, range * 0)
	alphaHA2 = newparameter(alphaHA2keep, range * 0)
	r0HA2 = newparameter(r0HA2keep, range * 0)
	
	D0CT2 = newparameter(D0CT2keep, range * 0)
	alphaCT2 = newparameter(alphaCT2keep, range * 0)
	r0CT2 = newparameter(r0CT2keep, range * 0)
	
	D0NH2 = newparameter(D0NH2keep, range * 0)
	alphaNH2 = newparameter(alphaNH2keep, range * 0)
	r0NH2 = newparameter(r0NH2keep, range * 0)
	
	# save and overwrite restart file every 100 loops
 
	if n % 100 == 0 :
		w5 = file("SA.restart","w")
		
		w5.write("n = %.4s \n\n" % (n))
		
		w5.write("s6 = %.4s \n\n" % (s6keep))		

		w5.write("fN = %.4s \n" % (fNkeep))
		w5.write("fH = %.4s \n\n" % (fHkeep))
		
		w5.write("D0HC = %.7s \n" % (D0HCkeep))
		w5.write("alphaHC = %.4s \n" % (alphaHCkeep))
		w5.write("r0HC = %.4s \n\n" % (r0HCkeep))
		
		w5.write("D0HA2 = %.7s \n" % (D0HA2keep))
		w5.write("alphaHA2 = %.4s \n" % (alphaHA2keep))
		w5.write("r0HA2 = %.4s \n\n" % (r0HA2keep))
		
		w5.write("D0CT2 = %.7s \n" % (D0CT2keep))
		w5.write("alphaCT2 = %.4s \n" % (alphaCT2keep))
		w5.write("r0CT2 = %.4s \n\n" % (r0CT2keep))
		
		w5.write("D0NH2 = %.7s \n" % (D0NH2keep))
		w5.write("alphaNH2 = %.4s \n" % (alphaNH2keep))
		w5.write("r0NH2 = %.4s \n\n" % (r0NH2keep))
		
		w5.write("Fkeep = %.15s" % (Fkeep[m]))
		
		w5.close()		

	n = n + 1

# ---ending of simulated annealing loop---	
	
# output SA temperature, Morse potential parameters, force and energy profiles, and Fcost values

w = file("SA.out","w")

w.write("n = %.4s \n\n" % (n - 1))

w.write("s6 = %.4s \n\n" % (s6keep))

w.write("fN = %.4s \n" % (fNkeep))
w.write("fH = %.4s \n\n" % (fHkeep))

w.write("D0HC = %.7s \n" % (D0HCkeep))
w.write("alphaHC = %.4s \n" % (alphaHCkeep))
w.write("r0HC = %.4s \n\n" % (r0HCkeep))

w.write("D0HA2 = %.7s \n" % (D0HA2keep))
w.write("alphaHA2 = %.4s \n" % (alphaHA2keep))
w.write("r0HA2 = %.4s \n\n" % (r0HA2keep))

w.write("D0CT2 = %.7s \n" % (D0CT2keep))
w.write("alphaCT2 = %.4s \n" % (alphaCT2keep))
w.write("r0CT2 = %.4s \n\n" % (r0CT2keep))

w.write("D0NH2 = %.7s \n" % (D0NH2keep))
w.write("alphaNH2 = %.4s \n" % (alphaNH2keep))
w.write("r0NH2 = %.4s \n\n" % (r0NH2keep))

w.write("F100FF = \n")
i = 72
while i < 360 :
	w.write("%.15s  %.15s  %.15s \n" % (F100FF[i][0], F100FF[i][1], F100FF[i][2]))
	i = i + 1

w.write("\nF111FF = \n")
i = 72
while i < 360 :
	w.write("%.15s  %.15s  %.15s \n" % (F111FF[i][0], F111FF[i][1], F111FF[i][2]))
	i = i + 1
	
w.write("\nE100FF = \n")
i = 0
while i < 15 : 
	w.write("%.15s \n" % (E100FF[i]))
	i = i + 1
	
w.write("\nE111FF = \n")
i = 0
while i < 15 : 
	w.write("%.15s \n" % (E111FF[i]))
	i = i + 1
	
w.write("\nFcost = \n")
i = 0
while i <= (n - 1) : 
	w.write("%.15s \n" % (Fcost[i]))
	i = i + 1
	
w.write("\nNkeep = \n")
i = 0
while i <= (m - 1) : 
	w.write("%.7s \n" % (Nkeep[i]))
	i = i + 1

w.write("\nFkeep = \n")
i = 1
while i <= m : 
	w.write("%.15s \n" % (Fkeep[i]))
	i = i + 1
	
w.close()

print "\n Job finished \n"