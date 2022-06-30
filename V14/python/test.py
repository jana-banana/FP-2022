import matplotlib.pyplot as plt
import numpy as np
import scipy.constants as const
from scipy.optimize import curve_fit
from uncertainties import ufloat
import uncertainties.unumpy as unp
from numpy import linalg 

import os

if os.path.exists("../build") == False:
    os.mkdir("../build")

#Für MAKE das ../ entfernen

a = 1 #cm

#Durchgäng - Dicke
d_i1 = 2*np.sqrt(2)*a #cm
d_i2 = 3*np.sqrt(2)*a #cm
d_i3 = 2*np.sqrt(2)*a #cm
d_i4 = 3*a #cm
d_i5 = 3*a #cm
d_i6 = 3*a #cm
d_i7 = 2*np.sqrt(2)*a #cm
d_i8 = 3*np.sqrt(2)*a #cm
d_i9 = 2*np.sqrt(2)*a #cm
d_i10 = 3*a #cm
d_i11 = 3*a #cm
d_i12 = 3*a #cm

#Time
t_0 = 300 #sec

# ---------------------------------Nullmessung/Spektrum----------------------------------------------------------------------------------
pulses = np.genfromtxt('../data/spektrum.txt', unpack=True)
x_axis = np.arange(len(pulses))

r_0 = np.sum(pulses)/t_0 #/per/sec
sig_0 = np.sqrt(np.sum(pulses))/t_0

R_0 = unp.uarray([r_0], [sig_0])

plt.figure()
plt.bar(x_axis, pulses)
plt.xlim(8,220)
plt.xlabel("Channel")
plt.ylabel("Anzahl der Ereignisse")
plt.grid()
plt.tight_layout()
plt.savefig('../build/spektrum.pdf')

#--------------------------------Würfel 1/ Aluminiumgehäuse-------------------------------------------------------------------------------

w1_hauptdiag = np.genfromtxt('../data/w1_hauptdiag.txt', unpack=True)

# N_1h = np.sum(w1_hauptdiag)
N_1h = np.amax(w1_hauptdiag, axis=0)
sigN_1h = np.sqrt(N_1h)

r_1h = N_1h/t_0
sig_1h = np.sqrt(N_1h)/t_0


w1_nebendiag = np.genfromtxt('../data/w1_nebendiag.txt', unpack=True)

# N_1n = np.sum(w1_nebendiag)
N_1n = np.amax(w1_nebendiag, axis=0)
sigN_1n = np.sqrt(N_1n)

r_1n = N_1n/t_0
sig_1n = np.sqrt(N_1n)/t_0


w1_durch = np.genfromtxt('../data/w1_durch.txt', unpack=True)

# N_1d = np.sum(w1_durch)
N_1d = np.amax(w1_durch, axis=0)
sigN_1d = np.sqrt(N_1d)

r_1d = N_1d/t_0
sig_1d = np.sqrt(N_1d)/t_0


print(f'Würfel 1 - Aluminiumhülle zählraten \n Hauptdiagonale :{r_1h:.6f} \\pm {sig_1h:.6f} \n Nebendiagonale: {r_1n:.6f} \\pm {sig_1n:.6f} \n Gerade durch:{r_1d:.6f} \\pm {sig_1d:.6f} \n')
print(f'Würfel 1 - Aluminiumhülle counts \n Hauptdiagonale :{N_1h:.1f} \\pm {sigN_1h:.1f} \n Nebendiagonale: {N_1n:.1f} \\pm {sigN_1n:.1f} \n Gerade durch:{N_1d:.1f} \\pm {sigN_1d:.1f} \n')


R_1h = unp.uarray([r_1h],[sig_1h])
R_1n = unp.uarray([r_1n],[sig_1n])
R_1d = unp.uarray([r_1d],[sig_1d])

#--------------------------------Würfel 2--------------------------------------------------------------------------------------------------
print(f'_____________________Würfel 2 - Gleiches Mu für alle_______________________')

names = ['1','2','7','8','11','12']
numbers = np.arange(len(names))
dicke_w2 = [d_i1, d_i2, d_i7, d_i8, d_i11, d_i12]

I_w2 = np.zeros(len(names))
sigI_w2 = np.zeros(len(names))
N_w2 = np.zeros(len(names))
sigN_w2 = np.zeros(len(names))
mu_w2 = np.zeros(len(names))
sig_w2 = np.zeros(len(names))

for x, n, d in zip(names, numbers, dicke_w2):

    w2_ix = np.genfromtxt('../data/w2_i'+x+'.txt', unpack=True)

    # N_ix = np.sum(w2_ix)
    N_ix = np.amax(w2_ix, axis=0)
    sigN_ix = np.sqrt(N_ix)

    r_w2_ix = N_ix/t_0
    sig_w2_ix = sigN_ix/t_0

    R_w2_ix = unp.uarray([r_w2_ix],[sig_w2_ix])

    if d == 3*np.sqrt(2)*a:
        mu_w2_ix = unp.log(R_1h/R_w2_ix)/ d
    
    elif d == 2*np.sqrt(2)*a:
        mu_w2_ix = unp.log(R_1n/R_w2_ix)/ d

    elif d == 3*a:
        mu_w2_ix = unp.log(R_1d/R_w2_ix)/ d

    N_w2[n] = N_ix
    sigN_w2[n] = sigN_ix 
    
    I_w2[n] = unp.nominal_values(R_w2_ix)
    sigI_w2[n] = unp.std_devs(R_w2_ix)
    
    mu_w2[n] = unp.nominal_values(mu_w2_ix)
    sig_w2[n] = unp.std_devs(mu_w2_ix)

print("Counts")
for n in range(6):
    print(f'{n+1}: {N_w2[n]:.1f} \\pm {sigN_w2[n]:.1f}')  

print("Intensitäten")
for n in range(6):
    print(f'{n+1}: {I_w2[n]:.6f} \\pm {sigI_w2[n]:.6f}')

print("Mu's")
for n in range(6):
    print(f'\n I{n+1}: mu = {mu_w2[n]:.6f} \\pm {sig_w2[n]:.6f} 1/cm')
    

mu_w2_mean = np.mean(unp.nominal_values(mu_w2))
sig_w2_mean = np.std(unp.nominal_values(mu_w2))/len(names)
mu_w2_m = unp.uarray(mu_w2_mean, sig_w2_mean)

abw = 1 - mu_w2_m/0.096

print(f'\n Würfel 2 \n mu mean = {mu_w2_mean:.6f} \\pm {sig_w2_mean:.6f} \n')
print(f'Abweichung zu Delrin: {abw}')



#--------------------------------Würfel 3--------------------------------------------------------------------------------------------------
print(f'_____________________Würfel 3 - Gleiches Mu für alle_______________________')
t_w3 = 2*300 #sec

names = ['5','7','8']
numbers = np.arange(len(names))
dicke_w3 = [d_i5, d_i7, d_i8]
I_w3 = np.zeros(len(names))
sigI_w3 = np.zeros(len(names))
N_w3 = np.zeros(len(names))
sigN_w3 = np.zeros(len(names))
mu_w3 = np.zeros(len(names))
sig_w3 = np.zeros(len(names))

for x, n, d in zip(names, numbers, dicke_w3):

    w3_ix_1 = np.genfromtxt('../data/w3_i'+x+'_1.txt', unpack=True)
    w3_ix_2 = np.genfromtxt('../data/w3_i'+x+'_2.txt', unpack=True)

    # N_ix = np.sum(w3_ix_1) + np.sum(w3_ix_2)
    N_ix = np.amax(w3_ix_1, axis=0) + np.amax(w3_ix_2, axis=0)
    sigN_ix = np.sqrt(N_ix)

    r_w3_ix = N_ix/t_w3
    sig_w3_ix = sigN_ix/t_w3

    R_w3_ix = unp.uarray([r_w3_ix],[sig_w3_ix])

    if d == 3*np.sqrt(2)*a:
        mu_w3_ix = unp.log(2*R_1h/R_w3_ix)/ d
    
    elif d == 2*np.sqrt(2)*a:
        mu_w3_ix = unp.log(2*R_1n/R_w3_ix)/ d

    elif d == 3*a:
        mu_w3_ix = unp.log(2*R_1d/R_w3_ix)/ d

    N_w3[n] = N_ix
    sigN_w3[n] = sigN_ix 
    
    I_w3[n] = unp.nominal_values(R_w3_ix)
    sigI_w3[n] = unp.std_devs(R_w3_ix)
    
    mu_w3[n] = unp.nominal_values(mu_w3_ix)
    sig_w3[n] = unp.std_devs(mu_w3_ix)

print("Counts")
for n in range(3):
    print(f'{n+1}: {N_w3[n]:.1f} \\pm {sigN_w3[n]:.1f}')  

print("Intensitäten")
for n in range(3):
    print(f'{n+1}: {I_w3[n]:.2f} \\pm {sigI_w3[n]:.2f}')

print("Mu's")
for n in range(3):
    print(f'\n I{n+1}: mu = {mu_w3[n]:.6f} \\pm {sig_w3[n]:.6f} 1/cm')

mu_w3_mean = np.mean(unp.nominal_values(mu_w3))
sig_w3_mean = np.std(unp.nominal_values(mu_w3))/len(names)
mu_w3_m = unp.uarray(mu_w3_mean, sig_w3_mean)

abw = 1 - mu_w3_m/1.415

print(f'\n Würfel 3 \n mu mean = {mu_w3_mean:.6f} \\pm {sig_w3_mean:.6f} \n')
print(f'Abweichung zu Blei: {abw}')


#--------------------------------Würfel 4--------------------------------------------------------------------------------------------------
print(f'_____________________Würfel 4 - unbekannte Bestimmug_______________________')

##Matrix A
w = np.sqrt(2)
a1 = np.array([(0,w,0,w,0,0,0,0,0),(0,0,w,0,w,0,w,0,0),(0,0,0,0,0,w,0,w,0)])
a2 = np.array([(1,1,1,0,0,0,0,0,0),(0,0,0,1,1,1,0,0,0),(0,0,0,0,0,0,1,1,1)])
a3 = np.array([(0,w,0,0,0,w,0,0,0),(w,0,0,0,w,0,0,0,w),(0,0,0,w,0,0,0,w,0)])
a4 = np.array([(0,0,1,0,0,1,0,0,1),(0,1,0,0,1,0,0,1,0),(1,0,0,1,0,0,1,0,0)])
A = a*np.vstack((a1,a2,a3,a4))

##Vektor I
names = ['1','2','3','4','5','6','7','8','9','10','11','12']
numbers = np.arange(len(names))
I_w4 = np.zeros(len(names))
sig_w4 = np.zeros(len(names))
N_w4 = np.zeros(len(names))
sigN_w4 = np.zeros(len(names))

for x, n in zip(names, numbers):

    w4_ix = np.genfromtxt('../data/w4_i'+x+'.txt', unpack=True)

    # N_ix = np.sum(w4_ix)
    N_ix = np.amax(w4_ix, axis=0)
    sigN_ix = np.sqrt(N_ix)

    r_w4_ix = N_ix/t_0
    sig_w4_ix = sigN_ix/t_0

    R_w4_ix = unp.uarray([r_w4_ix],[sig_w4_ix])

    if x in ['2','8']:
        I_w4_ix = unp.log(R_1h/R_w4_ix)
    
    elif x in ['1','3','7','9']:
        I_w4_ix = unp.log(R_1n/R_w4_ix)

    elif x in ['4','5','6','10','11','12']:
        I_w4_ix = unp.log(R_1d/R_w4_ix)

    N_w4[n] = N_ix
    sigN_w4[n] = sigN_ix 
    
    I_w4[n] = unp.nominal_values(I_w4_ix)
    sig_w4[n] = unp.std_devs(I_w4_ix)

print("Counts")
for n in range(12):
    print(f'{n+1}: {N_w4[n]:.1f} \\pm {sigN_w4[n]:.1f}')  

print("Intensitäten")
for n in range(12):
    print(f'{n+1}: {I_w4[n]:.2f} \\pm {sig_w4[n]:.2f}')


mu_w4 = np.linalg.inv(A.T @ A) @ A.T @ I_w4
    
##Kovarianzmatrix 
V_I = np.diag(sig_w4**2, k=0)
B = np.linalg.inv(A.T @ A) @ A.T
V_mu = B @ V_I @ B.T

sig_mu_w4 = np.sqrt(np.diag(V_mu))

print("mu mit fehler")
for n in range(9):
    print(f'{n+1}: {mu_w4[n]:.6f} \\pm {sig_mu_w4[n]:.6f}')

## Abweichung von Literaturwert
mu_4_a = unp.uarray(mu_w4, sig_mu_w4)
lit = np.array([0.211, 1.415, 0.211, 0.060, 1.415, 0.060, 0.096, 1.415, 0.060])
abw = 1 - mu_4_a/lit

for n in range(9):
    print(f'{n+1}: {abw[n]*100}')