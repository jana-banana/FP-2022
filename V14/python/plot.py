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

# ---------------------------------Nullmessung/Spektrum----------------------------------------------------------------------------------
pulses = np.genfromtxt('../data/spektrum.txt', unpack=True)
t_0 = 300 #sec
x_axis = np.arange(len(pulses))

r_0 = np.sum(pulses)/t_0 #/per/sec
sig_0 = np.sqrt(np.sum(pulses))/t_0

R_0 = unp.uarray([r_0], [sig_0])

plt.figure()
plt.bar(x_axis, pulses)
plt.xlabel("Channel")
plt.ylabel("Anzahl der Ereignisse")
plt.tight_layout()
plt.savefig('../build/spektrum.pdf')

#--------------------------------Würfel 1/ Aluminiumgehäuse-------------------------------------------------------------------------------

w1_hauptdiag = np.genfromtxt('../data/w1_hauptdiag_f.txt', unpack=True)
t_w1_hauptdiag = 120 #sec

r_1h = np.sum(w1_hauptdiag)/t_w1_hauptdiag
sig_1h = np.sqrt(np.sum(w1_hauptdiag))/t_w1_hauptdiag


# w1_nebendiag = np.genfromtxt('../data/w1_nebendiag_f.txt', unpack=True)
t_w1_nebendiag = 479.40 #sec

# r_1n = np.sum(w1_nebendiag)/t_w1_nebendiag
r_1n = 102568/t_w1_nebendiag
# sig_1n = np.sqrt(np.sum(w1_nebendiag))/t_w1_nebendiag
sig_1n = np.sqrt(102568)/t_w1_nebendiag


w1_durch = np.genfromtxt('../data/w1_durch_f.txt', unpack=True)
t_w1_durch = 120 #sec

r_1d = np.sum(w1_durch)/t_w1_durch
sig_1d = np.sqrt(np.sum(w1_durch))/t_w1_durch


print(f'Würfel 1 - Aluminiumhülle \n Hauptdiagonale :{r_1h} \\pm {sig_1h} \n Nebendiagonale: {r_1n} \\pm {sig_1n} \n Gerade durch:{r_1d} \\pm {sig_1d} \n')

R_1h = unp.uarray([r_1h],[sig_1h])
R_1n = unp.uarray([r_1n],[sig_1n])
R_1d = unp.uarray([r_1d],[sig_1d])

#--------------------------------Würfel 2--------------------------------------------------------------------------------------------------
print(f'_____________________Würfel 2 - Gleiches Mu für alle_______________________')
t_w2_i1 = 180 #sec

t_w2_i2 = 180 #sec

t_w2_i3 = 180 #sec

t_w2_i4 = 180 #sec

t_w2_i5 = 180 #sec

t_w2_i6 = 180 #sec

t_w2_i7 = 180 #sec

t_w2_i8 = 180 #sec

t_w2_i9 = 180 #sec

t_w2_i10 = 180 #sec

t_w2_i11 = 180 #sec

t_w2_i12 = 180 #sec


names = ['1','2','3','4','5','6','7','8','9','10','11','12']
numbers = np.arange(len(names))
times_w2 = [t_w2_i1, t_w2_i2, t_w2_i3, t_w2_i4, t_w2_i5, t_w2_i6, t_w2_i7, t_w2_i8, t_w2_i9, t_w2_i10, t_w2_i11, t_w2_i12]
dicke_w2 = [d_i1, d_i2, d_i3, d_i4, d_i5, d_i6, d_i7, d_i8, d_i9, d_i10, d_i11, d_i12]
mu_w2 = np.zeros(len(names))
sig_w2 = np.zeros(len(names))

for x, n, t, d in zip(names, numbers, times_w2, dicke_w2):

    w2_ix = np.genfromtxt('../data/w2_i'+x+'.txt', unpack=True)

    r_w2_ix = np.sum(w2_ix)/t
    sig_w2_ix = np.sqrt(np.sum(w2_ix))/t

    R_w2_ix = unp.uarray([r_w2_ix],[sig_w2_ix])

    if d == 3*np.sqrt(2)*a:
        mu_w2_ix = unp.log(R_1h/R_w2_ix)/ d
    
    elif d == 2*np.sqrt(2)*a:
        mu_w2_ix = unp.log(R_1n/R_w2_ix)/ d

    elif d == 3*a:
        mu_w2_ix = unp.log(R_1d/R_w2_ix)/ d
    
    mu_w2[n] = unp.nominal_values(mu_w2_ix)
    sig_w2[n] = unp.std_devs(mu_w2_ix)
    
    print(f'\n I{x}: mu = {mu_w2_ix} 1/cm')

# mu_w2_a = unp.uarray(mu_w2, sig_w2)
mu_w2_mean = unp.uarray([np.sum(mu_w2)/len(names)], [np.sum(sig_w2)/len(names)])

print(f'\n Würfel 2 \n mu mean = {mu_w2_mean} \n')

#--------------------------------Würfel 3--------------------------------------------------------------------------------------------------
print(f'_____________________Würfel 3 - Gleiches Mu für alle_______________________')
t_w3_i1 = 180 #sec

t_w3_i2 = 180 #sec

t_w3_i3 = 180 #sec

t_w3_i4 = 180 #sec

t_w3_i5 = 180 #sec

t_w3_i6 = 180 #sec

t_w3_i7 = 180 #sec

t_w3_i8 = 180 #sec

t_w3_i9 = 180 #sec

t_w3_i10 = 180 #sec

t_w3_i11 = 180 #sec

t_w3_i12 = 180 #sec


names = ['1','2','3','4','5','6','7','8','9','10','11','12']
numbers = np.arange(len(names))
times_w3 = [t_w3_i1, t_w3_i2, t_w3_i3, t_w3_i4, t_w3_i5, t_w3_i6, t_w3_i7, t_w3_i8, t_w3_i9, t_w3_i10, t_w3_i11, t_w3_i12]
dicke_w3 = [d_i1, d_i2, d_i3, d_i4, d_i5, d_i6, d_i7, d_i8, d_i9, d_i10, d_i11, d_i12]
mu_w3 = np.zeros(len(names))
sig_w3 = np.zeros(len(names))

for x, n, t, d in zip(names, numbers, times_w3, dicke_w3):

    w3_ix = np.genfromtxt('../data/w3_i'+x+'.txt', unpack=True)

    r_w3_ix = np.sum(w3_ix)/t
    sig_w3_ix = np.sqrt(np.sum(w3_ix))/t

    R_w3_ix = unp.uarray([r_w3_ix],[sig_w3_ix])

    if d == 3*np.sqrt(2)*a:
        mu_w3_ix = unp.log(R_1h/R_w3_ix)/ d
    
    elif d == 2*np.sqrt(2)*a:
        mu_w3_ix = unp.log(R_1n/R_w3_ix)/ d

    elif d == 3*a:
        mu_w3_ix = unp.log(R_1d/R_w3_ix)/ d
    
    mu_w3[n] = unp.nominal_values(mu_w3_ix)
    sig_w3[n] = unp.std_devs(mu_w3_ix)
    
    print(f'\n I{x}: mu = {mu_w3_ix} 1/cm')

# mu_w3_a = unp.uarray(mu_w3, sig_w3)
mu_w3_mean = unp.uarray([np.sum(mu_w3)/len(names)], [np.sum(sig_w3)/len(names)])

print(f'\n Würfel 3 \n mu mean = {mu_w3_mean} \n')


#--------------------------------Würfel 4--------------------------------------------------------------------------------------------------
print(f'_____________________Würfel 4 - unbekannte Bestimmug_______________________')
t_w4_i1 = 180 #sec
t_w4_i2 = 180 #sec
t_w4_i3 = 180 #sec
t_w4_i4 = 180 #sec
t_w4_i5 = 180 #sec
t_w4_i6 = 180 #sec
t_w4_i7 = 180 #sec
t_w4_i8 = 180 #sec
t_w4_i9 = 180 #sec
t_w4_i10 = 180 #sec
t_w4_i11 = 180 #sec
t_w4_i12 = 180 #sec

# Matrix A
w = np.sqrt(2)
a1 = np.array([0,w,0,w,0,0,0,0,0],[0,0,w,0,w,0,w,0,0],[0,0,0,0,0,w,0,w,0])
a2 = np.array([1,1,1,0,0,0,0,0,0],[0,0,0,1,1,1,0,0,0],[0,0,0,0,0,0,1,1,1])
a3 = np.array([0,w,0,0,0,w,0,0,0],[w,0,0,0,w,0,0,0,w],[0,0,0,w,0,0,0,w,0])
a4 = np.array([0,0,1,0,0,1,0,0,1],[0,1,0,0,1,0,0,1,0],[1,0,0,1,0,0,1,0,0])
A = a*np.vstack(a1,a2,a3,a4)

#Vektor I
names = ['1','2','3','4','5','6','7','8','9','10','11','12']
numbers = np.arange(len(names))
times_w4 = [t_w4_i1, t_w4_i2, t_w4_i3, t_w4_i4, t_w4_i5, t_w4_i6, t_w4_i7, t_w4_i8, t_w4_i9, t_w4_i10, t_w4_i11, t_w4_i12]
mu_w4 = np.zeros(len(names))
sig_w4 = np.zeros(len(names))

for x, n, t in zip(names, numbers, times):

    w4_ix = np.genfromtxt('../data/w4_i'+x+'.txt', unpack=True)

    r_w4_ix = np.sum(w4_ix)/t
    sig_w4_ix = np.sqrt(np.sum(w4_ix))/t

    R_w4_ix = unp.uarray([r_w4_ix],[sig_w4_ix])

    if x in ['2','8']:
        I_w4_ix = unp.log(R_1h/R_w4_ix)
    
    elif x in ['1','3','7','9']:
        I_w4_ix = unp.log(R_1n/R_w4_ix)

    elif x in ['4','5','6','10','11','12']:
        I_w4_ix = unp.log(R_1d/R_w4_ix)
    
    I_w4[n] = unp.nominal_values(I_w4_ix)
    sig_w4[n] = unp.std_devs(I_w4_ix)

mu_w4 = np.linalg.inv(A.T @ A) @ A.T @ I_w4

for n in numbers:
    print(f'{mu_w4[n]}')
    
#Kovarianzmatrix 
V_I = np.diag(sig_w4**2, k=0)
B = np.linalg.inv(A.T @ A) @ A.T
V_mu = B @ V_I @ B.T

sig_mu_w4 = np.sqrt(np.diag(V_mu))

for n in numbers:
    print(f'{sig_mu_w4[n]}')