import numpy as np
import math
import matplotlib.pyplot as plt

#Sample program, showing how to numerically make little time steps:
#$start by adding in all my necessary constants:
pi = math.pi
me = 9.1*10**(-31)
c = 2.998*10**(8)
qe = 1.602*10**(-19)
a0 = 5.29177* 10**(-11)

#declare an initial speed
v_0 =  [0.1*c, 0, 0]
#declare an initial position:
x_0 =  [0, 0, 0]
#declare a magnetic field, in Tesla:
B_0 = [0, 0, 50*10**(-4)]

#Calculate cyclotron radius
R_cyc = me*np.linalg.norm(v_0)/(np.linalg.norm(B_0)*qe)

#let's make a plot with this circle:
Circ = np.zeros([2,100]);
#increment:
inc = np.arange(100)

Circ[0,:] = R_cyc*np.cos(inc*2*pi/99)
Circ[1,:] = R_cyc*np.sin(inc*2*pi/99)-R_cyc

#First plot is our calculated 'ideal' circle
fig1 = plt.figure()
cp1 = plt.plot(Circ[0,:], Circ[1,:],'r')
plt.show()


#now let's use our time-incremented force equation
#remember F = qv x B
#initialize some arrays:
Motion_r = np.zeros([30000,3], dtype=np.longdouble);
Motion_v = np.zeros([300000,3], dtype=np.longdouble);
Motion_a = np.zeros([300000,3], dtype=np.longdouble);

#initial position
Motion_r[0,:] = x_0;
Motion_v[0,:] = v_0;

##Timestep is here:
delt = 1*10**(-10);
for i in range(0,29999):
    Bhere = [0,0, 50*10**(-4)]
    #Your job to calculated the current acceleration here:
    ai =  qe * np.cross(Motion_v[i,:], Bhere)/me
    #Then use this to change your speed and position:  vnew = vold + a*delt, xnew = xold + vold*delt + 1/2*a*delt^2
    Motion_v[i+1,:] = Motion_v[i,:] + ai * delt;
    Motion_r[i+1,:] = Motion_r[i,:] + Motion_v[i,:] * delt + 0.5* ai * delt**2;

"""
What went wrong is that our electron spirals into itself, when it should have a constant radius.
The timestep is too long compared to the speed of the particle, which lets the 
particle travel much further than it needs to before changing its direction
"""

#second plot is our calculated motion
fig2 = plt.figure()
cp2 = plt.plot(Motion_r[:,0], Motion_r[:,1],'b')
plt.show()

#better timestep
delt = 1*10**(-12);
for i in range(0,29999):
    Bhere = [0,0, 50*10**(-4)]
    ai =  qe * np.cross(Motion_v[i,:], Bhere)/me
    Motion_v[i+1,:] = Motion_v[i,:] + ai * delt;
    Motion_r[i+1,:] = Motion_r[i,:] + Motion_v[i,:] * delt + 0.5* ai * delt**2;


fig3 = plt.figure()
cp3 = plt.plot(Motion_r[:,0], Motion_r[:,1],'b')
plt.show()

#Problem 3
delt = 1*10**(-12);
for i in range(0,29999):
    Bhere = [0,0, 1000*10**(-4)*np.abs(Motion_r[i,1])]
    ai =  qe * np.cross(Motion_v[i,:], Bhere)/me
    Motion_v[i+1,:] = Motion_v[i,:] + ai * delt;
    Motion_r[i+1,:] = Motion_r[i,:] + Motion_v[i,:] * delt + 0.5* ai * delt**2;


fig4 = plt.figure()
cp4 = plt.plot(Motion_r[:,0], Motion_r[:,1],'b')
plt.show()



delt = 1*10**(-12);
for i in range(0,29999):
    Bhere = [0,0, 1000*10**(-4)*abs(Motion_r[i,0]**2 + Motion_r[i,1]**2)]
    ai =  qe * np.cross(Motion_v[i,:], Bhere)/me
    Motion_v[i+1,:] = Motion_v[i,:] + ai * delt;
    Motion_r[i+1,:] = Motion_r[i,:] + Motion_v[i,:] * delt + 0.5* ai * delt**2;


fig4 = plt.figure()
cp4 = plt.plot(Motion_r[:,0], Motion_r[:,1],'b')
plt.show()


delt = 1*10**(-12);
for i in range(0,29999):
    Bhere = [0,0, 1000*10**(-4)*25*np.linalg.norm(Motion_r[i,:])]
    ai =  qe * np.cross(Motion_v[i,:], Bhere)/me
    Motion_v[i+1,:] = Motion_v[i,:] + ai * delt;
    Motion_r[i+1,:] = Motion_r[i,:] + Motion_v[i,:] * delt + 0.5* ai * delt**2;


fig5 = plt.figure()
cp5 = plt.plot(Motion_r[:,0], Motion_r[:,1],'b')
plt.show()
