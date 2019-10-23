N=6  #NUmber of layers of molecules

import math
import time
import numpy as np
import matplotlib.pyplot as plt
#import matplotlib as plt
import scipy
from scipy.optimize import basinhopping
angle = np.random.uniform(0,2*math.pi,(N,1))

#angle = np.arange(2, 3, 0.1)
len(angle)
r=np.arange(N)




##define POLAR .987
#define J1 2
#define J2 .897
#define J3 .66
#define T 1400

start_time = time.time()

J1=3.0E-4
J2=1.0E-5
J3=2.0E-5
TPU=187
TAF = 175
T=181
#P=(TPU-T)/B
P=0.987
field=100

def hamiltonian(angle):
    cost = 0
    i=0
    z=1
    last_row=len(angle)
    while i <(len(angle)-1):
        cost = cost + J1*P**2 *np.cos(angle[i+1]-angle[i] ) - J2*P**4* np.cos(angle[i+1]-angle[i])  + J3*np.cos(angle[i+1]-angle[i-1]) - P*field*np.cos(angle[i])
        i=i+1
    energy=cost
    return (energy)
sol = scipy.optimize.basinhopping(hamiltonian, angle, niter=100, T=1.0, stepsize=0.5, minimizer_kwargs=None, take_step=None, accept_test=None, callback=None, interval=50, disp=False, niter_success=None, seed=None)
final_config = sol.x
final_config=math.pi-final_config%(2*math.pi)
print("Total Runtime is--- %s seconds ---" % (time.time() - start_time))
ay=plt.plot(r, final_config, '-ro')
#plt.show(block=True)

#ax = plt.subplot(111, projection='polar')
#ax.plot(final_config*180/math.pi,r)
plt.show(block=True)
