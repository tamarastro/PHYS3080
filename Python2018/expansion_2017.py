# ----------Import necessary packages----------
import numpy as np
import matplotlib
matplotlib.use('TkAgg')
import matplotlib, matplotlib.pyplot as plt
from scipy.integrate import *


# ----------Define the function we want to integrate----------

# Calculate 1/(adot/H0) where adot=da/dt for the chosen 
# scalefactor, matter density, and cosmological constant denstiy
def getadotinv(a,om,ol):
	ok=1.0-om-ol
	adot = np.sqrt(ok + om/a + ol*a*a)
	return 1./adot
 
# Create a vectorized version of it that allows an array of a values
vecGetadotinv = np.vectorize(getadotinv, excluded=['om', 'ol'])

#-------------------
# INITIALISE ARRAYS
#-------------------

# Initialise an array of scalefactors
a_lo   = 0.01
a_hi   = 2
a_nstep = 101  # by making an odd number of elements, you hit the centre element exactly
aarr  = np.linspace(a_lo, a_hi, a_nstep)  # Makes an array from a_lo to a_hi with the number of steps being a_nsteps

# Find which array index corresponds to a=1, i.e. corresponds to today 
iToday = np.argmin(np.abs(aarr - 1.0))

# Initialise values of OM to loop over
om_lo   = 0;
om_hi   = 2;
om_nstep = 6;
oms   = np.linspace(om_lo, om_hi, om_nstep)

# Set a value of OL
ol=0.7


#--------------------------------------------
# CALCULATE Conversion factors to get into Gyr
#--------------------------------------------

H0kmsmpc = 70                 # Choose a value of H0 in km/s/Mpc
H0s = H0kmsmpc * 3.24e-20     # Convert to inverse seconds
H0y = H0s* 3.156e16           # Convert to inverse Giga-years 


#--------------------------------------------
# CALCULATE time as a function of scaelfactor
#--------------------------------------------
# To calculate time we need to integrate the function 1/(da/dt), with respect to a 
lookback = []
for om in oms:
    adotinvs = vecGetadotinv(aarr, om, ol)
    ages = cumtrapz(adotinvs, x=aarr, initial=0)
    lookback.append((ages - ages[iToday])/H0y)
    



colours = matplotlib.cm.rainbow(np.linspace(0,1,oms.size))
fig = plt.figure()
ax = fig.add_subplot(1,1,1)
for (om, lookback, colour) in zip(oms, lookback, colours):
    ax.plot(lookback, aarr, label=r"$\Omega_m=%0.2f$"%om, color=colour)
ax.axvline(0, color="k", ls=":")
ax.axhline(1, color="k", ls=":")
ax.legend(loc="upper left",frameon=False)#,bbox_to_anchor=(0, 0.5))
plt.figtext(0.88,0.1,"$\Omega_\Lambda=%0.2f$"%ol,ha='right',va='bottom',weight='roman', size='large')
ax.set_xlabel("$t - t_0$ (Gyr)", fontsize=16)
ax.set_ylabel("$a$", fontsize=16)
plt.show()
#fig.savefig("test.pdf", bbox_inches="tight")