# ----------Import necessary packages----------
import numpy as np
import matplotlib, matplotlib.pyplot as plt
from scipy.integrate import cumtrapz


# ----------Define the function we want to integrate----------

# Calculate 1/(adot/H0) where adot=da/dt for the chosen
# scalefactor (a), matter density (om), and cosmological constant density (ol).
# Note this is dimensionless, (it's Friedmann's equation divided by H0).
# adot should be multiplied by H0 to get units.


def getadotinv(a,om,ol):
    ok=1.0-om-ol                         # curvature (ok=Omega_k)
    adot = np.sqrt(ok + om/a + ol*a*a)
    return 1./adot

# Create a vectorized version of it that allows an array of a values
vecGetadotinv = np.vectorize(getadotinv, excluded=['om', 'ol'])

#-------------------
# INITIALISE ARRAYS
#-------------------

# Initialise an array of scalefactors
a_lo    = 1.e-10
a_hi    = 2
a_nstep = 1001  # by making an odd number of elements, you hit the centre element exactly
aarr    = np.linspace(a_lo, a_hi, a_nstep)  # Makes an array from a_lo to a_hi with the number of steps being a_nsteps

# Find which array index corresponds to a=1, i.e. corresponds to today
iToday = np.argmin(np.abs(aarr - 1.0))

# Initialise values of Omega_M to loop over
om_lo    = 0.0
om_hi    = 5
om_nstep = 6
oms      = np.linspace(om_lo, om_hi, om_nstep)

# Set a value of Omega_Lambda
ol=0.7


#--------------------------------------------
# CALCULATE Conversion factors to get into Gyr
#--------------------------------------------

H0kmsmpc = 70                 # Choose a value of H0 in km/s/Mpc
#H0s = H0kmsmpc * 3.24e-20     # Convert from km/s/Mpc to inverse seconds (3.24e-20 Mpc/km)
#H0y = H0s* 3.156e16           # Convert from inverse seconds to inverse Giga-years (3.15e16 s/Gyr)
H0y = H0kmsmpc * 0.001022544  # Convert directly to inverse Giga-years


#--------------------------------------------
# CALCULATE time as a function of scalefactor
#--------------------------------------------
# To calculate time we need to integrate the function 1/(da/dt), with respect to a
lookback = []
for om in oms:
    adotinvs = vecGetadotinv(aarr, om, ol)       # Get a vector of adotinvs as a function of a
    ages = cumtrapz(adotinvs, x=aarr, initial=0) # To a trapezoidal integral over that array.
    lookback.append((ages - ages[iToday])/H0y)   # Shift the arrays by subtracting the current age of the universe
                                                 # and convert to billions of years by dividing by H0y.

# Example analytic solutions for testing
#tom0_analytic = 1.0/H0y/np.sqrt(ol)*np.log(aarr)  # om=0, ol=1.0
#tom0_analytic2 = 1.0/H0y/np.sqrt(ol) * ( np.log(aarr+np.sqrt(aarr*aarr+(1.0-ol)/ol)) - np.log(1.0+np.sqrt(1.0+(1.0-ol)/ol)) ) # om=0, ol=ol, ok=1-ol

#--------------------------------------------
# PLOT scalefactor as a function of time
#--------------------------------------------
colours = matplotlib.cm.rainbow(np.linspace(0,1,oms.size))
fig = plt.figure()
ax = fig.add_subplot(1,1,1)
for (om, lookback, colour) in zip(oms, lookback, colours):
    ax.plot(lookback, aarr, label=r"$\Omega_m=%0.2f$"%om, color=colour)
#ax.plot(tom0_analytic,aarr,'--',color='black') # Example analytic solutions
#ax.plot(tom0_analytic2,aarr,'--',color='black')

# plt.xlim([-20,20])  # How to change axis limits if you need to

ax.axvline(0, color="k", ls=":")
ax.axhline(1, color="k", ls=":")
ax.legend(loc="upper left",frameon=False)#,bbox_to_anchor=(0, 0.5))
plt.figtext(0.88,0.1,"$\Omega_\Lambda=%0.2f$"%ol,ha='right',va='bottom',weight='roman', size='large')
ax.set_xlabel("$t - t_0$ (Gyr)", fontsize=16)
ax.set_ylabel("$a$", fontsize=16)
plt.show()
#fig.savefig("test.pdf", bbox_inches="tight")