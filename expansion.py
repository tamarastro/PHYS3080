# ----------Import necessary packages----------
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from scipy.integrate import cumtrapz


def get_adot_inv(a, omega_matter, omega_lambda):
    """ Calculate 1/(adot/H0) where adot=da/dt
    
    Note this is dimensionless, (it's Friedmann's equation divided by H0).
    adot should be multiplied by H0 to get units.

    """
    omega_k = 1.0 - omega_matter - omega_lambda  # curvature (ok=Omega_k)
    adot = np.sqrt(omega_k + omega_matter / a + omega_lambda * a * a)
    return 1.0 / adot


# ---------- Initialise our variables ---------------------------
a_arr = np.linspace(1e-10, 2, 1001)  # Makes an array of scalefactors. 1001 steps from 1e-10 to 2
omega_matter_arr = np.linspace(0.0, 5.0, 6)  # Makes an array of matter densities. 6 steps from 0.0 to 5.0

# Set up other helpful variables
omega_lambda = 0.7
ind = np.argmin(np.abs(a_arr - 1.0))  # Find which array index corresponds to a=1 (i.e. today)

# Set up conversion factors
H0kmsmpc = 70  # Choose a value of H0 in km/s/Mpc
H0y = H0kmsmpc * 0.001022544  # Convert directly to inverse Giga-years


# ---------- Do the math ---------------------------
# Now for the fun part! For every matter density we want to calculate the lookback time, by integrating 1/adot wrt a
lookbacks = []
for om in omega_matter_arr:
    adotinvs = get_adot_inv(a_arr, om, omega_lambda)  # Get a vector of adotinvs as a function of a
    ages = cumtrapz(adotinvs, x=a_arr, initial=0)  # To a trapezoidal integral over that array.
    lookbacks.append((ages - ages[ind]) / H0y)  # Subtracting the current age of the universe and get into gigayears


# ---------- Plot everything ---------------------------

colours = matplotlib.cm.rainbow(np.linspace(0, 1, omega_matter_arr.size))  # Get us some colours

fig, ax = plt.subplots(nrows=1, ncols=1)
for (om, lookback, colour) in zip(omega_matter_arr, lookbacks, colours):
    ax.plot(lookback, a_arr, label=r"$\Omega_m=%0.2f$" % om, color=colour)

# If we want to plot some analytic solutions, uncomment the relevant lines
# t_om0_analytic = 1.0/H0y/np.sqrt(omega_lambda)*np.log(a_arr)  # om=0, ol=1.0
# ax.plot(t_om0_analytic, a_arr, '--', color='k', label=r"$\Omega_m=0, \ \Omega_\Lambda=1$")
# t_om0_analytic2 = 1.0/H0y/np.sqrt(omega_lambda) * (np.log(a_arr+np.sqrt(a_arr*a_arr+(1.0-omega_lambda)/omega_lambda)) - np.log(1.0+np.sqrt(1.0+(1.0-omega_lambda)/omega_lambda)))  # om=0, ol=ol, ok=1-ol
# ax.plot(t_om0_analytic2, a_arr, '--', color='gray', label=r"Flat $\Omega_m=0, \ \Omega_\Lambda=%0.2f$" % omega_lambda)

ax.set_xlim([-20, 20])  # How to change axis limits if you need to
ax.axvline(0, color="k", ls=":")
ax.axhline(1, color="k", ls=":")
ax.legend(loc="upper left", frameon=False)
plt.figtext(0.88, 0.12, r"$\Omega_\Lambda=%0.2f$" % omega_lambda, ha='right', va='bottom', weight='roman', size='large')
ax.set_xlabel("$t - t_0$ (Gyr)", fontsize=16)
ax.set_ylabel("$a$", fontsize=16)
fig.savefig("plots/expansion.png", bbox_inches="tight", transparent=True)
# plt.show()  # Show the plot interactively if you want
# For reports, save a pdf plot. For quick viewing, save a png
# fig.savefig("plots/expansion.pdf", bbox_inches="tight", transparent=True)
