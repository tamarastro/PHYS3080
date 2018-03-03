import os
import numpy as np
from scipy.integrate import quad
from chainconsumer import ChainConsumer


def Hz_inverse(z, om, ol):
    """ Calculate 1/H(z). Will integrate this function. """
    Hz = np.sqrt((1 + z) ** 2 * (om * z + 1) - ol * z * (z + 2))
    return 1.0 / Hz


def dist_mod(zs, om, ol):
    """ Calculate the distance modulus, correcting for curvature"""
    ok = 1.0 - om - ol
    x = np.array([quad(Hz_inverse, 0, z, args=(om, ol))[0] for z in zs])
    if ok < 0.0:
        R0 = 1 / np.sqrt(-ok)
        D = R0 * np.sin(x / R0)
    elif ok > 0.0:
        R0 = 1 / np.sqrt(ok)
        D = R0 * np.sinh(x / R0)
    else:
        D = x
    lum_dist = D * (1 + zs)
    dist_mod = 5 * np.log10(lum_dist)
    return dist_mod


# ---------- Import data -------------------------------------
data = np.genfromtxt("data/data.txt", names=True)
zz = data["Z"]
mu = data["MU"]
mu_error = data["MUERR"] + 0.02
mu_error2 = mu_error ** 2
# Define cosntants
H0 = 70.0
c_H0 = 3.e5 / H0

# ---------- Set up fitting ranges ---------------------------
n = 10  # Increase this for a finer grid
oms = np.linspace(0, 0.7, n)
ols = np.linspace(0, 1.0, n)
n_marg = 200
mscr_guess = 5.0 * np.log10(c_H0) + 25
mscr = np.linspace(mscr_guess - 0.5, mscr_guess + 0.5, n_marg)

chi2 = np.ones((n, n)) * np.inf  # Our chi2 values, set initially to super large values
mscr_used = np.zeros((n, n))

# ---------- Do the fit ---------------------------
saved_data_filename = "data/saved_grid_%d.txt" % n

if os.path.exists(saved_data_filename):  # Load the last run with n grid if we can find it
    print("Loading saved data. Delete %s to regenerate the points\n" % saved_data_filename)
    chi2 = np.loadtxt(saved_data_filename)
else:
    for i, om in enumerate(oms):
        for j, ol in enumerate(ols):
            mu_model = dist_mod(zz, om, ol)
            for k, m in enumerate(mscr):
                mu_model_norm = mu_model + m
                chi2_test = np.sum((mu_model_norm - mu) ** 2 / mu_error2)
                if chi2_test < chi2[i, j]:
                    chi2[i, j] = chi2_test
                    mscr_used[i, j] = k
        print("Done %d out of %d" % (i, oms.size))
    np.savetxt(saved_data_filename, chi2, fmt="%10.4f")

likelihood = np.exp(-0.5 * chi2)
c = ChainConsumer()
c.add_chain([oms, ols], weights=likelihood, grid=True, parameters=[r"$\Omega_m$", r"$\Omega_\Lambda$"], name="Supernova")
c.plotter.plot(filename="plots/fit.png", figsize="column")
print(c.analysis.get_latex_table())
