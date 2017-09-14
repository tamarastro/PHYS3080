import numpy as np
from scipy.integrate import quad


# ----------Define the function we want to integrate----------
def Hz_inverse(z, OM, OL):
    Hz = np.sqrt((1 + z) ** 2 * (OM * z + 1) - OL * z * (z + 2));
    return 1.0 / Hz;


# ----------Function to calculate the distance modulus--------
# ----------correcting for curvature of the universe----------
def dist_mod(zs, om, ol):
    ok = 1.0 - om - ol
    x = np.array([quad(Hz_inverse, 0, z, args=(om, ol))[0] for z in zs])
    if ok < 0.0:
        R0 = 1 / np.sqrt(-ok)
        D = R0 * np.sin(x / R0)

    if ok > 0.0:
        R0 = 1 / np.sqrt(ok)
        D = R0 * np.sinh(x / R0)

    if ok == 0.0:
        D = x

    LumDist = D * (1 + zs)
    DM = 5 * np.log10(LumDist)
    return DM


# ---------- Import data -------------------------------------
data = np.loadtxt("data.txt")
zz = data[:, 0]
mu = data[:, 2]
mu_error = data[:, 3] + 0.02
mu_error2 = mu_error ** 2
# Define cosntants
H0 = 70.0
c_H0 = 3.e5 / H0

# ---------- Set up fitting ranges ---------------------------
n = 6  # Increase this for a finer grid
oms = np.linspace(0, 0.7, n)
ols = np.linspace(0, 1.0, n)
n_marg = 200
mscr_guess = 5.0 * np.log10(c_H0) + 25
mscr = np.linspace(mscr_guess - 0.5, mscr_guess + 0.5, n_marg)

chi2 = np.ones((n, n)) * 9e9  # Our chi2 values, set initially to super large values
mscr_used = np.zeros((n, n))

# ---------- Do the fit ---------------------------
for i, om in enumerate(oms):
    for j, ol in enumerate(ols):
        mu_model = dist_mod(zz, om, ol)
        for k, m in enumerate(mscr):
            mu_model_norm = mu_model + m;
            chi2_test = np.sum((mu_model_norm - mu) ** 2 / mu_error2);
            if chi2_test < chi2[i, j]:
                chi2[i, j] = chi2_test
                mscr_used[i, j] = k
    print("Done %d out of %d" % (i, oms.size))

np.savetxt("myfit.txt", chi2, fmt="%10.4f")  # Back this up or something yo
likelihood = np.exp(-0.5 * chi2)

from chainconsumer import ChainConsumer

c = ChainConsumer()
c.add_chain([oms, ols], weights=likelihood, grid=True, parameters=[r"$\Omega_m$", r"$\Omega_\Lambda$"])
c.plot(filename="output.png", figsize="page")
