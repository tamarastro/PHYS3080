import numpy as np
from scipy.integrate import quad


def Hz_inverse(z,OM,OL):
    Hz = np.sqrt((1+z)**2 * (OM*z+1)-OL*z*(z+2));
    return 1.0/Hz;
        
 
def dist_mod(zs, om, ol):
    ok = 1.0 - om - ol
    x = np.array([quad(Hz_inverse, 0, z, args=(om, ol))[0] for z in zs])
    if ok < 0.0:
        R0 = 1/np.sqrt(-ok)
        D = R0*np.sin(x/R0)
    
    if ok > 0.0:
        R0 = 1/np.sqrt(ok);
        D = R0*np.sinh(x/R0)
    
    if ok == 0.0:
        D = x

    LumDist = D*(1+zs);
    DM = 5*np.log10(LumDist);
    return DM
    
# Import data
data = np.genfromtxt("jla_lcparams_simple_new.txt",names=True,comments='#',dtype=None,skip_header=14)
zz = data['zcmb']
mu = data['mb']
mu_error = data['dmb']  
mu_error2 = mu_error**2 # squared for ease of use later
ids = data['set'] # Source ID [3 = Nearby; 2 = SDSS; 1 = SNLS; 4 = HST]

# Define cosntants
H0 = 70.0
c_H0 = 3.e5/H0 

# Set up fitting ranges
n = 21 # Increase this for a finer grid
oms = np.linspace(0, 0.7, n)
ols = np.linspace(0, 1.0, n)
n_marg = 200
mscr_guess = 43.1 # Approx 5.0 * np.log10(c_H0) + 25 - 19
mscr = np.linspace(mscr_guess - 0.5, mscr_guess + 0.5, n_marg)

chi2 = np.ones((n,n)) * np.infty # Our chi2 values, set initially to super large values
mscr_used = np.zeros((n,n))

# Do the fit
for i, om in enumerate(oms):
	for j, ol in enumerate(ols):
		mu_model = dist_mod(zz,om,ol)
		for k, m in enumerate(mscr):
			mu_model_norm = mu_model + m;
			chi2_test = np.sum( (mu_model_norm - mu)**2 / mu_error2 );
			if chi2_test < chi2[i,j]:
				chi2[i,j] = chi2_test
				mscr_used[i,j] = k
	print("Done %d out of %d" % (i, oms.size))
	
np.savetxt("myfit.txt", chi2, fmt="%10.4f") # Back this up or something yo

likelihood = np.exp(-0.5 * (chi2-chi2.min()))

print('min chi2 = ',chi2.min()/(zz.size - 3))
ijbest = np.unravel_index(np.argmin(chi2),(n,n))
print('Best fit (Omega_M,Omega_Lambda) = (',oms[ijbest[0]],',' ,ols[ijbest[1]],')' )

from chainconsumer import ChainConsumer
c = ChainConsumer()
c.add_chain([oms, ols], weights=likelihood, grid=True, parameters=[r"$\Omega_m$", r"$\Omega_\Lambda$"])
c.plot(filename="output.png", figsize="page")