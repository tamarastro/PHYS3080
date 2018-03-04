import numpy as np
from scipy.integrate import quad
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt


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
#data = np.loadtxt("jla_lcparams_simple.txt", comments='#', usecols=(1,2,3,4))
data = np.genfromtxt("data/jla_lcparams_simple_new.txt",names=True,comments='#',dtype=None,skip_header=11)
zz = data['zcmb']
mu = data['mb']
mu_error = data['dmb'] 
mu_error2 = mu_error**2 # squared for ease of use later
ids = data['set'] # Source ID [3 = Nearby; 2 = SDSS; 1 = SNLS; 4 = HST]
idnums  = [3, 2, 1, 4]
idnames = ['Nearby','SDSS','SNLS','HST']

#Calculate a model prediction
zz_model = np.logspace(-2,0.2,50) #Make logarithmic redshift array to better sample low-z
mu_model=dist_mod(zz_model,0.3,0.7) #Calculate the distance modulus corresponding to the model redshifts

# To calculate the arbitrary offset we could use:
# total((mu_data-mu_theory)/muerr^2)/total(1./muerr^2)
# But for now I'm just going to set it manually
mscr=43.1
mu_model = mu_model+mscr 



for i, idnum in enumerate(idnums):
    subset = np.where(ids==idnum)
    plt.errorbar(zz[subset],mu[subset],yerr=mu_error[subset],fmt='.',label=idnames[i])
plt.plot(zz_model,mu_model,color='black')
plt.legend(loc='lower right',frameon=False)
plt.xlim(0,1.4)
plt.xlabel('redshift')
plt.ylabel('magnitude')
plt.savefig("plots/hubble_diagram.png")


