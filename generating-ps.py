import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate
from scipy import integrate
import camb

#Defining the Cosmology
omegab,omegac,NS,tau,H0,As = [0.022383,0.12011,0.96605,0.0543,
70.0,2.1005829617*1e-9]

KMAX = 10
kMIN = 1e-2
NPOINTS = 100

#CAMB parameters

#Linear
params1 = camb.CAMBparams()
params1.set_cosmology(ombh2 = omegab, omch2 = omegac, omk = 0, 
H0 = H0, TCMB = 2.7255, tau=tau)
params1.WantTransfer = True
params1.NonLinear = 'NonLinear_none'
params1.set_matter_power(kmax=KMAX)
params1.InitPower.set_params(ns = NS,As = As)

#Halo Fit
params2 = camb.CAMBparams()
params2.set_cosmology(ombh2 = omegab, omch2 = omegac, omk = 0, 
H0 = H0, TCMB = 2.7255, tau=tau)
params2.WantTransfer = True
params2.NonLinear = 'NonLinear_both'
params2.set_matter_power(kmax=KMAX)
params2.InitPower.set_params(ns = NS,As = As)
params2.InitPower.set_params(ns = NS,As = As)


results1 = camb.get_results(params1)
results2 = camb.get_results(params2)


#Asking for the Power Spectrum
print("Calculating the Power Spectrum...\n")
kh,z,pk = results1.get_matter_power_spectrum(minkh=kMIN, maxkh=KMAX, npoints = NPOINTS)
pk2 = results2.get_matter_power_spectrum(minkh=kMIN, maxkh=KMAX, npoints = NPOINTS)[2]
print("Done!")

#Saving Data
p_lin = pk[0,:]
p_nlin = pk2[0,:]
data = np.column_stack((kh,p_lin,p_nlin))
np.savetxt("/home/thiago/Desktop/PS.txt",data)

#Plot
plt.loglog(kh,p_lin,label="Linear")
plt.loglog(kh,p_nlin,label = "Halo-Fit")
plt.legend(loc='best')
plt.ylabel("$P_{lin}$")
plt.xlabel("K (h/Mpc)")
plt.grid()
plt.show()
