import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate
from scipy import integrate
import camb
from tqdm import tqdm

#Reading a PS
data = np.loadtxt("PS.txt")
kh = data[:,0]
p_lin = data[:,1]
p_nlin = data[:,2]
p_linear = interpolate.interp1d(kh,p_lin)

#Integrand of P13
def f13(q,k):
	r = np.exp(q)/k
	#F=-42.0*r**2.0+100.0-158.0/r**2.0+12.0/r**4.0+3.0/r**5.0*(r**2.0-1.0)**3.0*(2.0+7.0*r**2.0)*np.log((r+1.0)/np.abs(r-1.0))
	#return 1.0/252.0/(2.0*np.pi)**2.0*p_linear(np.exp(q))*np.exp(3*q)*F
	if r < (100.0):
		F=-42.0*r**2.0+100.0-158.0/r**2.0+12.0/r**4.0+3.0/r**5.0*(r**2.0-1.0)**3.0*(2.0+7.0*r**2.0)*np.log((r+1.0)/np.abs(r-1.0))
		return 1.0/252.0/(2.0*np.pi)**2.0*p_linear(np.exp(q))*np.exp(3*q)*F
	else:
		return -61.0/210.0/(6*np.pi**2)*np.exp(3*q)/r**2*p_linear(np.exp(q))

#Limits of the angular integral in P22
def mulow(x):
  return max(-1.0,(kh[-1]**2.0-kk**2.0-np.exp(x)**2.0)/(-2.0*kk*np.exp(x)))

def muhigh(x):
  return min(1.0,(kh[0]**2.0-kk**2.0-np.exp(x)**2.0)/(-2.0*kk*np.exp(x)))


#Integrand of P22
def upper_mu(x):
	return min(1.0,(kk**2 + np.exp(2*x))/(2*kk*np.exp(x)))

def lower_mu(x):
	return max(-1.0,-(kk**2+np.exp(x))/(2*kk*np.exp(x)))


def f22(mu,q,k):
	r = np.exp(q)/k
	F = (7.0*mu+(3.0-10.0*mu**2)*r)/(14.0*r*(r**2-2.0*mu*r+1.0))
	psik = (k**2+np.exp(2*q)-2.0*k*mu*np.exp(q))**0.5
	
	if (psik>kh[0] and psik<kh[-1]):
		return 1.0/2.0/np.pi**2.0*np.exp(3*q)*p_linear(np.exp(q))*p_linear(psik)*F**2
	else:
		return 0
	
	'''
	if (r<100.0):
		if (psik>kh[0] and psik<kh[-1]):
			return 1.0/2.0/np.pi**2.0*np.exp(3*q)*p_linear(np.exp(q))*p_linear(psik)*F**2
		else:
			return 0
	else:
		return (9.0/(98.0*2.0*np.pi**2))*p_linear(np.exp(q))**2*np.exp(3*q)/r**4
	'''
P22 = np.zeros_like(kh)
P13 = np.zeros_like(kh)

for i in tqdm(range(0,P22.shape[0])):
	kk = kh[i]
	P22[i] = integrate.dblquad(f22,np.log(kh[0]),np.log(kh[-1]),lower_mu,upper_mu,args=(kh[i],),epsrel=1e-4,epsabs=50)[0]
	p13 = integrate.quad(f13,np.log(kh[0]),np.log(kh[-1]),args=(kh[i],),limit=200,epsrel=1e-4)
	P13[i] = p_linear(kh[i])*p13[0]

result = np.column_stack((kh,P22,P13,p_lin,p_nlin))
np.savetxt("/home/thiago/Desktop/1-loop10.txt",result)

plt.loglog(kh,abs(P13),'--',label="$|P_{13}|$")
plt.loglog(kh,P22,'--',label="$P_{22}$")
plt.loglog(kh,abs(P22+P13),label="$|P_{22}+P_{13}|$")
plt.loglog(kh,p_lin, label="$P_{linear}$")
plt.loglog(kh,p_lin+P22+P13,'.',label="$P_{1-loop}$")
plt.loglog(kh,p_nlin,':',label="$P_{non-linear}$")
#plt.xlim((kh[0],2))
plt.ylim((1,10e4))
plt.grid()
plt.legend(loc='best')
plt.savefig("1-loop10.png")
plt.show()	









			
