import matplotlib.pyplot as plt
import numpy as np
import astropy.units as u
from scipy.integrate import quad

def Salpeter(mass):
	value, skip = quad(lambda x: x**-1.35, mass[0], mass[-1])
	xi = 1./value
	eta = xi * mass**-1.35
	return eta, xi

def func_core(alpha, M_star, phi_star, M):
	A = 0.4*np.log(10)
	B = np.power(10.,(0.4*(alpha + 1)*(M_star - M)))
	C = -np.power(10.,(0.4*(M_star - M)))
	value = A*phi_star*B*np.exp(C)
	return value

def Schechter(alpha, M_star, phi_star):
	xaxis_M = np.arange(-26., -17., 0.1)
	yaxis_phi = func_core(alpha, M_star, phi_star, xaxis_M)
	return xaxis_M, yaxis_phi

test_phi_star = 1.66*1e-2*(1/(u.Mpc * u.Mpc * u.Mpc))
test_alpha = -0.81
test_M_star = -23.19
data_x, data_y = Schechter(test_alpha, test_M_star, test_phi_star)
virgo_x, virgo_y = Schechter(-1.35, test_M_star, test_phi_star)
other_x, other_y = Schechter(-0.6, test_M_star, test_phi_star)


fig = plt.figure()
fig.suptitle("Galaxy Luminosity Function")
plt.semilogy(data_x, data_y, color='black', linestyle='-', label=r'Schechter $\alpha = -0.81$')
plt.semilogy(virgo_x, virgo_y, color = "r", linestyle = "--", label=r"Virgo Cluster $\alpha = -1.35$")
plt.semilogy(other_x, other_y, color = "b", linestyle = ":", label=r"Others $\alpha = -0.6$")
plt.xlabel(r"$M_k$-5log h")
plt.ylabel(r"$\phi$/($h^3$ $Mpc^-3$ $mag^-1$)")
plt.legend(loc = "lower right")
plt.savefig("Schechter_Function.png")

mass_bin = np.linspace(0.1, 120., 50)
eta_bin, xi = Salpeter(mass_bin)
frac_1, skip = quad(lambda x: xi * x**-1.35, 1, 120)
frac_2, skip = quad(lambda x: xi * x**-1.35, .1, 120)
print("xi: %f\nfrac > 1: %.3f\nfrac total: %.3f\n"%(xi, frac_1, frac_2))