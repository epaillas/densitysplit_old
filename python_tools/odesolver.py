import numpy as np
from scipy.integrate import odeint, quad
from scipy.interpolate import InterpolatedUnivariateSpline
import matplotlib.pyplot as plt
import sys
plt.rcParams['text.usetex']=True

# The equation to solve reads
# y'' + (1/2 - 3/2 w om_l) y' + om_m/2 (y^{-3} - 1) y = 0

# let h = y'
# then the equation becomes
# h' + (1/2 - 3/2 w om_l) h + om_m/2 (y^{-3} - 1) y = 0

def TimeIntegral(a):
    return (1./(om_m*a**(-1) + om_l*a**(-1-3*w))**(1/2))

def CosmologicalTime(a):
    H0 = 100
    t = []
    for i in a:
	    ans_i = quad(TimeIntegral, 0, i)[0]
	    t.append((1/(H0))*ans_i)
    t = np.asarray(t)

    t = InterpolatedUnivariateSpline(a, t, k=3)

    return t

def sphe(g, a, om_m0, om_l0, w):
    om_m = om_m0 / (om_m0 + om_l0 * a**(-3*w))
    om_l = om_l0 / (om_m0 * a**(3*w) + om_l0)
    y, h = g
    dgda = [h, -(1/2 - 3/2*om_l*w)*h - om_m/2*(y**(-3) - 1)*y]
    return dgda

def Hubble(a, om_m, om_l):
    return 100 * np.sqrt(om_m * a ** -3 + om_l)


delta_r_file = '/Volumes/BlackIce/density_split/den_cats/Real/TopHatProfiles/20Mpc/\
Galaxies_HOD_z0.57_Real_den2.CCF_DM_monopole'

vr_file = '/Volumes/BlackIce/density_split/den_cats/Real/TopHatProfiles/20Mpc/\
Galaxies_HOD_z0.57_Real_den2.CCF_gal_vr'

z = np.linspace(10000, 0.57, 1000)
a = 1/(1 + z)

# read void-matter correlation function
data = np.genfromtxt(delta_r_file)
r_for_delta = data[:,0]
delta_r = data[:,-2]
delta_r = InterpolatedUnivariateSpline(r_for_delta, delta_r, k=3, ext=3)

integral = np.zeros_like(r_for_delta)
for i in range(len(integral)):
    integral[i] = quad(lambda x: delta_r(x) * x ** 2, 0, r_for_delta[i], full_output=1)[0]
Delta_r = 3 * integral / r_for_delta ** 3
Deltap_r = np.gradient(Delta_r)
Delta_r = InterpolatedUnivariateSpline(r_for_delta, Delta_r, k=3, ext=3)
Deltap_r = InterpolatedUnivariateSpline(r_for_delta, Deltap_r, k=3, ext=3)

data = np.genfromtxt(vr_file)
r_for_vr = data[:,0]
vr = data[:,-2]
vr = InterpolatedUnivariateSpline(r_for_vr, vr, k=3, ext=3)


om_m = 0.285
om_l = 1 - om_m
w = -1
t = CosmologicalTime(a)(a)

delta_li = np.linspace(-1,2.5, 1000)

sols1 = []
sols2 = []

for dl in delta_li:

    g0 = [1 - dl/3, -dl/3]

    sol = odeint(sphe, g0, t, args=(om_m, om_l,w))
    y = sol[:,0]
    yprime = sol[:,1]

    sols1.append(y[-1]**-3 - 1)
    sols2.append(yprime[-1])

sols1 = np.asarray(sols1)
sols2 = np.asarray(sols2)

interp_den = InterpolatedUnivariateSpline(sols1, delta_li, k=3)
interp_vel = InterpolatedUnivariateSpline(delta_li, sols2, k=3)

matched_dls = interp_den(Delta_r(r_for_delta))
matched_vpecs = interp_vel(matched_dls)

H = Hubble(a=a[-1], om_m=om_m, om_l=om_l)
q = r_for_delta * a[-1] * (1 + Delta_r(r_for_delta))**(1/3) 
vpec = matched_vpecs * 100 * q

fig, ax = plt.subplots(figsize=(6,5))
f = 0.7596841096514576
vpec_lin = -1/3 *f  * H * r_for_delta * a[-1]* Delta_r(r_for_delta)
ax.plot(r_for_vr, vr(r_for_vr), label='Measurement',marker='o', mfc='none', mew=1.5, ms=6.0, ls='',
color='k', alpha=0.7)
ax.plot(r_for_delta, vpec, label='Solving ODE')
ax.plot(r_for_delta, vpec_lin, label='Linear approx.')
ax.set_xlabel(r'$r\ [\mathrm{Mpc}/h]$', fontsize=17)
ax.set_ylabel(r'$v_{\mathrm{pec}}\ [\mathrm{km/s}]$', fontsize=17)
ax.tick_params(labelsize=13)
ax.legend(fontsize=17)
plt.tight_layout()
plt.show()