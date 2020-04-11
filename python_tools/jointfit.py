
import numpy as np
import sys
import os
from astropy.io import fits
from cosmology import Cosmology
from scipy.integrate import quad, simps
from scipy.interpolate import RectBivariateSpline, InterpolatedUnivariateSpline, interp1d
from scipy.optimize import fsolve

class Model1:
    '''
    Nadathur & Percival RSD model (full)
    '''

    def __init__(self, delta_r_files, xi_r_files, sv_files, xi_smu_files,
                 covmat_file, smins, smaxs, full_fit=1):

        delta_r_files = delta_r_files.split(',')
        xi_r_files = xi_r_files.split(',')
        sv_files = sv_files.split(',')
        xi_smu_files = xi_smu_files.split(',')
        smins = [int(i) for i in smins.split(',')]
        smaxs = [int(i) for i in smaxs.split(',')]

        self.ndenbins = len(delta_r_files)
        delta_r_file = {}
        xi_r_file = {}
        sv_file = {}
        xi_smu_file = {}
        smin = {}
        smax = {}

        for j in range(self.ndenbins):
            delta_r_file['den{}'.format(j)] = delta_r_files[j]
            xi_r_file['den{}'.format(j)] = xi_r_files[j]
            sv_file['den{}'.format(j)] = sv_files[j]
            xi_smu_file['den{}'.format(j)] = xi_smu_files[j]
            smin['den{}'.format(j)] = smins[j]
            smax['den{}'.format(j)] = smaxs[j]
            

        # full fit (monopole + quadrupole)
        self.full_fit = bool(full_fit)

        print("Setting up redshift-space distortions model.")

        # cosmology for Minerva
        self.om_m = 0.285
        self.s8 = 0.828
        self.cosmo = Cosmology(om_m=self.om_m, s8=self.s8)
        self.nmocks = 120

        self.eff_z = 0.57
        self.b = 2.05

        self.growth = self.cosmo.get_growth(self.eff_z)
        self.f = self.cosmo.get_f(self.eff_z)
        self.s8norm = self.s8 * self.growth 

        eofz = np.sqrt((self.om_m * (1 + self.eff_z) ** 3 + 1 - self.om_m))
        self.iaH = (1 + self.eff_z) / (100. * eofz) 

        # read covariance matrix
        if os.path.isfile(covmat_file):
            print('Reading covariance matrix: ' + covmat_file)
            self.cov = np.load(covmat_file)
            self.icov = np.linalg.inv(self.cov)
        else:
            sys.exit('Covariance matrix not found.')

        self.r_for_xi = {}
        self.r_for_delta = {}
        self.r_for_sv = {}
        self.s_for_xi = {}
        self.mu_for_xi = {}
        self.delta_r = {}
        self.Delta_r = {}
        self.xi_r = {}
        self.sv = {}
        self.xi0_s = {}
        self.xi2_s = {}

        self.datavec = np.array([])

        for j in range(self.ndenbins):
            denbin = 'den{}'.format(j)
            # read real-space monopole
            data = np.genfromtxt(xi_r_file[denbin])
            self.r_for_xi[denbin] = data[:,0]
            xi_r = data[:,-1]
            self.xi_r[denbin] = InterpolatedUnivariateSpline(self.r_for_xi[denbin], xi_r, k=3, ext=3)

            # read void-matter correlation function
            data = np.genfromtxt(delta_r_file[denbin])
            self.r_for_delta[denbin] = data[:,0]
            delta_r = data[:,-1]
            self.delta_r[denbin] = InterpolatedUnivariateSpline(self.r_for_delta[denbin], delta_r, k=3, ext=3)

            integral = np.zeros_like(self.r_for_delta[denbin])
            for i in range(len(integral)):
                integral[i] = quad(lambda x: self.delta_r[denbin](x) * x ** 2, 0, self.r_for_delta[denbin][i], full_output=1)[0]
            Delta_r = 3 * integral / self.r_for_delta[denbin] ** 3
            self.Delta_r[denbin] = InterpolatedUnivariateSpline(self.r_for_delta[denbin], Delta_r, k=3, ext=3)

            # read los velocity dispersion profile
            data = np.genfromtxt(sv_file[denbin])
            self.r_for_sv[denbin] = data[:,0]
            sv = data[:,-1] / data[-1, -1]
            self.sv[denbin ]= InterpolatedUnivariateSpline(self.r_for_sv[denbin], sv, k=3, ext=3)

            # read redshift-space correlation function
            self.s_for_xi[denbin], self.mu_for_xi[denbin], xi_smu_obs = self.readCorrFile(xi_smu_file[denbin])
            s, self.xi0_s[denbin] = self._getMonopole(self.s_for_xi[denbin], self.mu_for_xi[denbin], xi_smu_obs)
            s, self.xi2_s[denbin] = self._getQuadrupole(self.s_for_xi[denbin], self.mu_for_xi[denbin], xi_smu_obs)

            # restrict measured vectors to the desired fitting scales
            scales = (s >= smin[denbin]) & (s <= smax[denbin])

            self.s_for_xi[denbin] = self.s_for_xi[denbin][scales]
            self.r_for_xi[denbin] = self.r_for_xi[denbin][scales]
            self.r_for_delta[denbin] = self.r_for_delta[denbin][scales]
            self.r_for_sv[denbin] = self.r_for_sv[denbin][scales]
            self.xi0_s[denbin] = self.xi0_s[denbin][scales]
            self.xi2_s[denbin] = self.xi2_s[denbin][scales]

            if self.full_fit:
                self.datavec = np.concatenate((self.datavec, self.xi0_s[denbin], self.xi2_s[denbin]))
            else:
                self.datavec = np.concatenate((self.datavec, self.xi2_s[denbin]))

        

    def log_likelihood(self, theta):
        fs8, sigma_v, epsilon = theta
        alpha = 1.0
        alpha_para = alpha * epsilon ** (-2/3)
        alpha_perp = epsilon * alpha_para

        modelvec = np.array([])

        for j in range(self.ndenbins):
            denbin = 'den{}'.format(j)

            xi0, xi2 = self.theory_multipoles(fs8, sigma_v,
                                            alpha_perp, alpha_para,
                                            self.s_for_xi[denbin], self.mu_for_xi[denbin], denbin)

            if self.full_fit:
                modelvec = np.concatenate((modelvec, xi0, xi2))
            else:
                modelvec = np.concatenate((modelvec, xi2))

        chi2 = np.dot(np.dot((modelvec - self.datavec), self.icov), modelvec - self.datavec)
        loglike = -self.nmocks/2 * np.log(1 + chi2/(self.nmocks-1))
        return loglike

    def log_prior(self, theta):
        fs8, sigma_v, epsilon = theta


        if 0.1 < fs8 < 2.0 and 50 < sigma_v < 500 and 0.8 < epsilon < 1.2:
            return 0.0
        
        return -np.inf


    def theory_multipoles(self, fs8, sigma_v, alpha_perp, alpha_para, s, mu, denbin):

        monopole = np.zeros(len(s))
        quadrupole = np.zeros(len(s))
        true_mu = np.zeros(len(mu))
        xi_model = np.zeros(len(mu))
        scaled_fs8 = fs8 / self.s8norm

        # rescale input monopole functions to account for alpha values
        mus = np.linspace(0, 1., 101)
        r = self.r_for_delta[denbin]
        rescaled_r = np.zeros_like(r)
        for i in range(len(r)):
            rescaled_r[i] = np.trapz((r[i] * alpha_para) * np.sqrt(1. + (1. - mus ** 2) *
                            (alpha_perp ** 2 / alpha_para ** 2 - 1)), mus)

        x = rescaled_r
        y1 = self.xi_r[denbin](r)
        y2 = self.delta_r[denbin](r)
        y3 = self.Delta_r[denbin](r)
        y4 = self.sv[denbin](r)

        # build rescaled interpolating functions using the relabelled separation vectors
        rescaled_xi_r = InterpolatedUnivariateSpline(x, y1, k=3)
        rescaled_delta_r = InterpolatedUnivariateSpline(x, y2, k=3, ext=3)
        rescaled_Delta_r = InterpolatedUnivariateSpline(x, y3, k=3, ext=3)
        rescaled_sv = InterpolatedUnivariateSpline(x, y4, k=3, ext=3)
        sigma_v = alpha_para * sigma_v
 

        for i in range(len(s)):
            for j in range(len(mu)):
                true_sperp = s[i] * np.sqrt(1 - mu[j] ** 2) * alpha_perp
                true_spar = s[i] * mu[j] * alpha_para
                true_s = np.sqrt(true_spar ** 2. + true_sperp ** 2.)
                true_mu[j] = true_spar / true_s

                rpar = true_spar + true_s * scaled_fs8 * rescaled_Delta_r(true_s) * true_mu[j] / 3.
                sy_central = sigma_v * rescaled_sv(np.sqrt(true_sperp**2 + rpar**2)) * self.iaH
                y = np.linspace(-3 * sy_central, 3 * sy_central, 100)

                rpar = true_spar + true_s * scaled_fs8 * rescaled_Delta_r(true_s) * true_mu[j] / 3. - y
                rr = np.sqrt(true_sperp ** 2 + rpar ** 2)
                sy = sigma_v * rescaled_sv(rr) * self.iaH

                integrand = (1 + rescaled_xi_r(rr)) * \
                            (1 + (scaled_fs8 * rescaled_Delta_r(rr) / 3. - y * true_mu[j] / rr) * (1 - true_mu[j]**2) +
                             scaled_fs8 * (rescaled_delta_r(rr) - 2 * rescaled_Delta_r(rr) / 3.) * true_mu[j]**2)
                integrand = integrand * np.exp(-(y**2) / (2 * sy**2)) / (np.sqrt(2 * np.pi) * sy)
                xi_model[j] = np.trapz(integrand, y) - 1


            # build interpolating function for xi_smu at true_mu
            mufunc = InterpolatedUnivariateSpline(true_mu, xi_model, k=3)
            
            # get multipoles
            xaxis = np.linspace(-1, 1, 1000)

            yaxis = mufunc(xaxis) / 2
            monopole[i] = simps(yaxis, xaxis)

            yaxis = mufunc(xaxis) * 5 / 2 * (3 * xaxis**2 - 1) / 2
            quadrupole[i] = simps(yaxis, xaxis)
            
            
        return monopole, quadrupole


    def readCorrFile(self, fname):
        data = np.genfromtxt(fname)
        s = np.unique(data[:,0])
        mu = np.unique(data[:,1])

        xi_smu = np.zeros([len(s), len(mu)])
        counter = 0
        for i in range(len(s)):
            for j in range(len(mu)):
                xi_smu[j, i] = data[counter, -1]
                counter += 1

        return s, mu, xi_smu

    def _getMonopole(self, s, mu, xi_smu):
        monopole = np.zeros(xi_smu.shape[0])
        for i in range(xi_smu.shape[0]):
            mufunc = InterpolatedUnivariateSpline(mu, xi_smu[i, :], k=3)
            xaxis = np.linspace(-1, 1, 1000)
            yaxis = mufunc(xaxis) / 2
            monopole[i] = simps(yaxis, xaxis)
        return s, monopole

    def _getQuadrupole(self, s, mu, xi_smu):
        quadrupole = np.zeros(xi_smu.shape[0])
        for i in range(xi_smu.shape[0]):
            mufunc = InterpolatedUnivariateSpline(mu, xi_smu[i, :], k=3)
            xaxis = np.linspace(-1, 1, 1000)
            yaxis = mufunc(xaxis) * 5 / 2 * (3 * xaxis**2 - 1) / 2
            quadrupole[i] = simps(yaxis, xaxis)

        return s, quadrupole



            



