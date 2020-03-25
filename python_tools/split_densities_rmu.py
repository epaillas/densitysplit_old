import numpy as np
from scipy.io import FortranFile
from scipy.interpolate import InterpolatedUnivariateSpline
import matplotlib.pyplot as plt
from scipy.integrate import simps, quad
import click

def get_xi_smu(data, nrbins, nmubins):

        xi_smu = np.zeros([nrbins, nmubins])
        counter = 0
        for j in range(nrbins):
            for i in range(nmubins):
                #xi_smu[i, j] = data[counter]
                xi_smu[i, j] = data[counter]
                counter += 1

        return xi_smu

def getMonopole(s, mu, xi_smu):
        monopole = np.zeros(xi_smu.shape[0])
        for i in range(xi_smu.shape[0]):
            mufunc = InterpolatedUnivariateSpline(mu, xi_smu[i, :], k=3)
            xaxis = np.linspace(-1, 1, 1000)
            yaxis = mufunc(xaxis) / 2
            monopole[i] = simps(yaxis, xaxis)
        return s, monopole

def getQuadrupole(s, mu, xi_smu):
        quadrupole = np.zeros(xi_smu.shape[0])
        for i in range(xi_smu.shape[0]):
            mufunc = InterpolatedUnivariateSpline(mu, xi_smu[i, :], k=3)
            xaxis = np.linspace(-1, 1, 1000)
            yaxis = mufunc(xaxis) * 5 / 2 * (3 * xaxis**2 - 1) / 2
            quadrupole[i] = simps(yaxis, xaxis)
        return s, quadrupole


# read file name from command line
@click.command()
@click.option('--gal_den_file', type=str, required=True)
@click.option('--handle', type=str, required=True)

def split_densities(gal_den_file, handle):


    print('\nSplitting densities for the following arguments:')
    print('gal_den_file: {}'.format(gal_den_file))
    print('handle: {}'.format(handle))

    # open galaxy density file and get dimensions
    f = FortranFile(gal_den_file, 'r')
    ncentres = f.read_ints()[0]
    nrbins = f.read_ints()[0]
    nmubins = f.read_ints()[0]
    print('ncentres, nrbins, nmubins = ({}, {}, {})'.format(ncentres, nrbins, nmubins))

    # read raw data and close file
    rbin = f.read_reals(dtype=np.float64)
    mubin = f.read_reals(dtype=np.float64)
    gal_dd = f.read_reals(dtype=np.float64).reshape(nmubins* nrbins,  ncentres)
    gal_delta = f.read_reals(dtype=np.float64).reshape(nmubins* nrbins, ncentres)
    gal_cumdelta = f.read_reals(dtype=np.float64).reshape(nmubins*nrbins, ncentres)

    f.close()

    gal_xi0 = np.zeros([nrbins, ncentres])
    gal_xi2 = np.zeros([nrbins, ncentres])
    gal_xibar = np.zeros([nrbins, ncentres])

    for i in range(ncentres):
        xi_smu = get_xi_smu(data=gal_delta[:, i], nrbins=nrbins, nmubins=nmubins)

        rbin, xi0 = getMonopole(rbin, mubin, xi_smu)
        rbin, xi2 = getQuadrupole(rbin, mubin, xi_smu)
        xi0_func = InterpolatedUnivariateSpline(rbin, xi0, k=3)

        integral = np.zeros_like(rbin)
        for ii in range(len(integral)):
            xaxis = np.linspace(0, rbin[ii], 1000)
            yaxis = xi0_func(xaxis) * xaxis ** 2
            integral[ii] = simps(yaxis, xaxis)
            #integral[ii] = quad(lambda x: xi0_func(x) * x ** 2, 0, rbin[ii], full_output=1)[0]
        xibar = 3 * integral / rbin ** 3

        gal_xi0[:, i] = xi0
        gal_xi2[:, i] = xi2
        gal_xibar[:, i] = xibar


    profiles = [np.c_[gal_xi0[:,i], gal_xi2[:,i], gal_xibar[:,i]] for i in range(ncentres)]
    profiles = np.asarray(profiles)

    # sort profiles according to Delta(r=20mpc/h)
    sorted_profiles = profiles[np.argsort(profiles[:, 10, -1])]

    # divide profiles by their Delta(r=20mpc/h)
    profiles = {}
    for i in range(1, 6):
        profiles['den{}'.format(i)] = sorted_profiles[int((i-1)*ncentres/5):int(i*ncentres/5)]

    # get average densities
    gal_xi0 = {}
    gal_xi2 = {}

    for i in range(1, 6):
        den = profiles['den{}'.format(i)]
        gal_xi0['den{}'.format(i)] = np.mean(den, axis=0)[:,-3]
        gal_xi2['den{}'.format(i)] = np.mean(den, axis=0)[:,-2]


    for i in range(1,6):
        gal_xi0_file = handle + '_den{}'.format(i) + '.CCF_monopole'
        gal_xi2_file = handle + '_den{}'.format(i) + '.CCF_quadrupole'

        xi0_out = np.c_[rbin, gal_xi0['den{}'.format(i)]]
        xi2_out = np.c_[rbin, gal_xi2['den{}'.format(i)]]

        fmt = 2*'%10.3f '

        np.savetxt(gal_xi0_file, xi0_out, fmt=fmt)
        np.savetxt(gal_xi2_file, xi2_out, fmt=fmt)

        
if __name__ == '__main__':
    split_densities()




