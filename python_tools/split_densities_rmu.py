import numpy as np
from scipy.io import FortranFile
from scipy.interpolate import InterpolatedUnivariateSpline
import matplotlib.pyplot as plt
from scipy.integrate import simps, quad
import click


# read file name from command line
@click.command()
@click.option('--gal_den_monopole', type=str, required=True)
@click.option('--gal_den_rmu', type=str, required=True)
@click.option('--handle', type=str, required=True)

def split_densities(gal_den_monopole, gal_den_rmu, handle):


    print('\nSplitting densities for the following arguments:')
    print('gal_den_monopole: {}'.format(gal_den_monopole))
    print('gal_den_rmu: {}'.format(gal_den_rmu))
    print('handle: {}'.format(handle))

    # open galaxy density (r-mu bins)
    f = FortranFile(gal_den_rmu, 'r')
    ncentres = f.read_ints()[0]
    nrbins = f.read_ints()[0]
    nmubins = f.read_ints()[0]
    print('ncentres, nrbins, nmubins = ({}, {}, {})'.format(ncentres, nrbins, nmubins))

    # read raw data and close file
    rbin = f.read_reals(dtype=np.float64).T
    mubin = f.read_reals(dtype=np.float64).T
    dd = f.read_reals(dtype=np.float64).reshape(nmubins* nrbins,  ncentres).T
    xi_rmu = f.read_reals(dtype=np.float64).reshape(nmubins* nrbins, ncentres).T
    xibar_rmu = f.read_reals(dtype=np.float64).reshape(nmubins*nrbins, ncentres).T
    f.close()

    # open galaxy density (monopole)
    f = FortranFile(gal_den_monopole, 'r')
    ncentres = f.read_ints()[0]
    nrbins = f.read_ints()[0]
    print('ncentres, nrbins = ({}, {})'.format(ncentres, nrbins))

    # read raw data and close file
    rbin = f.read_reals(dtype=np.float64)
    dd = f.read_reals(dtype=np.float64).reshape(nrbins, ncentres).T
    xi_r = f.read_reals(dtype=np.float64).reshape(nrbins, ncentres).T
    xibar_r = f.read_reals(dtype=np.float64).reshape(nrbins, ncentres).T
    f.close()

    # get indices of sorted array
    sorted_idx = np.argsort(xibar_r[:, 10])

    # sort arrays
    xi_rmu = xi_rmu[sorted_idx]

    sorted_profiles = np.c_[xi_rmu]
    
    profiles = {}
    # divide profiles in quantiles
    for i in range(1, 6):
        profiles['den{}'.format(i)] = sorted_profiles[int((i-1)*ncentres/5):int(i*ncentres/5)]

    # get mean split profiles
    xi_rmu_mean = {}

    for i in range(1, 6):
        den = profiles['den{}'.format(i)]
        xi_rmu_mean['den{}'.format(i)] = np.mean(den, axis=0)

    # save split profiles to file
    bins = np.zeros([nrbins * nmubins, 2])
    count = 0
    for i in range(nrbins):
        for j in range(nmubins):

            bins[count, 0] = rbin[i]
            bins[count, 1] = mubin[j]
            count += 1


    for i in range(1,6):
        xi_rmu_file = handle + '_den{}'.format(i) + '.CCF_rmu'
        xi_rmu_out = np.c_[bins, xi_rmu_mean['_den{}'.format(i)]]

        fmt = 3*'%10.3f '
        np.savetxt(xi_rmu_file, xi_rmu_out, fmt=fmt)

        
if __name__ == '__main__':
    split_densities()




