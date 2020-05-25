import numpy as np
from scipy.io import FortranFile
from scipy.interpolate import InterpolatedUnivariateSpline
import matplotlib.pyplot as plt
from scipy.integrate import simps, quad
import click


# read file name from command line
@click.command()
@click.option('--gal_den_spi', type=str, required=True)
@click.option('--filter_file', type=str, required=True)
@click.option('--handle', type=str, required=True)
@click.option('--ndenbins', type=int, required=True)

def split_densities(gal_den_spi,
                    filter_file,
                    handle,
                    ndenbins):


    print('\nSplitting densities for the following arguments:')
    print('gal_den_spi: {}'.format(gal_den_spi))
    print('filter_file: {}'.format(filter_file))
    print('handle: {}'.format(handle))

    # open galaxy density (r-mu bins)
    f = FortranFile(gal_den_spi, 'r')
    ncentres = f.read_ints()[0]
    nrbins = f.read_ints()[0]
    nmubins = f.read_ints()[0]
    print('ncentres, nrbins, nmubins = ({}, {}, {})'.format(ncentres, nrbins, nmubins))

    # read raw data and close file
    rbin = f.read_reals(dtype=np.float64).T
    mubin = f.read_reals(dtype=np.float64).T
    xi_spi = f.read_reals(dtype=np.float64).reshape(nmubins* nrbins, ncentres).T
    f.close()

    # open filter file
    f = FortranFile(filter_file, 'r')
    ncentres = f.read_ints()[0]
    print('ncentres: {}'.format(ncentres))
    smoothed_delta = f.read_reals(dtype=np.float64)
    idx = np.argsort(smoothed_delta)

    # sort arrays
    xi_spi = xi_spi[idx]

    sorted_profiles = np.c_[xi_spi]
    
    profiles = {}
    # divide profiles in quantiles
    for i in range(1, ndenbins + 1):
        profiles['den{}'.format(i)] = sorted_profiles[int((i-1)*ncentres/ndenbins):int(i*ncentres/ndenbins)]

    # get mean split profiles
    xi_spi_mean = {}

    for i in range(1, ndenbins + 1):
        den = profiles['den{}'.format(i)]
        xi_spi_mean['den{}'.format(i)] = np.mean(den, axis=0)

    # save split profiles to file
    bins = np.zeros([nrbins * nmubins, 2])
    count = 0
    for j in range(nmubins):
        for i in range(nrbins):

            bins[count, 0] = rbin[i]
            bins[count, 1] = mubin[j]
            count += 1


    for i in range(1, ndenbins + 1):
        xi_spi_file = handle + '_den{}'.format(i) + '.CCF_gal_spi'
        xi_spi_out = np.c_[bins, xi_spi_mean['den{}'.format(i)]]

        fmt = 3*'%10.3f '
        np.savetxt(xi_spi_file, xi_spi_out, fmt=fmt)

        
if __name__ == '__main__':
    split_densities()




