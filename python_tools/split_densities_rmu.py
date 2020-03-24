import numpy as np
from scipy.io import FortranFile
from scipy.interpolate import InterpolatedUnivariateSpline
import matplotlib.pyplot as plt
import click

def get_xi_smu(data, nrbins, nmubins):

        xi_smu = np.zeros([nrbins, nmubins])
        counter = 0
        for j in range(nrbins):
            for i in range(nmubins):
                xi_smu[i, j] = data[counter]
                counter += 1

        return xi_smu


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
    gal_dd = f.read_reals(dtype=np.float64).reshape(nmubins*nrbins,  ncentres)
    gal_delta = f.read_reals(dtype=np.float64).reshape(nmubins*nrbins, ncentres)
    gal_cumdelta = f.read_reals(dtype=np.float64).reshape(nmubins*nrbins, ncentres)

    f.close()

    profiles = [np.c_[gal_delta[:,i], gal_cumdelta[:,i]] for i in range(ncentres)]
    profiles = np.asarray(profiles)

    # sort profiles according to Delta(r=20mpc/h)
    sorted_profiles = profiles[np.argsort(profiles[:, 10, -1])]

    # divide profiles by their Delta(r=20mpc/h)
    profiles = {}
    for i in range(1, 6):
        profiles['den{}'.format(i)] = sorted_profiles[int((i-1)*ncentres/5):int(i*ncentres/5)]

    # get average densities
    gal_delta = {}

    for i in range(1, 6):
        den = profiles['den{}'.format(i)]
        gal_delta['den{}'.format(i)] = np.mean(den, axis=0)[:,-2]

    for i in range(1,6):
        gal_delta_file = handle + '_den{}'.format(i) + '.CCF_gal_rmu'

        cout = np.zeros([nrbins*nmubins, 3])

        fmt = 3*'%10.3f '

        counter = 0
        for jj in range(nrbins):
            for ii in range(nmubins):
                cout[counter, 0] = rbin[jj]
                cout[counter, 1] = mubin[ii]
                cout[counter, 2] = gal_delta['den{}'.format(i)][counter]
                
                counter += 1

        np.savetxt(gal_delta_file, cout, fmt=fmt)

        
if __name__ == '__main__':
    split_densities()




