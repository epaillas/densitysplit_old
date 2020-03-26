import numpy as np
from scipy.io import FortranFile
import click
import matplotlib.pyplot as plt


# read file name from command line
@click.command()
@click.option('--gal_den_file', type=str, required=True)
@click.option('--dm_den_file', type=str, required=True)
@click.option('--has_velocity', type=bool, required=True)
@click.option('--handle', type=str, required=True)

def split_densities(gal_den_file,
                         dm_den_file,
                         has_velocity,
                         handle):


    print('\nSplitting densities for the following arguments:')
    print('gal_den_file: {}'.format(gal_den_file))
    print('dm_den_file: {}'.format(dm_den_file))
    print('has_velocity: {}'.format(has_velocity))
    print('handle: {}'.format(handle))

    # open galaxy density file and get dimensions
    f = FortranFile(gal_den_file, 'r')
    ncentres = f.read_ints()[0]
    nbins = f.read_ints()[0]
    print('ncentres, nbins = ({}, {})'.format(ncentres, nbins))

    # read raw data and close file
    rbin = f.read_reals()
    gal_dd = f.read_reals(dtype=np.float64).reshape(nbins, ncentres)
    gal_delta = f.read_reals(dtype=np.float64).reshape(nbins, ncentres)
    gal_cumdelta = f.read_reals(dtype=np.float64).reshape(nbins, ncentres)

    if has_velocity:
        gal_vr = f.read_reals(dtype=np.float64).reshape(nbins, ncentres)
        gal_sv_los = f.read_reals(dtype=np.float64).reshape(nbins, ncentres)

    f.close()

    # open DM density file
    f = FortranFile(dm_den_file, 'r')
    ncentres = f.read_ints()[0]
    nbins = f.read_ints()[0]
    print('ncentres, nbins = ({}, {})'.format(ncentres, nbins))

    rbin = f.read_reals()
    dm_dd = f.read_reals(dtype=np.float64).reshape(nbins, ncentres)
    dm_delta = f.read_reals(dtype=np.float64).reshape(nbins, ncentres)
    dm_cumdelta = f.read_reals(dtype=np.float64).reshape(nbins, ncentres)

    if has_velocity:
        dm_vr = f.read_reals(dtype=np.float64).reshape(nbins, ncentres)
        dm_sv_los = f.read_reals(dtype=np.float64).reshape(nbins, ncentres)

    f.close()


    # assemble individual profiles into an array
    if has_velocity:
        profiles = [np.c_[gal_dd[:,i], gal_delta[:,i], gal_cumdelta[:,i],
                    dm_delta[:,i], dm_cumdelta[:,i],
                    gal_vr[:,i], gal_sv_los[:,i]] for i in range(ncentres)]
    else:
        profiles = [np.c_[gal_dd[:,i], gal_delta[:,i], gal_cumdelta[:,i],
                    dm_cumdelta[:,i]]
                    for i in range(ncentres)]

    profiles = np.asarray(profiles)

    # sort profiles according to Delta(r=20mpc/h)
    sorted_profiles = profiles[np.argsort(profiles[:, 10, -5])]

    # divide profiles by their Delta(r=20mpc/h)
    profiles = {}
    for i in range(1, 6):
        profiles['den{}'.format(i)] = sorted_profiles[int((i-1)*ncentres/5):int(i*ncentres/5)]

    # get average densities
    gal_delta = {}
    dm_delta = {}
    gal_cumdelta = {}
    dm_cumdelta = {}

    for i in range(1, 6):
        den = profiles['den{}'.format(i)]
        gal_delta['den{}'.format(i)] = np.mean(den, axis=0)[:,-6]
        gal_cumdelta['den{}'.format(i)] = np.mean(den, axis=0)[:,-5]
        dm_delta['den{}'.format(i)] = np.mean(den, axis=0)[:,-4]
        dm_cumdelta['den{}'.format(i)] = np.mean(den, axis=0)[:,-3]

    # get average vr
    gal_vr = {}
    gal_sv_los = {}

    for i in range(1,6):
        gal_vr['den{}'.format(i)] = np.zeros(nbins)
        gal_sv_los['den{}'.format(i)] = np.zeros(nbins)
        den = profiles['den{}'.format(i)]

        for j in range(len(den)):
            weights = den[j, :, 0]
            gal_vr['den{}'.format(i)] += den[j, :, -2] * weights
            gal_sv_los['den{}'.format(i)] += den[j, :, -1] * weights

        sumweights = np.sum(den, axis=0)[:,0]
        gal_vr['den{}'.format(i)] /= sumweights
        gal_sv_los['den{}'.format(i)] /= sumweights


    for i in range(1,6):
        dm_delta_file = handle + '_den{}'.format(i) + '.CCF_DM_monopole'
        gal_delta_file = handle + '_den{}'.format(i) + '.CCF_gal_monopole'
        gal_vr_file = handle + '_den{}'.format(i) + '.CCF_gal_vr'
        gal_sv_los_file = handle + '_den{}'.format(i) + '.CCF_gal_sv_los'

        fmt = 2*'%10.3f '
        
        np.savetxt(dm_delta_file, np.c_[rbin, dm_delta['den{}'.format(i)]], fmt=fmt)
        np.savetxt(gal_delta_file, np.c_[rbin, gal_delta['den{}'.format(i)]], fmt=fmt)
        np.savetxt(gal_vr_file, np.c_[rbin, gal_vr['den{}'.format(i)]], fmt=fmt)
        np.savetxt(gal_sv_los_file, np.c_[rbin, gal_sv_los['den{}'.format(i)]], fmt=fmt)



if __name__ == '__main__':
    split_densities()




