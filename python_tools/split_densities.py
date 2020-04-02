import numpy as np
from scipy.io import FortranFile
import click
import matplotlib.pyplot as plt


# read file name from command line
@click.command()
@click.option('--gal_den_monopole', type=str, required=True)
@click.option('--dm_den_monopole', type=str, required=True)
@click.option('--has_velocity', type=bool, required=True)
@click.option('--handle', type=str, required=True)

def split_densities(gal_den_monopole,
                    dm_den_monopole,
                    has_velocity,
                    handle):


    print('\nSplitting densities for the following arguments:')
    print('gal_den_monopole: {}'.format(gal_den_monopole))
    print('dm_den_monopole: {}'.format(dm_den_monopole))
    print('has_velocity: {}'.format(has_velocity))
    print('handle: {}'.format(handle))

    # open galaxy density file and get dimensions
    f = FortranFile(gal_den_monopole, 'r')
    ncentres = f.read_ints()[0]
    nbins = f.read_ints()[0]
    print('ncentres, nbins = ({}, {})'.format(ncentres, nbins))

    # read raw data and close file
    rbin = f.read_reals()
    gal_dd = f.read_reals(dtype=np.float64).reshape(nbins, ncentres)
    xi_r = f.read_reals(dtype=np.float64).reshape(nbins, ncentres)
    xibar_r = f.read_reals(dtype=np.float64).reshape(nbins, ncentres)


    if has_velocity:
        gal_vr = f.read_reals(dtype=np.float64).reshape(nbins, ncentres)
        gal_sv_los = f.read_reals(dtype=np.float64).reshape(nbins, ncentres)

    f.close()

    # open DM density file
    f = FortranFile(dm_den_monopole, 'r')
    ncentres = f.read_ints()[0]
    nbins = f.read_ints()[0]
    print('ncentres, nbins = ({}, {})'.format(ncentres, nbins))

    rbin = f.read_reals()
    dm_dd = f.read_reals(dtype=np.float64).reshape(nbins, ncentres)
    delta_r = f.read_reals(dtype=np.float64).reshape(nbins, ncentres)
    Delta_r = f.read_reals(dtype=np.float64).reshape(nbins, ncentres)

    if has_velocity:
        dm_vr = f.read_reals(dtype=np.float64).reshape(nbins, ncentres)
        dm_sv_los = f.read_reals(dtype=np.float64).reshape(nbins, ncentres)

    f.close()


    # assemble individual profiles into an array
    if has_velocity:
        profiles = [np.c_[gal_dd[:,i], xi_r[:,i], xibar_r[:,i],
                    delta_r[:,i], Delta_r[:,i],
                    gal_vr[:,i], gal_sv_los[:,i]] for i in range(ncentres)]
    else:
        profiles = [np.c_[gal_dd[:,i], xi_r[:,i], xibar_r[:,i],
                    Delta_r[:,i]]
                    for i in range(ncentres)]

    profiles = np.asarray(profiles)

    # sort profiles according to Delta(r=20mpc/h)
    sorted_profiles = profiles[np.argsort(profiles[:, 10, -5])]

    # divide profiles by their Delta(r=20mpc/h)
    profiles = {}
    for i in range(1, 6):
        profiles['den{}'.format(i)] = sorted_profiles[int((i-1)*ncentres/5):int(i*ncentres/5)]

    # get average densities
    xi_r = {}
    delta_r = {}
    xibar_r = {}
    Delta_r = {}

    for i in range(1, 6):
        den = profiles['den{}'.format(i)]
        xi_r['den{}'.format(i)] = np.mean(den, axis=0)[:,-6]
        xibar_r['den{}'.format(i)] = np.mean(den, axis=0)[:,-5]
        delta_r['den{}'.format(i)] = np.mean(den, axis=0)[:,-4]
        Delta_r['den{}'.format(i)] = np.mean(den, axis=0)[:,-3]

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
        delta_r_file = handle + '_den{}'.format(i) + '.CCF_dm_den_monopole'
        xi_r_file = handle + '_den{}'.format(i) + '.CCF_gal_den_monopole'
        gal_vr_file = handle + '_den{}'.format(i) + '.CCF_gal_vr'
        gal_sv_los_file = handle + '_den{}'.format(i) + '.CCF_gal_sv_los'

        fmt = 2*'%10.3f '
        
        np.savetxt(delta_r_file, np.c_[rbin, delta_r['den{}'.format(i)]], fmt=fmt)
        np.savetxt(xi_r_file, np.c_[rbin, xi_r['den{}'.format(i)]], fmt=fmt)
        np.savetxt(gal_vr_file, np.c_[rbin, gal_vr['den{}'.format(i)]], fmt=fmt)
        np.savetxt(gal_sv_los_file, np.c_[rbin, gal_sv_los['den{}'.format(i)]], fmt=fmt)



if __name__ == '__main__':
    split_densities()




