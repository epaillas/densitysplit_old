import numpy as np
from scipy.io import FortranFile
import click
import matplotlib.pyplot as plt


# read file name from command line
@click.command()
@click.option('--gal_den_monopole', type=str, required=True)
@click.option('--dm_den_monopole', type=str, required=True)
@click.option('--filter_file', type=str, required=True)
@click.option('--has_velocity', type=bool, required=True)
@click.option('--handle', type=str, required=True)
@click.option('--ndenbins', type=int, required=True)

def split_densities(gal_den_monopole,
                    dm_den_monopole,
                    filter_file,
                    has_velocity,
                    handle,
                    ndenbins):


    print('\nSplitting densities for the following arguments:')
    print('gal_den_monopole: {}'.format(gal_den_monopole))
    print('dm_den_monopole: {}'.format(dm_den_monopole))
    print('filter_file: {}'.format(filter_file))
    print('has_velocity: {}'.format(has_velocity))
    print('handle: {}'.format(handle))

    # open galaxy density file and get dimensions
    f = FortranFile(gal_den_monopole, 'r')
    ncentres = f.read_ints()[0]
    nbins = f.read_ints()[0]
    print('ncentres, nbins = ({}, {})'.format(ncentres, nbins))

    # read raw data and close file
    rbin = f.read_reals(dtype=np.float64)
    gal_dd = f.read_reals(dtype=np.float64).reshape(nbins, ncentres)
    xi_r = f.read_reals(dtype=np.float64).reshape(nbins, ncentres)
    xibar_r = f.read_reals(dtype=np.float64).reshape(nbins, ncentres)


    if has_velocity:
        gal_vv_r = f.read_reals(dtype=np.float64).reshape(nbins, ncentres)
        gal_vv_los = f.read_reals(dtype=np.float64).reshape(nbins, ncentres)
        gal_vv2_los = f.read_reals(dtype=np.float64).reshape(nbins, ncentres)

    f.close()

    # open DM density file
    f = FortranFile(dm_den_monopole, 'r')
    ncentres = f.read_ints()[0]
    nbins = f.read_ints()[0]
    print('ncentres, nbins = ({}, {})'.format(ncentres, nbins))

    rbin = f.read_reals(dtype=np.float64)
    dm_dd = f.read_reals(dtype=np.float64).reshape(nbins, ncentres)
    delta_r = f.read_reals(dtype=np.float64).reshape(nbins, ncentres)
    Delta_r = f.read_reals(dtype=np.float64).reshape(nbins, ncentres)

    if has_velocity:
        dm_vv_r = f.read_reals(dtype=np.float64).reshape(nbins, ncentres)
        dm_vv_los = f.read_reals(dtype=np.float64).reshape(nbins, ncentres)
        dm_vv2_los = f.read_reals(dtype=np.float64).reshape(nbins, ncentres)

    f.close()

    # open filter file
    f = FortranFile(filter_file, 'r')
    ncentres = f.read_ints()[0]
    print('ncentres: {}'.format(ncentres))
    smoothed_delta = f.read_reals(dtype=np.float64)
    idx = np.argsort(smoothed_delta)


    # assemble individual profiles into an array
    if has_velocity:
        profiles = [np.c_[gal_dd[:,i],
                          xi_r[:,i],
                          xibar_r[:,i],
                          delta_r[:,i],
                          Delta_r[:,i],
                          gal_vv_r[:,i],
                          gal_vv_los[:,i],
                          gal_vv2_los[:,i]] for i in range(ncentres)]
    else:
        profiles = [np.c_[gal_dd[:,i],
                          xi_r[:,i],
                          xibar_r[:,i],
                          Delta_r[:,i]] for i in range(ncentres)]

    profiles = np.asarray(profiles)

    # sort profiles according to Delta(r=20mpc/h)
    sorted_profiles = profiles[idx]

    # divide profiles by their Delta(r=20mpc/h)
    profiles = {}
    for i in range(1, ndenbins + 1):
        profiles['den{}'.format(i)] = sorted_profiles[int((i-1)*ncentres/ndenbins):int(i*ncentres/ndenbins)]

    # get average densities
    xi_r = {}
    delta_r = {}
    xibar_r = {}
    Delta_r = {}
    dd = {}
    gal_vr = {}
    gal_sv_los = {}


    for i in range(1, ndenbins + 1):
        den = profiles['den{}'.format(i)]
        dd = np.sum(den, axis=0)[:,-8]
        xi_r['den{}'.format(i)] = np.mean(den, axis=0)[:,-7]
        xibar_r['den{}'.format(i)] = np.mean(den, axis=0)[:,-6]
        delta_r['den{}'.format(i)] = np.mean(den, axis=0)[:,-5]
        Delta_r['den{}'.format(i)] = np.mean(den, axis=0)[:,-4]
        vv_r = np.sum(den, axis=0)[:,-3]
        vv_los = np.sum(den, axis=0)[:,-2]
        vv2_los = np.sum(den, axis=0)[:,-1]

        gal_vr['den{}'.format(i)] = vv_r / dd
        gal_sv_los['den{}'.format(i)] = np.sqrt((vv2_los - (vv_los**2 / dd)) / (dd - 1))


    for i in range(1, ndenbins + 1):
        delta_r_file = handle + '_den{}'.format(i) + '.CCF_DM_monopole'
        xi_r_file = handle + '_den{}'.format(i) + '.CCF_gal_monopole'
        gal_vr_file = handle + '_den{}'.format(i) + '.CCF_gal_vr'
        gal_sv_los_file = handle + '_den{}'.format(i) + '.CCF_gal_svlos'

        fmt = 2*'%10.3f '
        
        np.savetxt(delta_r_file, np.c_[rbin, delta_r['den{}'.format(i)]], fmt=fmt)
        np.savetxt(xi_r_file, np.c_[rbin, xi_r['den{}'.format(i)]], fmt=fmt)
        np.savetxt(gal_vr_file, np.c_[rbin, gal_vr['den{}'.format(i)]], fmt=fmt)
        np.savetxt(gal_sv_los_file, np.c_[rbin, gal_sv_los['den{}'.format(i)]], fmt=fmt)



if __name__ == '__main__':
    split_densities()




