import numpy as np
from jointfit import Model1
import os
import sys
import argparse
import emcee
from multiprocessing import Pool

def log_probability(theta):
        lp = model.log_prior(theta)
        if not np.isfinite(lp):
            return -np.inf
        return lp + model.log_likelihood(theta)

parser = argparse.ArgumentParser(description='MCMC for the void-galaxy correlation function.')

parser.add_argument('--ncores', type=int)
parser.add_argument('--xi_smu', type=str)
parser.add_argument('--xi_r', type=str)
parser.add_argument('--delta_r', type=str)
parser.add_argument('--sv_r', type=str)
parser.add_argument('--covmat', type=str)
parser.add_argument('--full_fit', type=int)
parser.add_argument('--smin', type=str)
parser.add_argument('--smax', type=str)
parser.add_argument('--model', type=int)

args = parser.parse_args()  

os.environ["OMP_NUM_THREADS"] = "1"

if args.model == 1:
    model = Model1(delta_r_files=args.delta_r, xi_r_files=args.xi_r, sv_files=args.sv_r,
                    xi_smu_files=args.xi_smu, covmat_file=args.covmat,
                    full_fit=args.full_fit, smins=args.smin, smaxs=args.smax)

    if args.full_fit == 1:
        backend_name = 'Model1_Joint_FullFit_{}-{}.h5'.format(args.smin, args.smax)
    else:
        backend_name = 'Model1_Joint_QuadFit_{}-{}.h5'.format(args.smin, args.smax)
    ndim = 4
    nwalkers = 32
    niter = 5000

    fs8 = 0.472
    sigma_v1 = 300
    sigma_v2 = 300
    epsilon = 1.0

    start_params = np.array([fs8, sigma_v1, sigma_v2, epsilon])
    scales = [1, 1000, 1000, 1]

    p0 = [start_params + 1e-2 * np.random.randn(ndim) * scales for i in range(nwalkers)]

    print('Running emcee with the following parameters:')
    print('nwalkers: ' + str(nwalkers))
    print('ndim: ' + str(ndim))
    print('niter: ' + str(niter))
    print('backend: ' + backend_name)
    print('Running in {} CPUs'.format(args.ncores))

    backend = emcee.backends.HDFBackend(backend_name)
    backend.reset(nwalkers, ndim)

    with Pool(processes=args.ncores) as pool:

        sampler = emcee.EnsembleSampler(nwalkers, ndim,
                                        log_probability,
                                        backend=backend,
                                        pool=pool)
        sampler.run_mcmc(p0, niter, progress=True)


if args.model == 2:
    model = Model1(delta_r_files=args.delta_r, xi_r_files=args.xi_r,
                    xi_smu_files=args.xi_smu, covmat_file=args.covmat,
                    full_fit=args.full_fit, smins=args.smin, smaxs=args.smax)

    if args.full_fit == 1:
        backend_name = 'Model2_Joint_FullFit_{}-{}.h5'.format(args.smin, args.smax)
    else:
        backend_name = 'Model2_Joint_QuadFit_{}-{}.h5'.format(args.smin, args.smax)
    ndim = 2
    nwalkers = 32
    niter = 5000

    fs8 = 0.472
    epsilon = 1.0

    start_params = np.array([fs8, epsilon])
    scales = [1, 1]

    p0 = [start_params + 1e-2 * np.random.randn(ndim) * scales for i in range(nwalkers)]

    print('Running emcee with the following parameters:')
    print('nwalkers: ' + str(nwalkers))
    print('ndim: ' + str(ndim))
    print('niter: ' + str(niter))
    print('backend: ' + backend_name)
    print('Running in {} CPUs'.format(args.ncores))

    backend = emcee.backends.HDFBackend(backend_name)
    backend.reset(nwalkers, ndim)

    with Pool(processes=args.ncores) as pool:

        sampler = emcee.EnsembleSampler(nwalkers, ndim,
                                        log_probability,
                                        backend=backend,
                                        pool=pool)
        sampler.run_mcmc(p0, niter, progress=True)


