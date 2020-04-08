import numpy as np
import matplotlib.pyplot as plt
from rsdmodel import NadathurAndPercival
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
parser.add_argument('--smin', type=float)
parser.add_argument('--smax', type=float)

args = parser.parse_args()  

os.environ["OMP_NUM_THREADS"] = "1"


model = NadathurAndPercival(delta_r_file=args.delta_r, xi_r_file=args.xi_r, sv_file=args.sv_r,
                 xi_smu_file=args.xi_smu, covmat_file=args.covmat,
                 full_fit=args.full_fit, smin=args.smin, smax=args.smax)

if args.full_fit == 1:
    backend_name = args.xi_smu + '_RSD_FullFit_{}-{}.h5'.format(args.smin, args.smax)
else:
    backend_name = args.xi_smu + '_RSD_QuadFit_{}-{}.h5'.format(args.smin, args.smax)
ndim = 3
nwalkers = 32
niter = 5000

fs8 = 0.472
sigma_v = 300
epsilon = 1.0

start_params = np.array([fs8, sigma_v, epsilon])
scales = [1, 1000, 1]

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
