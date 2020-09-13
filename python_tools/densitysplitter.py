import numpy as np
import sys
import os
import glob
import subprocess
from astropy.io import fits
from scipy.io import FortranFile


class DensitySplitter:

    def __init__(self,
                 handle,
                 tracer_file,
                 centres_file,
                 nrandoms,
                 boxsize_x,
                 boxsize_y,
                 boxsize_z,
                 dmin,
                 dmax,
                 ngrid_x,
                 ngrid_y,
                 ngrid_z,
                 get_filter=True,
                 filter_type='tophat',
                 filter_size=20,
                 randoms_from_tracers=False,
                 qperp=1.0,
                 qpara=1.0):

        # file names
        self.handle = handle
        self.tracer_file = tracer_file
        self.centres_file = centres_file

        self.get_filter = get_filter
        self.nrandoms = int(nrandoms)
        self.boxsize_x = boxsize_x
        self.boxsize_y = boxsize_y
        self.boxsize_z = boxsize_z
        self.dmin = dmin
        self.dmax = dmax
        self.ngrid_x = ngrid_x
        self.ngrid_y = ngrid_y
        self.ngrid_z = ngrid_z
        self.filter_size = filter_size
        self.filter_type = filter_type
        self.randoms_from_tracers = randoms_from_tracers
        self.qperp = qperp
        self.qpara = qpara

        print('handle: {}'.format(self.handle))
        print('tracer_file: {}'.format(self.tracer_file))
        print('centres_file: {}'.format(self.centres_file))
        print('nrandoms: {}'.format(self.nrandoms))

        if os.path.isfile(self.centres_file):
            print('Centres file found. Skipping random centre generation.')
        else:
            self.random_points()

        if self.get_filter:
            getattr(self, filter_type + '_filter')()

    def random_points(self):
        '''
        Generates random points on a box
        of length and writes them down
        to an unformatted Fortran 90 file.
        '''
        np.random.seed(0)

        if self.randoms_from_tracers:
            print('Randoms will be generated from galaxy positions.')
            fin = FortranFile(self.tracer_file, 'r')
            nrows = fin.read_ints()[0]
            ncols = fin.read_ints()[0]
            pos = fin.read_reals(dtype=np.float64).reshape(nrows, ncols)
            idx = np.random.choice(nrows, size=self.nrandoms, replace=False)
            cout = pos[idx]

        else:
            print('Randoms will be generated from a uniform distribution.')
            x = np.random.uniform(0, self.boxsize_x, self.nrandoms)
            y = np.random.uniform(0, self.boxsize_y, self.nrandoms)
            z = np.random.uniform(0, self.boxsize_z, self.nrandoms)
            cout = np.c_[x, y, z]

        cout = cout.astype('float64')
        f = FortranFile(self.centres_file, 'w')
        nrows, ncols = np.shape(cout)
        f.write_record(nrows)
        f.write_record(ncols)
        f.write_record(cout)
        f.close()

        return

    def gaussian_filter(self):
        '''
        Computes Gaussian smoothed Delta
        for a given filter size.
        '''
        print('Calculating filter using Gaussian smoothing.')
        fout = self.handle + '.gal_GaussianDelta.unf'
        logfile = self.handle + '.gal_GaussianDelta.log'

        binpath = sys.path[0] + '/bin/'
        cmd = [binpath + 'gaussian_filter.exe',
               self.tracer_file,
               self.centres_file,
               fout,
               str(self.boxsize_x),
               str(self.boxsize_y),
               str(self.boxsize_z),
               str(self.dmin),
               str(self.dmax),
               str(self.filter_size),
               str(self.ngrid_x),
               str(self.ngrid_y),
               str(self.ngrid_z),
               str(self.qperp),
               str(self.qpara)]

        log = open(logfile, "w+")
        subprocess.call(cmd, stdout=log, stderr=log)

    def tophat_filter(self):
        '''
        Computes top-hat Delta
        for a given filter size.
        '''
        print('Calculating filter using Top-hat smoothing.')
        fout = self.handle + '.gal_TopHatDelta.unf'
        logfile = self.handle + '.gal_TopHatDelta.log'

        binpath = sys.path[0] + '/bin/'
        cmd = [binpath + 'tophat_filter.exe',
               self.tracer_file,
               self.centres_file,
               fout,
               str(self.boxsize_x),
               str(self.boxsize_y),
               str(self.boxsize_z),
               str(self.dmin),
               str(self.dmax),
               str(self.filter_size),
               str(self.ngrid_x),
               str(self.ngrid_y),
               str(self.ngrid_z),
               str(self.qperp),
               str(self.qpara)]

        log = open(logfile, "w+")
        subprocess.call(cmd, stdout=log, stderr=log)
