import numpy as np
import sys
import os
import glob
import subprocess
from astropy.io import fits
from scipy.io import FortranFile

class DensitySplitter:

    def __init__(self, handle, tracer_file, centres_file, nrandoms, box_size,
                 dmin, dmax, nrbins, is_box=True, ngrid=100, is_matter=False):

        steps = [int(i) for i in steps.split(',')]

        # file names
        self.handle = handle
        self.tracer_file = tracer_file
        self.centres_file = centres_file

        self.is_matter = is_matter
        self.is_box = is_box
        self.nrandoms = int(nrandoms)
        self.box_size = box_size
        self.dmin = dmin
        self.dmax = dmax
        self.nrbins = nrbins
        self.ngrid = ngrid
        self.steps = steps

        if os.path.isfile(self.centres_file):
            print('Centres file found. Skipping random centre generation.')
        else:
            self.GenerateRandomPoints()

        print('Calculating density/velocity profiles.')
        self.DensityProfiles()

    def GenerateRandomPoints(self):
        '''
        Generates random points on a box
        of length and writes them down
        to an unformatted Fortran 90 file.
        '''

        x = np.random.uniform(0, self.box_size, self.nrandoms)
        y = np.random.uniform(0, self.box_size, self.nrandoms)
        z = np.random.uniform(0, self.box_size, self.nrandoms)

        cout = np.c_[x, y, z]

        f = FortranFile(self.centres_file, 'w')
        nrows, ncols = np.shape(cout)
        f.write_record(nrows)
        f.write_record(ncols)
        f.write_record(cout)
        f.close()

        return

    def DensityProfiles(self):
        '''
        Computes delta(r), Delta(r),
        vel(r), DD(R) profiles from
        the random centres.
        '''
        if self.is_matter:
            fout = self.handle + '.DM_den.unf'
        else:
            fout = self.handle + '.gal_den.unf'

        if self.is_box:
            binpath = sys.path[0] + '/bin/'
            cmd = [binpath + 'density_profiles.exe',
                   self.tracer_file,
                   self.centres_file,
                   fout,
                   str(self.box_size),
                   str(self.dmin),
                   str(self.dmax),
                   str(self.nrbins),
                   str(self.ngrid)]
        
        logfile = self.handle + '.den.log'
        log = open(logfile, "w+")
        subprocess.call(cmd, stdout=log, stderr=log)

