import numpy as np
import sys
import os
import glob
import subprocess
from astropy.io import fits
from scipy.io import FortranFile

class DensitySplitter:

    def __init__(self, handle, tracer_file, centres_file, nrandoms, box_size,
                 dmin, dmax, nrbins, is_box=True, ngrid=100, is_matter=False,
                 get_monopole=True, get_rmu=False, randoms_from_gal=False):

        # file names
        self.handle = handle
        self.tracer_file = tracer_file
        self.centres_file = centres_file

        self.get_monopole = get_monopole
        self.get_rmu = get_rmu
        self.is_matter = is_matter
        self.is_box = is_box
        self.nrandoms = int(nrandoms)
        self.box_size = box_size
        self.dmin = dmin
        self.dmax = dmax
        self.nrbins = nrbins
        self.ngrid = ngrid

        print('handle: {}'.format(self.handle))
        print('tracer_file: {}'.format(self.tracer_file))
        print('centres_file: {}'.format(self.centres_file))
        print('is_matter: {}'.format(self.is_matter))
        print('nrandoms: {}'.format(self.nrandoms))

        if os.path.isfile(self.centres_file):
            print('Centres file found. Skipping random centre generation.')
        else:
            self.GenerateRandomPoints()

        print('Calculating density/velocity profiles.')
        if self.get_monopole:
            self.CCF_monopole()
        if self.get_rmu:
            self.CCF_rmu()


    def GenerateRandomPoints(self, randoms_from_gal=False):
        '''
        Generates random points on a box
        of length and writes them down
        to an unformatted Fortran 90 file.
        '''

        if randoms_from_gal:
            fin = FortranFile(self.tracer_file, 'r')
            nrows = fin.read_ints()[0]
            ncols = fin.read_ints()[0]
            pos = fin.read_reals().reshape(nrows, ncols)
            cout = np.random.choice(pos, size=self.nrandoms, replace=False)

        else:
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

    def CCF_monopole(self):
        '''
        Computes delta(r), Delta(r),
        vel(r), DD(R) profiles from
        the random centres.
        '''
        if self.is_matter:
            fout = self.handle + '.CCF_DM_monopole.unf'
            logfile = self.handle + '.CCF_DM_monopole.log'
        else:
            fout = self.handle + '.CCF_gal_monopole.unf'
            logfile = self.handle + '.CCF_gal_monopole.log'

        if self.is_box:
            binpath = sys.path[0] + '/bin/'
            cmd = [binpath + 'CCF_monopole.exe',
                   self.tracer_file,
                   self.centres_file,
                   fout,
                   str(self.box_size),
                   str(self.dmin),
                   str(self.dmax),
                   str(self.nrbins),
                   str(self.ngrid)]
        
        log = open(logfile, "w+")
        subprocess.call(cmd, stdout=log, stderr=log)

    def CCF_rmu(self):
        '''
        Computes delta(r), Delta(r),
        DD(R) profiles from
        the random centres.
        '''
        if self.is_matter:
            fout = self.handle + '.CCF_DM_rmu.unf'
            logfile = self.handle + '.CCF_DM_rmu.log'
        else:
            fout = self.handle + '.CCF_gal_rmu.unf'
            logfile = self.handle + '.CCF_gal_rmu.log'

        if self.is_box:
            binpath = sys.path[0] + '/bin/'
            cmd = [binpath + 'CCF_rmu.exe',
                   self.tracer_file,
                   self.centres_file,
                   fout,
                   str(self.box_size),
                   str(self.dmin),
                   str(self.dmax),
                   str(self.nrbins),
                   str(self.ngrid)]
        
        log = open(logfile, "w+")
        subprocess.call(cmd, stdout=log, stderr=log)

