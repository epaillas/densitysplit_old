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
                 get_monopole=True, get_rmu=False, get_spi=True,
                 get_filter=True, filter_type='tophat', 
                 filter_size=20, randoms_from_gal=False):

        # file names
        self.handle = handle
        self.tracer_file = tracer_file
        self.centres_file = centres_file

        self.get_monopole = get_monopole
        self.get_rmu = get_rmu
        self.get_spi = get_spi
        self.get_filter = get_filter
        self.is_matter = is_matter
        self.is_box = is_box
        self.nrandoms = int(nrandoms)
        self.box_size = box_size
        self.dmin = dmin
        self.dmax = dmax
        self.nrbins = nrbins
        self.nmubins = 80
        self.ngrid = ngrid
        self.filter_size = filter_size
        self.filter_type = filter_type
        self.randoms_from_gal = randoms_from_gal

        print('handle: {}'.format(self.handle))
        print('tracer_file: {}'.format(self.tracer_file))
        print('centres_file: {}'.format(self.centres_file))
        print('is_matter: {}'.format(self.is_matter))
        print('nrandoms: {}'.format(self.nrandoms))

        if os.path.isfile(self.centres_file):
            print('Centres file found. Skipping random centre generation.')
        else:
            self.GenerateRandomPoints()

        if self.get_monopole:
            print('Calculating CCF monopole.')
            self.CCF_monopole()
        if self.get_rmu:
            print('Calculating CCF in r-mu')
            self.CCF_rmu()
        if self.get_spi:
            print('Calculating CCF in sigma-pi')
            self.CCF_spi()
        if self.get_filter:
            if self.filter_type == 'gaussian':
                print('Calculating Gaussian smoothed Delta.')
                self.gaussian_filter()
            elif self.filter_type == 'tophat':
                print('Calculating top-hat Delta.')
                self.tophat_filter()
            else:
                sys.exit('Filter type not recognized. Aborting...')  



    def GenerateRandomPoints(self):
        '''
        Generates random points on a box
        of length and writes them down
        to an unformatted Fortran 90 file.
        '''
        np.random.seed(0)

        if self.randoms_from_gal:
            print('Randoms will be generated from galaxy positions.')
            fin = FortranFile(self.tracer_file, 'r')
            nrows = fin.read_ints()[0]
            ncols = fin.read_ints()[0]
            pos = fin.read_reals(dtype=np.float32).reshape(nrows, ncols)
            idx = np.random.choice(nrows, size=self.nrandoms, replace=False)
            cout = pos[idx]

        else:
            print('Randoms will be generated from a uniform distribution.')
            x = np.random.uniform(0, self.box_size, self.nrandoms)
            y = np.random.uniform(0, self.box_size, self.nrandoms)
            z = np.random.uniform(0, self.box_size, self.nrandoms)
            cout = np.c_[x, y, z]

        cout = cout.astype('float32')
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

    def gaussian_filter(self):
        '''
        Computes Gaussian smoothed Delta
        for a given filter size.
        '''
        if self.is_matter:
            fout = self.handle + '.DM_SmoothedDelta.unf'
            logfile = self.handle + '.DM_SmoothedDelta.log'
        else:
            fout = self.handle + '.gal_SmoothedDelta.unf'
            logfile = self.handle + '.gal_SmoothedDelta.log'

        if self.is_box:
            binpath = sys.path[0] + '/bin/'
            cmd = [binpath + 'gaussian_filter.exe',
                   self.tracer_file,
                   self.centres_file,
                   fout,
                   str(self.box_size),
                   str(self.dmin),
                   str(self.dmax),
                   str(self.filter_size),
                   str(self.ngrid)]
        
        log = open(logfile, "w+")
        subprocess.call(cmd, stdout=log, stderr=log)

    def tophat_filter(self):
        '''
        Computes top-hat Delta
        for a given filter size.
        '''
        if self.is_matter:
            fout = self.handle + '.DM_TopHatDelta.unf'
            logfile = self.handle + '.DM_TopHatDelta.log'
        else:
            fout = self.handle + '.gal_TopHatDelta.unf'
            logfile = self.handle + '.gal_TopHatDelta.log'

        if self.is_box:
            binpath = sys.path[0] + '/bin/'
            cmd = [binpath + 'tophat_filter.exe',
                   self.tracer_file,
                   self.centres_file,
                   fout,
                   str(self.box_size),
                   str(self.dmin),
                   str(self.dmax),
                   str(self.filter_size),
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
                   str(self.nmubins),
                   str(self.ngrid)]
        
        log = open(logfile, "w+")
        subprocess.call(cmd, stdout=log, stderr=log)

    def CCF_spi(self):
        '''
        Computes delta(r) profiles from
        the random centres in bins of s and pi.
        '''
        if self.is_matter:
            fout = self.handle + '.CCF_DM_spi.unf'
            logfile = self.handle + '.CCF_DM_spi.log'
        else:
            fout = self.handle + '.CCF_gal_spi.unf'
            logfile = self.handle + '.CCF_gal_spi.log'

        if self.is_box:
            binpath = sys.path[0] + '/bin/'
            cmd = [binpath + 'CCF_spi.exe',
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

