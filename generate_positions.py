import numpy as np
import sys
from scipy.io import FortranFile
import click

@click.command()
@click.option('--data_filename', type=str, required=True)
@click.option('--centres_filename', type=str, required=True)
@click.option('--nrandoms', type=int, required=True)
@click.option('--sampling', type=str, required=True)
@click.option('--boxsize', type=float, required=True)
def generate_positions(data_filename,
                     centres_filename,
                     nrandoms,
                     sampling,
                     boxsize):
        '''
        Generates random points on a box
        writes them to an unformatted
        Fortran 90 file.
        '''
        np.random.seed(0)

        if sampling == 'tracers':
            print('Randoms will be generated from galaxy positions.')
            fin = FortranFile(data_filename, 'r')
            nrows = fin.read_ints()[0]
            ncols = fin.read_ints()[0]
            pos = fin.read_reals(dtype=np.float64).reshape(nrows, ncols)
            idx = np.random.choice(nrows, size=nrandoms, replace=False)
            cout = pos[idx]

        elif sampling == 'uniform':
            print('Randoms will be generated from a uniform distribution.')
            x = np.random.uniform(0, boxsize, nrandoms)
            y = np.random.uniform(0, boxsize, nrandoms)
            z = np.random.uniform(0, boxsize, nrandoms)
            cout = np.c_[x, y, z]
        else:
            sys.exit('Sampling type not recognized')

        cout = cout.astype('float64')
        f = FortranFile(centres_filename, 'w')
        nrows, ncols = np.shape(cout)
        f.write_record(nrows)
        f.write_record(ncols)
        f.write_record(cout)
        f.close()

if __name__ == '__main__':
    generate_positions()
