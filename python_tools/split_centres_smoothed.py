import numpy as np
from scipy.io import FortranFile
import click

# read file name from command line
@click.command()
@click.option('--centres_file', type=str, required=True)
@click.option('--filter_file', type=str, required=True)
@click.option('--handle', type=str, required=True)
@click.option('--ndenbins', type=int, required=True)

def split_densities(centres_file,
                    filter_file,
                    handle,
                    ndenbins):


    print('\nSplitting centres for the following arguments:')
    print('centres: {}'.format(centres_file))
    print('filter_file: {}'.format(filter_file))
    print('handle: {}'.format(handle))

    # open centres file and get dimensions
    f = FortranFile(centres_file, 'r')
    nrows = f.read_ints()[0]
    ncols = f.read_ints()[0]
    print('nrows, ncols= ({}, {})'.format(nrows, ncols))
    # read raw data and close file
    centres = f.read_reals(dtype=np.float32).reshape(nrows, ncols)
    f.close()

    # open filter file
    f = FortranFile(filter_file, 'r')
    ncentres = f.read_ints()[0]
    print('ncentres: {}'.format(ncentres))
    smoothed_delta = f.read_reals(dtype=np.float32)
    idx = np.argsort(smoothed_delta)
    f.close()

    # sort profiles according to Delta(r=20mpc/h)
    sorted_centres = centres[idx]

    # divide profiles by their Delta(r=20mpc/h)
    binned_centres = {}
    for i in range(1, ndenbins + 1):
        binned_centres['den{}'.format(i)] = sorted_centres[int((i-1)*ncentres/ndenbins):int(i*ncentres/ndenbins)]

        out_file = handle + '_den{}'.format(i) + '.cen.unf'
        cout = binned_centres['den{}'.format(i)] 
        f = FortranFile(out_file, 'w')
        f.write_record(np.shape(cout)[0])
        f.write_record(np.shape(cout)[1])
        f.write_record(cout)
        f.close()

        
if __name__ == '__main__':
    split_densities()




