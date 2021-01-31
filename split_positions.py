import numpy as np
from scipy.io import FortranFile
import click

# read file name from command line
@click.command()
@click.option('--centres_filename', type=str, required=True)
@click.option('--filter_filename', type=str, required=True)
@click.option('--quantiles', type=int, required=True)
def split_positions(centres_filename,
                    filter_filename,
                    quantiles):

    '''
    Splits a set of input positions in different quantiles,
    according to the local galaxy density.
    '''
    print('\nSplitting positions using the following arguments:')
    print('centres: {}'.format(centres_filename))
    print('filter_filename: {}'.format(filter_filename))
    print('quantiles: {}'.format(quantiles))

    # open centres file and get dimensions
    f = FortranFile(centres_filename, 'r')
    nrows = f.read_ints()[0]
    ncols = f.read_ints()[0]
    print('nrows, ncols= ({}, {})'.format(nrows, ncols))
    # read raw data and close file
    centres = f.read_reals(dtype=np.float64).reshape(nrows, ncols)
    f.close()

    # open filter file
    f = FortranFile(filter_filename, 'r')
    ncentres = f.read_ints()[0]
    print('ncentres: {}'.format(ncentres))
    smoothed_delta = f.read_reals(dtype=np.float64)
    idx = np.argsort(smoothed_delta)
    f.close()

    # sort profiles according to Delta(r=20mpc/h)
    sorted_centres = centres[idx]

    # divide profiles by their Delta(r=20mpc/h)
    binned_centres = {}
    for i in range(1, quantiles + 1):
        binned_centres['den{}'.format(i)] = sorted_centres[int((i-1)*ncentres/quantiles):int(i*ncentres/quantiles)]

        output_file = centres_filename + '.DS{}'.format(i)
        cout = binned_centres['den{}'.format(i)]
        print('Shape of cout: {}'.format(np.shape(cout)))
        f = FortranFile(output_file, 'w')
        f.write_record(np.shape(cout)[0])
        f.write_record(np.shape(cout)[1])
        f.write_record(cout)
        f.close()

        
if __name__ == '__main__':
    split_positions()




