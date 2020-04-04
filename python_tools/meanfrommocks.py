import numpy as np
import glob
import click
import sys

@click.command()

@click.option('--handle_in', type=str, required=True)
@click.option('--handle_out', type=str, required=True)
@click.option('--is_velocity', type=bool, required=True)
@click.option('--is_matter', type=bool, required=True)
@click.option('--corr_type', type=str, default='monopole')

def mean_from_mocks(handle_in,
                    handle_out,
                    is_velocity,
                    corr_type):

    print('\nAveraging mean from mocks for the following arguments:')
    print('handle_in: {}'.format(handle_in))
    print('handle_out: {}'.format(handle_out))
    print('is_velocity: {}'.format(is_velocity))

    # possible file extensions
    if corr_type == 'monopole':
        file_ext = ['CCF_gal_monopole', 'CCF_DM_monopole']
    elif corr_type =='rmu':
        file_ext = ['CCF_gal_rmu']
    else:
        sys.exit('Correlation type not recognized.')

    if is_velocity:
        file_ext = ['CCF_gal_vr', 'CCF_gal_svlos']

    # loop over all mocks and calculate mean
    for ext in file_ext:
        print('\nAveraging extension: {}'.format(ext))
        hin = handle_in + '.{}'.format(ext)
        hout = handle_out + '.{}'.format(ext)

        mock_files = sorted(glob.glob(hin))
        data_list = []

        for mock_file in mock_files:
            data = np.genfromtxt(mock_file)
            data[np.isnan(data)] = -1
            data[data == np.inf] = -1
            data_list.append(data)

        data_list = np.asarray(data_list)
        cout = np.mean(data_list, axis=0)

        print('np.shape(data_list): {}'.format(np.shape(data_list)))
        print('np.shape(cout): {}'.format(np.shape(cout)))

        fmt = np.shape(cout)[1] * '%10.3f '

        np.savetxt(hout, cout, fmt=fmt)

if __name__ == '__main__':
    mean_from_mocks()