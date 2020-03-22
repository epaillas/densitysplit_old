import numpy as np
import glob
import click

@click.command()

@click.option('--handle_in', type=str, required=True)
@click.option('--handle_out', type=str, required=True)
@click.option('--has_velocity', type=bool, required=True)

def mean_from_mocks(handle_in,
                    handle_out,
                    has_velocity):

    print('\nAveraging mean from mocks for the following arguments:')
    print('handle_in: {}'.format(handle_in))
    print('handle_out: {}'.format(handle_out))
    print('has_velocity: {}'.format(has_velocity))

    # possible file extensions
    corr_type = ['CCF_gal_monopole', 'CCF_DM_monopole']

    if has_velocity:
        corr_type = corr_type + ['CCF_gal_monopole', 'CCF_DM_monopole', 'CCF_gal_vr', 'CCF_gal_svlos']

    # loop over all mocks and calculate mean
    for ext in corr_type:
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

        np.savetxt(hout, cout)

if __name__ == '__main__':
    mean_from_mocks()