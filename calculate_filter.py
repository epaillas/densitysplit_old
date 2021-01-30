import sys
import click
import subprocess

@click.command()
@click.option('--data_filename', type=str, required=True)
@click.option('--centres_filename', type=str, required=True)
@click.option('--filter_filename', type=str, required=True)
@click.option('--boxsize', type=float, required=True)
@click.option('--ngrid', type=int, required=True)
@click.option('--filter_type', type=str, default='tophat')
@click.option('--filter_size', type=float, default=20)
@click.option('--dmin', type=float, default=0)
@click.option('--dmax', type=float)
@click.option('--qperp', type=float, default=1.0)
@click.option('--qpara', type=float, default=1.0)
def calculate_filter(data_filename,
                     centres_filename,
                     filter_filename,
                     boxsize,
                     ngrid,
                     dmin,
                     dmax,
                     get_filter,
                     filter_size,
                     filter_type,
                     qperp,
                     qpara):

        if dmax == None:
                if filter_type == 'tophat':
                        dmax = filter_size
                elif filter_type == 'gaussian':
                        dmax = 5 * filter_size

        logfile = filter_filename + '.log'
        binpath = sys.path[0] + '/bin/{}.exe'.format(filter_type)

        cmd = [binpath,
                data_filename,
                centres_filename,
                filter_filename,
                str(boxsize),
                str(dmin),
                str(dmax),
                str(filter_size),
                str(ngrid),
                str(qperp),
                str(qpara)]

        log = open(logfile, "w+")
        subprocess.call(cmd, stdout=log, stderr=log)

if __name__ == '__main__':
    calculate_filter()