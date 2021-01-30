import sys
import click

@click.command()
@click.option('--data_filename', type=str, required=True)
@click.option('--centres_filename', type=str, required=True)
@click.option('--filter_filename', type=str, required=True)
@click.option('--boxsize_x', type=float, required=True)
@click.option('--boxsize_y', type=float, required=True)
@click.option('--boxsize_z', type=float, required=True)
@click.option('--ngrid_x', type=int, required=True)
@click.option('--ngrid_y', type=int, required=True)
@click.option('--ngrid_z', type=int, required=True)
@click.option('--dmin', type=float, required=True)
@click.option('--dmax', type=float, required=True)
@click.option('--get_filter', type=bool, default=True)
@click.option('--filter_size', type=float, default=20)
@click.option('--filter_type', type=str, default='tophat')
@click.option('--qperp', type=float, default=1.0)
@click.option('--qpara', type=float, default=1.0)

def tophat_filter(data_filename,
                  centres_filename,
                  filter_filename,
                  boxsize_x,
                  boxsize_y,
                  boxsize_z,
                  ngrid_x,
                  ngrid_y,
                  ngrid_z,
                  dmin,
                  dmax,
                  get_filter,
                  filter_size,
                  filter_type,
                  qperp,
                  qpara)  

logfile = filter_filename + '.log'
binpath = sys.path[0] + '/bin/{}.exe'.format(filter_type)

cmd = [binpath,
        data_filename,
        centres_filename,
        filter_filename,
        str(boxsize_x),
        str(boxsize_y),
        str(boxsize_z),
        str(dmin),
        str(dmax),
        str(filter_size),
        str(ngrid_x),
        str(ngrid_y),
        str(ngrid_z),
        str(qperp),
        str(qpara)]

log = open(logfile, "w+")
subprocess.call(cmd, stdout=log, stderr=log)
