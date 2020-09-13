from python_tools.densitysplitter import DensitySplitter
import click


@click.command()
@click.option('--handle', type=str, required=True)
@click.option('--tracer_file', type=str, required=True)
@click.option('--centres_file', type=str, required=True)
@click.option('--nrandoms', type=int, required=True)
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
@click.option('--randoms_from_tracers', type=bool, default=False)
@click.option('--qperp', type=float, default=1.0)
@click.option('--qpara', type=float, default=1.0)
def run_density_splitter(handle,
                         tracer_file,
                         centres_file,
                         nrandoms,
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
                         randoms_from_tracers,
                         qperp,
                         qpara):

    ds = DensitySplitter(handle=handle,
                         tracer_file=tracer_file,
                         centres_file=centres_file,
                         nrandoms=nrandoms,
                         boxsize_x=boxsize_x,
                         boxsize_y=boxsize_y,
                         boxsize_z=boxsize_z,
                         ngrid_x=ngrid_x,
                         ngrid_y=ngrid_y,
                         ngrid_z=ngrid_z,
                         dmin=dmin,
                         dmax=dmax,
                         get_filter=get_filter,
                         filter_size=filter_size,
                         filter_type=filter_type,
                         randoms_from_tracers=randoms_from_tracers,
                         qperp=qperp,
                         qpara=qpara)


if __name__ == '__main__':
    run_density_splitter()
