from python_tools.densitysplitter import DensitySplitter
import click

@click.command()
@click.option('--handle', type=str, required=True)
@click.option('--tracer_file', type=str, required=True)
@click.option('--centres_file', type=str, required=True)
@click.option('--nrandoms', type=int, required=True)
@click.option('--box_size', type=float, required=True)
@click.option('--is_matter', type=bool, default=False)
@click.option('--dmin', type=float, default=0)
@click.option('--dmax', type=float, default=100)
@click.option('--nrbins', type=int, default=50)
@click.option('--get_monopole', type=bool, default=True)
@click.option('--get_rmu', type=bool, default=False)
@click.option('--get_spi', type=bool, default=False)
@click.option('--get_filter', type=bool, default=True)
@click.option('--filter_size', type=float, default=20)
@click.option('--filter_type', type=str, default='tophat')
@click.option('--randoms_from_gal', type=bool, default=False)

def run_density_splitter(handle,
                        tracer_file,
                        centres_file,
                        nrandoms,
                        box_size,
                        is_matter,
                        dmin,
                        dmax,
                        nrbins,
                        get_monopole,
                        get_rmu,
                        get_spi,
                        get_filter,
                        filter_size,
                        filter_type,
                        randoms_from_gal):

    ds = DensitySplitter(handle=handle,
                        tracer_file=tracer_file,
                        centres_file=centres_file,
                        nrandoms=nrandoms,
                        box_size=box_size,
                        is_matter=is_matter,
                        dmin=dmin,
                        dmax=dmax,
                        nrbins=nrbins,
                        get_monopole=get_monopole,
                        get_rmu=get_rmu,
                        get_spi=get_spi,
                        get_filter=get_filter,
                        filter_size=filter_size,
                        filter_type=filter_type,
                        randoms_from_gal=randoms_from_gal)


if __name__ == '__main__':
    run_density_splitter()