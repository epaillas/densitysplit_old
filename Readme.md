### Density Split

This set of codes can be used to split the galaxy density field in different quantiles, as in Paillas et al. 2021 (https://arxiv.org/abs/2101.09854)

**1) generate_positions.py**: generates a set of randomly position centres on a simulation box

Input arguments:

* centres_filename: type str, name of the output file where the random positions will be stored
* npositions: type int, number of random positions to be generated
* sampling: type str, how to sample the random positions; 'uniform' for a uniform disitribution across 
            the simulation volume, or 'tracers' to randomly sample from tracer positions 
* boxsize: type float, size of the simulation box
* data_filename: type str, optional, name of the file where tracer positions are stored. Must be a Fortran 90 unformatted file. To convert from an ASCII text file, use ascii_to_unformatted.py

Example of usage:
`python generate_positions.py --centres_filename centres.dat --npositions 10000 --sampling tophat`
    

**2) calculate_filter.py**: calculates the integrated galaxy density at a given distance for each random position

Input arguments:

    --data_filename: type str, name of the file where tracer positions are stored. Must be a Fortran 90 unformatted file. To convert from an ASCII text file, use ascii_to_unformatted.py
    --filter_filename: type str, name of the file where the filtered densities will be stored
    --boxsize: type float, size of the simulation box
    --ngrid: type int, number of cells to divide the simulation box to ease up calculations (100 cells for a 1 Gpc box is a good number)
    --filter_type: type str, can be either 'tophat', or 'gaussian'
    --filter_size: size of the filter (e.g. 15 Mpc/h in Paillas et al. 2021
    --dmin: type float, optional,  minimum distance for the filter calculation. Defaults to zero. I suggest to leave this unchanged unless you know what you're doing
    --dmax: type float, optional, maximum distance for the filter calculation. Defaults to filter_size for a tophat filter, or 5*filter_size for a Gaussian filter. I suggest to leave this unchanged unless you know what you're doing
    --qperp: type float, optional, defaults to 1, adds geometrical distortion parameter along the perpendicular direction. 
    --qpara: type float, optional, defaults to 1, adds geometrical distortion parameter along the line of sight


**3) split_positions.py**: splits the positions in different quantiles by using the chosen filter

Input arguments:

* centres_filename: type str, name of the file where random positions are stored.
* filter_filename: type str, name of the file where the filtered densities are stored (calculated from calculate_filter.py)
* quantiles: type int, number of quantiles for the split

Example of usage
`python split_positions.py --centres_filename centres.dat --filter_filename filter.dat --quantiles 5`


For inquiries, please contact me at epaillas@astro.puc.cl