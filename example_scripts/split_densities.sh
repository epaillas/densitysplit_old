#!/bin/bash

set -e

for NMOCK in $(seq -f "%03g" 1 2)
do

gal_den_file=den_cats/Real/Galaxies_HOD_"$NMOCK"_z0.57_Real.gal_den.unf
dm_den_file=den_cats/Real/Galaxies_HOD_"$NMOCK"_z0.57_Real.DM_den.unf
has_velocity=True
handle=den_cats/Galaxies_HOD_"$NMOCK"_z0.57_Real


python $HOME/code/density_splitter/python_tools/split_densities.py \
--gal_den_file "$gal_den_file" \
--dm_den_file "$dm_den_file" \
--has_velocity "$has_velocity" \
--handle "$handle"

done
