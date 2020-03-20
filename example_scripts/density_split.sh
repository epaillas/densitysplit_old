#!/bin/bash

#---- Universal settings ----
handle=
tracer_file=
nrandoms=
box_size=
steps=
is_matter=
dmin=
dmax=
nrbin=

#-------- Run (do not modify below this point) --------
python $HOME/code/density_splitter/density_split.py \
--handle "$handle" \
--tracer_file "$tracer_file" \
--nrandoms "$nrandoms" \
--box_size "$box_size" \
--steps "$steps" \
--is_matter "$is_matter" \
--dmin "$dmin" \
--dmax "$dmax" \
--nrbin "$nrbin"
