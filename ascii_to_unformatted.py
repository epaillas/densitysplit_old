import numpy as np
import sys
from scipy.io import FortranFile
import argparse

parser = argparse.ArgumentParser(description='')
parser.add_argument('--input_filename', type=str, required=True)
parser.add_argument('--output_filename', type=str, required=True)
args = parser.parse_args()

print('Converting ASCII file to an unformatted Fortran 90 file.')
print('input_filename: {}'.format(args.input_filename))
print('output_filename: {}'.format(args.output_filename))

data = np.genfromtxt(args.input_filename)
f = FortranFile(args.output_filename, 'w')

nrows, ncols = np.shape(data)

f.write_record(nrows)
f.write_record(ncols)
f.write_record(data)

f.close()
