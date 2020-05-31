import numpy as np
import argparse

parser = argparse.ArgumentParser(description='')

parser.add_argument('--handle', type=str)
parser.add_argument('--is_matter', type=int, default=0)
args = parser.parse_args()  

fmt = 2*'%15.5f '

# read galaxy monopole data
if args.is_matter:
    ext_in = '.CCF_DM_monopole'
else:
    ext_in = '.CCF_gal_monopole'
fname = args.handle + ext_in
data = np.genfromtxt(fname)

# xi_r
r_for_xi = data[:,0]
delta = data[:,1]
if args.is_matter:
    ext_out = '.DM_delta_r'
else:
    ext_out = '.gal_xi_r'
fout = args.handle + ext_out
cout = np.c_[r_for_xi, delta]
np.savetxt(fout, cout, fmt=fmt)

# integrated xi_r
r_for_xi = data[:,0]
int_delta = data[:,2]
if args.is_matter:
    ext_out = '.DM_int_delta_r'
else:
    ext_out = '.gal_int_xi_r'
fout = args.handle + ext_out
cout = np.c_[r_for_xi, int_delta]
np.savetxt(fout, cout, fmt=fmt)

# v_r
r_for_v = data[:,0]
v_r = data[:,3]
if args.is_matter:
    ext_out = '.DM_v_r'
else:
    ext_out = '.gal_v_r'
fout = args.handle + ext_out
cout = np.c_[r_for_v, v_r]
np.savetxt(fout, cout, fmt=fmt)

# sv_los
r_for_v = data[:,0]
sv_los = data[:,4]
if args.is_matter:
    ext_out = '.DM_sv_los'
else:
    ext_out = '.gal_sv_los'
fout = args.handle + ext_out
cout = np.c_[r_for_v, sv_los]
np.savetxt(fout, cout, fmt=fmt)





