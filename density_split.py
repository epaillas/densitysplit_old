from python_tools.densitysplitter import DensitySplitter
import argparse

parser = argparse.ArgumentParser(description='Density splitter.')

parser.add_argument('--handle', type=str, required=True)
parser.add_argument('--tracer_file', type=str, required=True)
parser.add_argument('--nrandoms', type=int, required=True)
parser.add_argument('--box_size', type=float, required=True)
parser.add_argument('--steps', type=str, default='1,2')
parser.add_argument('--is_matter', type=bool, default=False)
parser.add_argument('--dmin', type=float, default=0)
parser.add_argument('--dmax', type=float, default=100)
parser.add_argument('--nrbins', type=int, default=50)

args = parser.parse_args()  


ds = DensitySplitter(handle=args.handle, tracer_file=args.tracer_file,
                     nrandoms=args.nrandoms, box_size=args.box_size,
                     steps=args.steps, is_matter=args.is_matter,
                     dmin=args.dmin, dmax=args.dmax, nrbins=args.nrbins)

