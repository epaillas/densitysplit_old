import numpy as np
import matplotlib.pyplot as plt
from python_tools.densitysplitter import DensitySplitter
import os
import sys
import argparse


parser = argparse.ArgumentParser(description='Density splitter.')

parser.add_argument('--handle', type=str, required=True)
parser.add_argument('--tracer_file', type=str, required=True)
parser.add_argument('--nrandoms', type=int, required=True)
parser.add_argument('--box_size', type=float, required=True)
parser.add_argument('--steps', type=str, default='1,2')

args = parser.parse_args()  


ds = DensitySplitter(handle=args.handle, tracer_file=args.tracer_file,
                     nrandoms=args.nrandoms, box_size=args.box_size,
                     steps=args.steps)

