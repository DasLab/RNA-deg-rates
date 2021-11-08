#!/usr/bin/python

import numpy as np
from arnie.pfunc import pfunc
from arnie.free_energy import free_energy
from arnie.bpps import bpps
import arnie.utils as utils
import argparse

parser = argparse.ArgumentParser( description='output p_unp for all sequences in a file' )
parser.add_argument( 'filename',type=str,help='filename with sequences');
parser.add_argument( '--output','-o',type=str,help='csv-formatted p_unp values',default='output.txt');
parser.add_argument( '--pkg','-p',type=str,help='package to use, e.g., vienna or contrafold',default='vienna');
args = parser.parse_args()
print args.filename


fid = open( args.output,'w' );


pkg = args.pkg;
lines = open( args.filename ).readlines()
for i,line in enumerate(lines):
    sequence = line.rstrip()
    sequence = sequence.upper()
    sequence = sequence.replace( 'T','U');
    print sequence

    #sequence = sequence[1:101]
    
    bp_matrix = bpps(sequence, package=pkg)

    print 'Doing line ',i+1, 'of',len(lines)
    
    # We sum all values over one axis to get total base pairing probability per nucleotide.
    p_net_base_pairing = np.sum(bp_matrix, axis=0)
    p_net_base_pairing.reshape(1,-1)
    
    print 1.0-p_net_base_pairing

    for k in range(len(p_net_base_pairing)):
        print k, p_net_base_pairing[k]
        if k > 0: fid.write( ',' )
        fid.write( '%f' % (1.0-p_net_base_pairing[k]) )
    fid.write('\n')

fid.close()
