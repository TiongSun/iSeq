#!/usr/bin/env python3

import sys

#==============================================================================
# Input
#==============================================================================

input_csv = sys.argv[1] # columns = position, A, C, G, T
output_filename = input_csv + '_consensus_seq.txt'

#==============================================================================
# Write consensus
#==============================================================================

consensus = ''

with open(input_csv) as f:
    for line,row in enumerate(f):
        if line >2: # get rid of titles
            readrow = [int(float(x)) for x in row.split()[:5]]
            position = readrow[0]
            rowdic = {'A': readrow[1], 'C':readrow[2],'G':readrow[3],'T':readrow[4]}
            sumall = sum((rowdic.values()))
            base = max(rowdic.keys(), key = (lambda key: rowdic[key]))
            if sumall > 0:
                consensus += base

with open(output_filename,'w') as r:
    r.write(consensus)
