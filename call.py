#!/usr/bin/env python3

import sys
import pandas as pd
from Bio import SeqIO

#==============================================================================
# Input
#==============================================================================
input1 = sys.argv[1] # input_csv from igvtool wigs
input2 = sys.argv[2]# ref fasta

input_csv = input1 # columns = position, A, C, G, T
refseq = input2

output_filename = input2 + '_call.csv'
depth_checker_filename = input2 + '_depth_checker.csv'

#==============================================================================
# Calculation
#==============================================================================
refpos = [list(i.seq) for i in SeqIO.parse(refseq,'fasta')][0]

output, poslist = [],[]
min_read_depth  = 1
count = 0

with open(input_csv) as f:
    for num,row in enumerate(f):
        if num >2: # get rid of titles
            count += 1
            readrow = [int(float(x)) for x in row.split()[:5]]
            position = readrow[0]
            ref = refpos[position-1].upper() 
            
            rowdic = {'A': readrow[1], 'C':readrow[2],'G':readrow[3],'T':readrow[4]}
            sumall = sum((rowdic.values()))
            poslist.append([position, sumall])

# no read depth

            if position != count:
                try:
                    out = {'Position': str(poslist[0][-2]) + ' to ' +str(position),
                           'Mutation': 'Read Depth = 0',
                           'Ref Base Rate':  '0',
                           'Read Base Rate': '0',
                           'Read Depth': '0'}
                    output.append(out)
                    count = position
                except IndexError:
                    out = {'Position': '1 to ' +str(position),
                           'Mutation': 'Read Depth = 0',
                           'Ref Base Rate':  '0',
                           'Read Base Rate': '0',
                           'Read Depth': '0'}
                    output.append(out)
                    count = position
                    
                
            elif sumall == 0:
                out = {'Position': position,
                       'Mutation': 'Read Depth = 0',
                       'Ref Base Rate':  '0',
                       'Read Base Rate': '0',
                       'Read Depth': sumall}
                output.append(out)

# read depth more than 0
            
            else:           
                base = max(rowdic.keys(), key = (lambda key: rowdic[key]))             
                depth_ref = rowdic[ref]
                depth_base = rowdic[base]              
                try: 
                    percentage_base = depth_base/sumall
                except ZeroDivisionError:
                    percentage_base = 0                    
                try:
                    percentage_ref = depth_ref/sumall
                except ZeroDivisionError:
                    percentage_ref = 0

# low read depth
                    
                if sumall < min_read_depth: 
                    out = {'Position': position,
                           'Mutation': 'Read Depth <' + str(min_read_depth),
                           'Ref Base Rate':  ref + ' = ' + str(round(percentage_ref,2)),
                           'Read Base Rate': base + ' = ' + str(round(percentage_base,2)),
                           'Read Depth': sumall}
                    output.append(out)

# sufficient read depth, recode max base and mutation
                    
                elif ref != base:
                    mut = ref + ' to ' + base
                    out = {'Position': position,
                           'Mutation': mut,
                           'Ref Base Rate':  ref + ' = ' + str(round(percentage_ref,2)),
                           'Read Base Rate': base + ' = ' + str(round(percentage_base,2)),
                           'Read Depth': sumall}
                    output.append(out)

df = pd.DataFrame.from_dict(output)
try: 
    df = df[['Position','Mutation','Ref Base Rate', 'Read Base Rate','Read Depth']].set_index('Position')
except KeyError:
    pass
df.to_csv(output_filename)

# Depth checker
depth_checker = []
for i in range(1,len(poslist)):
    if poslist[i-1][1] > (poslist[i][1]* 5):
        percent = (1 - (poslist[i][1]/poslist[i-1][1]))*100
        depth_checker.append({'Position': poslist[i][0],
                              'Note': 'Read depth drops by {} percent'.format(round(percent,0))})
    elif (poslist[i-1][1])*2 < poslist[i][1]:
        try:
            percent = ((poslist[i][1]/poslist[i-1][1])-1)*100
            depth_checker.append({'Position': poslist[i][0],
                                  'Note': 'Read depth increases by {} percent'.format(round(percent,0))})
        except ZeroDivisionError:
            pass

df2 = pd.DataFrame.from_dict(depth_checker)
try: 
    df2 = df2[['Position', 'Note']].set_index('Position')
except KeyError:
    pass
df2.to_csv(depth_checker_filename)        
        
    
