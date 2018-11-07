#!/usr/bin/env python
'''
# Reading NGS Fastq files by codon rather than nucleotide.
# Parse through each individual reads
'''
import pysam  
import pandas as pd
from Bio import SeqIO
import sys
import os
import time
import itertools

#==============================================================================
# Input
#==============================================================================
# =============================================================================
# input_bam_file = sys.argv[1]
# ref_fasta = sys.argv[2] 
# inputgenbankfile = sys.argv[3]  # only for extracting recoded position and fasta length
# =============================================================================

Input_Bam_File = '/Volumes/ChinLab/Team_Genome/TS_Scripts/example_pool/20kb_wt_sC2_with_SI.fa_sorted.bam'
ref_fasta = '/Volumes/ChinLab/Team_Genome/TS_Scripts/example_pool/20kb_wt_sC2_with_SI.fa'
Input_GenBank_File = '/Volumes/ChinLab/Team_Genome/TS_Scripts/example_pool/20kb_wt_sC2_with_SI.gb'

# Example = v11
recoding_scheme = [{'target': 'TCG', 'recode': 'AGT'},
                   {'target': 'TCA', 'recode': 'AGT'},
                   {'target': 'TAG', 'recode': 'TAA'}]

# NGS filter
Min_Read_Length = 200 # filter out read with read length lower than indicated
MAPQ = 0 # filter out reads with mapping quality less than indicated, set to zero because of the long synthetic insertion

# =============================================================================
# Define output and functions
# =============================================================================
Fasta_Filename = os.path.splitext(os.path.basename(ref_fasta))[0] # isolate filename           
output_image_filename = ref_fasta + '.svg'
output_csv_name = ref_fasta + '.csv'

# get a list of codon properties: 0 = codon, 1 = codon start, 2 = codon end, 3 = codon position set; 4 = strand ('f' or 'r'), 5= new genbank annotation
def get_codon(genbank, target_codon):  
    codon_pos_list = []
    if type(target_codon) != list:
        target_codon = [target_codon] 
    for feature in genbank.features:
        if feature.type == "CDS":
            start=feature.location.start.position
            end=feature.location.end.position
            sense=feature.strand   
            if sense == 1:                                                  
                orf = genbank.seq[start:end].upper()   
                codon = [orf[i:i+3] for i in range(0,len(orf),3)]             
                codon_pos = [{'codon': str(x),
                              'start': (i*3)+start,
                              'end': (i*3)+start+3,
                              'position': set(range((i*3)+start, (i*3)+start+3)), # set(range(start, end))
                              'strand': 'f',
                              'recode':[recoding_scheme[a].get('recode') for a in range(len(recoding_scheme)) if str(x) == recoding_scheme[a].get('target')][0]} for i, x in enumerate(codon) if x in target_codon]     
            else:
                orf = genbank.seq[start:end].reverse_complement().upper()              
                codon = [orf[i:i+3] for i in range(0,len(orf),3)]                          
                codon_pos = [{'codon': str(x),
                              'start': ((len(codon)-i)*3)+start,
                              'end': ((len(codon)-i)*3)+start-3,
                              'position': set(range(((len(codon)-i)*3)+start-3, ((len(codon)-i)*3)+start)), 
                              'strand': 'r',
                              'recode':[recoding_scheme[a].get('recode') for a in range(len(recoding_scheme)) if str(x) == recoding_scheme[a].get('target')][0]} for i, x in enumerate(codon) if x in target_codon] 
            codon_pos_list.append(codon_pos)    
    return list(itertools.chain.from_iterable(codon_pos_list))

def read_bam_by_row(CodonPosition, TargetCodon, RecodedCodon, FastaFilename, InputBamFile, MinReadLength, mapq):
    ''' 
    return output that is a dictionary containing:
                       ({'Position'       : position investigated,
                        'ID list'         : ID of all corespponding ngs read in 'Seq list',
                        'Seq list'        : all the sequence containing the targeted codon,
                        'freq'            : frequency of read at this particular position,
                        'Depth'           : ngs read depth at this position,
                        'NoneRecodedCodon': NoneRecodedCodon})
    
    == Arg ==
    TargetList         : a list of targeted positions (remember to -1 for Python 0-indexing)
        Eg. targetlist = [29,30,31]
    NoneRecodedCodon   : string of the wt codon
    RecodedCodon       : string of the recoded codon
    MinReadLength      : read less than this number will be discarded from further analysis
    MAPQ               : MAPQ quality below this value will be discarded from further analysis

    '''
    codon_record = []
    with pysam.AlignmentFile(InputBamFile,"rb") as bamfile: 
        
        for read in bamfile.fetch(FastaFilename, CodonPosition[0], CodonPosition[-1]):  # take only the seq that contain targeted region, make sure targeted region is arranged by order
            try:  # file 1 sequence number 9 raised ValueError
                sequence = read.query_alignment_sequence # string
                position = read.get_reference_positions() # list
                quality = read.mapping_quality
                length = read.reference_length
            except ValueError:
                pass   
                     
            if len(position) == len(sequence) and length > MinReadLength: # no indel & correct position and min length = 200
                if quality > mapq and len(position) > 0: 
                    try:     
                        target_seq = [sequence[position.index(i)].upper() for i in CodonPosition] # exact index for targeted position for each read 
                        seq = ''.join(target_seq) # Eg ['T','C','G'] are joined as 'TCG' 
                        if seq in [TargetCodon,RecodedCodon]: # check that the seq is either nrc or rc, else is ngs errors and discarded
                            codon_record.append(seq)        
                    except ValueError:  # ValueError raised when there is a deletion in the sequence and targetlist position is not present
                        pass
        try: 
            TotalRecoded = codon_record.count(RecodedCodon)
            freq = round((TotalRecoded/len(codon_record)),2)
        except ZeroDivisionError:
            freq = 0
        output = {'Position'        : CodonPosition[0],
                  'Frequency'       : freq,
                  'Depth'           : len(codon_record),
                  'Target Codon'    : TargetCodon,
                  'Recoded Codon'   : RecodedCodon,
                  'Entries'         : codon_record,}
    return output

if __name__ == '__main__':
    start_time = time.time()
    
    #==============================================================================
    # Calculate codon frequency
    #==============================================================================
    mds = SeqIO.read(Input_GenBank_File, "genbank")                             
    target_codon = [i.get('target') for i in recoding_scheme]                     
    codon_list = get_codon(mds, target_codon)  
    
    codon_frequency = [read_bam_by_row(CodonPosition = sorted(list(i.get('position'))), 
                                       TargetCodon = i.get('codon'), 
                                       RecodedCodon = i.get('recode'), 
                                       FastaFilename = Fasta_Filename,
                                       InputBamFile = Input_Bam_File, 
                                       MinReadLength = Min_Read_Length, 
                                       mapq = MAPQ) for i in codon_list]
    
    # export data to a csv file
    df = pd.DataFrame.from_dict(codon_frequency)
    df.to_csv(output_csv_name)
    
    # backup for manual entry
    df = pd.read_csv(output_csv_name)
    
    #==============================================================================
    # Plot
    #==============================================================================
    position = list(df['Position'])
    depth = list(df['Depth'])
    frequency = list(df['Frequency'])
    
    # get the length of fasta to set x-axis limit
    for record in SeqIO.parse(Input_GenBank_File, 'genbank'):
        fasta_length = len(record.seq)
    
    # Arrange Data to start and end from 0
    position.insert(0,684) # start of ORF = 685 (684 in python)
    position.append(17480) # end = 18485 (18484 in python)
    depth.insert(0,0)
    depth.append(0)
    frequency.insert(0,1)
    frequency.append(1)
    
    #==============================================================================
    # Plot line chart 
    #==============================================================================
    import matplotlib.pyplot as plt
    
    plt.rcParams["figure.figsize"] = [8,2]
    plt.rcParams["figure.dpi"] = 300
    
    fig,(ax1,ax2) = plt.subplots(2,sharex = True)
    
    ax1.plot(position, depth, linestyle = "-", color = "blue", linewidth = 0.5, markerfacecolor="none")
    ax1.axis([684, 17480,0,max(depth)])    # (xmin, xmax, ymin, ymax) xmin = 658 because it is the start of the first gene; xmax = 18485 because it is the end of last gene
    
    ax2.plot(position, frequency,"s", linestyle = "-", color = "red", linewidth = 0.5, markerfacecolor="none", markeredgecolor = "none")  # Draw line connecting dots
    ax2.plot(position, frequency,"s", ms = 2, markerfacecolor="red", markeredgecolor = "none")  # Draw red dots then black dots 
    ax2.axis([684, 17480,0,1])  # (xmin, xmax, ymin, ymax) xmin = 658 because it is the start of the first gene; xmax = 18485 because it is the end of last gene 
    
    # Show and save image as
    fig.savefig(output_image_filename, dpi=300)
    
    # Print processing time
    end_time = time.time() - start_time
    print('Total processing time = {}'.format(time.strftime('%H:%M:%S', time.gmtime(end_time))))
