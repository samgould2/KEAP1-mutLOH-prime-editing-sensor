#importing necessary packages
from collections import defaultdict
import gzip
from Bio.Align import PairwiseAligner
import numpy as np
#import matplotlib.pyplot as plt
import Bio.Seq
import pandas as pd
#import matplotlib.pyplot as plt
#import seaborn as sns
import Bio.Seq
import Bio.SeqIO
from Bio.SeqRecord import SeqRecord
import os
import sys
from pathlib import Path


#-------Loading in system args
# Parse command-line arguments to read in files of interest
if len(sys.argv) != 5:
    print("Usage: python3 sensor_extraction.py <input_df> <R1_FILE> -o <folder_name>")
    sys.exit(1)

#example usage

input_df = pd.read_csv(Path(sys.argv[1]))
R1_FILE= Path(sys.argv[2])
folder_name = str(sys.argv[4])



#----------GLOBAL VARIABLES

GZ = False #True when files are gzipped

# Minimum average read quality
MIN_QUALITY = 30 #originally set to 30
tevo_front = 'CGCGGTTCTATCTAG'
scaff_front = 'GCACCGAGTCGGTGC'


#-------Functions

def fastq_reader(fname, gz=False):
    _open = gzip.open if gz else open
    proc_read = (lambda line: line.strip().decode()) if gz else (lambda line: line.strip())
    proc_qual = (lambda line: line.strip()) if gz else (lambda line: line.strip().encode('utf-8'))
    with _open(fname, 'r') as f:
        for i, line in enumerate(f):
            if i % 4 == 0:
                identifier = proc_read(line)
            elif i % 4 == 1:
                read = proc_read(line)
            elif i % 4 == 3:
                quality = proc_qual(line)
                yield identifier[-20:], read, quality #protospacer (index2), read, quality score

def joined_fastq_peg_identifier(r1_true, input_df):
    proto = r1_true[:19]


    if (scaff_front in r1_true) & (tevo_front in r1_true):
        s = r1_true.find(scaff_front) + len(scaff_front)
        e = r1_true.find(tevo_front)

        ext = r1_true[s:e]

        j = input_df[(input_df['Protospacer']==('G' + proto)) & (input_df['RTT_PBS']==ext)]
        if len(j)>0:
            id = j['pegRNA_id'].values[0]
            out_class = 'correct_id'
        else:
            id = None
            out_class = 'unidentified'
    
    else:
        out_class = 'unaligned'
        id = None

    return out_class, id

'''Prepare fastq records.'''
def to_IOSeq_rec(sensor_seq, i, q1):
    qname = 'read_' + str(i)
    record = SeqRecord(Bio.Seq.Seq(sensor_seq), id=qname, name=qname, description='', dbxrefs=[])

    #add quality score to the record
    qual = (np.frombuffer(q1, dtype=np.uint8) - 33)
    record.letter_annotations["phred_quality"] = qual[0:len(sensor_seq)]
    
    return record

def extraction_filtration(folder_name, input_df, R1_FILE, breakpoint=False):

    """
    Takes in reads and returns dataframe containing pegRNA counts AND sensor outcomes.

    Also returns summary of outcomes including:
    - low quality (and which of the reads are low quality (or if all are low quality))
    - no extension match
    - no protospacer match
    - decoupled extension-protospacer
    - correct identification
    """

    os.mkdir(folder_name) #make the new directory (needs to be a non-existent folder)

    #write 500 empty fastq files to a new folder
    for i in list(input_df['pegRNA_id']):
        f_name = i + '.fastq'
        with open(folder_name + '/' + f_name, 'w') as fp:
            pass

    #------initialize a dataframe for holding the pegRNA counts and sensor outcomes
    d1 = pd.DataFrame(dict(zip(['pegRNA_id'], [list(input_df['pegRNA_id'])])))
    cols = ['total', 'extracted', 'recombined']
    z = np.zeros((len(cols), len(input_df)))
    d2 = pd.DataFrame(dict(zip(cols, z)))
    df1 = pd.concat((d1, d2), axis=1).set_index('pegRNA_id')

    #------initialize a dataframe for holding metadata about the identification of sensors
    outcomes = ['good_quality', 'low_quality', 'unaligned', 'unidentified', 'correct_id']
    #unaligned means cant find tevopreQ1 or the end of the scaffold for identifying the extension
    #unidentified means that there's no match to the sensor or extension

    outcomes_count = np.zeros(len(outcomes))
    class_df = pd.DataFrame(dict(zip(['classification', 'count'],[outcomes, outcomes_count]))).set_index('classification')

    #iterating through the reads...
    for i, ((proto1, r1, q1)) in enumerate(fastq_reader(R1_FILE, gz=GZ), 1):
        #r1 = fastq-joind read (need to reverse complement to get in correct orientation)

        #first check the quality
        qual1 = (np.frombuffer(q1, dtype=np.uint8) - 33).mean()
        if qual1 > MIN_QUALITY:
            quality_out = 'good_quality'
        else:
            quality_out = 'low_quality'

        class_df.loc[quality_out, 'count']+=1

        #if quality is poor, abort
        if quality_out != 'good_quality':
            continue

        #else, if quality is good, continue
        else:

            #taking reverse complement to put it in correct orientation
            r1_true = str(Bio.Seq.Seq(r1).reverse_complement())
        

            #identify the pegRNA
            out_class, id = joined_fastq_peg_identifier(r1_true, input_df)
            class_df.loc[out_class, 'count']+=1

            if out_class=='correct_id':
                #count the pegRNA
                df1.loc[id, 'total'] += 1

                #and now try to pull out the sensor
                #first checking for recombination based on the last 8 nt of the read
                sens_anchor = input_df.loc[input_df['pegRNA_id']==id, 'sensor_wt'].values[0][-8:]
                sens_read_anchor = r1_true[-8:]
                if sens_anchor==sens_read_anchor:
                    df1.loc[id, 'extracted'] += 1

                    #not doing anything fancy; just taking the whole expected 60 nt sensor
                    sensor_seq = r1_true[-60:]
                    
                    #add the sensor read to the appropriate fastq file (numbered according to the pegRNA index)
                    out_file = folder_name + '/' + str(id) + '.fastq'

                    record = to_IOSeq_rec(sensor_seq, i, q1)

                    with open(out_file, 'a') as fq:
                        #and write it to the appropriate file
                        Bio.SeqIO.write(record, fq, 'fastq')


                else:
                    df1.loc[id, 'recombined'] += 1

        if breakpoint != False:
            if i>breakpoint:
                break

    return df1, class_df

count_df, class_df = extraction_filtration(folder_name, input_df, R1_FILE, breakpoint=False)

#save the counts and classification dataframes
count_df.to_csv(folder_name + '_count_df.csv')
class_df.to_csv(folder_name + '_classification_df.csv')

