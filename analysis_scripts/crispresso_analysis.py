import pandas as pd
import os
import sys
from pathlib import Path


#-------Loading in system args
# Parse command-line arguments to read in files of interest
if len(sys.argv) != 3:
    print("Usage: python3 crispresso_analysis.py <input_df> <sample_name>")
    sys.exit(1)

#example usage
df = pd.read_csv(Path(sys.argv[1]))
sample_name = Path(sys.argv[2])

#--------run it for all of the pegRNAs
#for i in range(len(p53)):

for i, val in df.iterrows():

    peg_id = val['pegRNA_id']
    wt = val['sensor_wt']
    edited = val['sensor_alt']
    wt_qwc = '5-55'
    edit_qwc = '5-55'
    proto = val['Protospacer']

    #NEED TO MODIFY TO AVOID THE non-targeting??
    #going to try without for now...should just break in some other way

    #generating crispresso command
    f_path = f"./{sample_name}/{peg_id}.fastq"
    start = f"CRISPResso --fastq_r1 {f_path} "
    #removed --discard_indel_reads; leaving this in there...
    p2 = f"--amplicon_seq {wt} --expected_hdr_amplicon_seq {edited} --guide_seq {proto} --quantification_window_coordinates {wt_qwc},{edit_qwc} --suppress_report --suppress_plots --plot_window_size 15 --exclude_bp_from_left 0 --exclude_bp_from_right 0" #--discard_indel_reads

    #output folder = crispresso_output
    out_path = f" -o ./crispresso_output/{sample_name}" 

    command = start+ p2 + out_path

    #running the crispresso command
    #apparently this is bad practice (but subprocess module didn't work...)
    os.system(command)