import pandas as pd
import os
import sys
from pathlib import Path


#-------Loading in system args
# Parse command-line arguments to read in files of interest
if len(sys.argv) != 4:
    print("Usage: python3 crispresso_analysis_aggregation.py <input_df> <folder_name> <quant_zero>")
    sys.exit(1)

#example usage
p53 = pd.read_csv(Path(sys.argv[1]))
sample_name = str(sys.argv[2])
quant_zero = pd.read_csv(Path(sys.argv[3]))

samp_name2 = sample_name.split("/")[1]

#runnnig it to create a large dataframe
fp = './crispresso_output' #path to output folders


#--------run it for all of the pegRNAs
rows = []
for i, val in p53.iterrows():

    peg_num = val['pegRNA_id']
    output_folder_x = os.listdir(fp + '/' + sample_name + f"/CRISPResso_on_{peg_num}")
    peg_id = peg_num

    if 'CRISPResso_quantification_of_editing_frequency.txt' in output_folder_x:
        quant = pd.read_csv(fp + '/' + sample_name + f"/CRISPResso_on_{peg_num}" + '/' + 'CRISPResso_quantification_of_editing_frequency.txt', sep='\t')
        quant['peg_id'] = peg_id

    else:
        quant = quant_zero.copy()
        quant['peg_id'] = peg_id
    
    rows.append(quant)

output = pd.concat(rows)  
output.to_csv(f"./crispresso_output/{samp_name2}_crispresso_aggregated.csv", index=False)