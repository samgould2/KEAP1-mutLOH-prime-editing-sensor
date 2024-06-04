# KEAP1-mutLOH-prime-editing-sensor

This repository includes code for generating a prime editing sensor library targeting VUS, annotated driver mutations, and germline variants in KEAP1.

1. The notebook used for the generation of the prime editing sensor library is *keap1_library_generation.ipynb*
    - For more information about using PEGG, visit https://pegg.readthedocs.io

2. The scripts used for analysis of raw sequencing files are located in the *analysis_scripts* folder
    1. *fastq_join_cluster.sh* is the script for performing a fastq-join on paired reads on a slurm computing cluster
    2. *sensor_extraction.py* is the python script for counting pegRNAs, and sorting sensor sites into different fastq files based on their pegRNA identity. There is an associated slurm .sh script as well.
    3. *crispress_analysis.py* is the python script for performing running crispresso2 on each sensor site in each sample. There is an associated slurm .sh script as well.
    4. *crispresso_analsysis_aggregation.py* is the script for aggregating the crispresso2 editing quantification into a single table for each sample.

3. *main_figure.ipynb* is a notebook for generating all of the main figure panels.
4. *supplementary_figures.ipynb* is a notebook for generating all of the supplementary figures.
    - All data in the above 2 notebooks are referenced locally so that they can be run upon download of the repository.
        - Crispresso editing data is stored in the *crispresso* folder
        - MAGeCK data is stored in the *mageck* folder
        - Raw counts data is stored in the *count_tables* folder
        - Reporter validation flow data is stored in the *ARE-GFP_and_PEAR_flow_data* folder
