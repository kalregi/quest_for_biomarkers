from os.path import *
import Preprocessing as Preprocessing
import FastqFileLister as lister
import PreprocessingConfig as config

# Datasets to preprocess
datasetList = config.get_dataset_names()

# Generating basic quality control profile for the datasets
for dataset in datasetList:
    Preprocessing.generate_quality_control_profile(dataset, threads=10)

# Listing fastq files to CSV
file_list_path = "/gfs/data/curated_metagenomes_qc/FileLists"
for dataset in datasetList:
    lister.print_file_list(dataset, dataset+"_filelist.csv", output_file_path=file_list_path)

# Trim fastq files based on config parameters
for dataset_name in datasetList:
    Preprocessing.run_trimmomatic(dataset_name)

# Generating quality control profile for the trimmed data
for dataset_name in datasetList:
    Preprocessing.generate_quality_control_profile(dataset_name=dataset_name, trimmed = True, threads=10)