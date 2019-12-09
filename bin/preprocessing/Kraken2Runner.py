import os

import time
import datetime
from multiprocessing import Pool


from os.path import *
from os import listdir
import PreprocessingConfig as config

def run_kraken2_program(kraken_savedir, input_file_pairs, threads=20, db_path='', other_flags='',
                        report_postfix='.report.out', output_postfix='.output.out', unclassified_postfix='.unclassified_#.fastq'):
    """
    Runs the kraken2 program and classifies each read in a pair-end sequencing.
    It is assumed that kraken2 and its dependencies are available.

    :param kraken_savedir:     output folder for the kraken results
    :param input_file_pairs:   prefix -> pair-end sequence paths
    :param threads:            number of cpu cores to be used for the classification
    :param db_path:            path to the preprocessed and indexed database
    :param other_flags:        other flags for the analyses
    :param output_postfix:     postfix for output files
    :param report_postfix:     postfix for report files
    :return:                   list of kraken report files
    """
    kraken_report_files = {}
    i = 0
    for act_prefix in input_file_pairs:
        kraken_report_files[act_prefix] = join(kraken_savedir, act_prefix + report_postfix)
        kraken_cmd = '/gfs/progs/kraken2/kraken2  \
                  --report {0} \
                  --output {1} \
                  --db {2} \
                  --threads {3} \
                  --unclassified-out {4} \
                  {5} \
                  --paired {6}'.format(join(kraken_savedir, act_prefix + report_postfix),
                                       join(kraken_savedir, act_prefix + output_postfix),
                                       db_path,
                                       threads,
                                       join(kraken_savedir, act_prefix + unclassified_postfix),
                                       ' '.join(other_flags),
                                       ' '.join(input_file_pairs[act_prefix]))

        try:
            os.mkdir(kraken_savedir)
            print("Directory " , kraken_savedir ,  " created")
        except FileExistsError:
            pass
            #print("Directory " , output_dir ,  " already exists")

        os.system("echo " + kraken_cmd)
        os.system(kraken_cmd)
    return kraken_report_files


def get_file_pairs(dataset_name):
    file_directory = config.get_qcreport_source_path(dataset_name, True)

    if os.path.exists(file_directory):
        file_list = [f for f in listdir(file_directory)
                        if isfile(join(file_directory, f)) and not 'junk' in f]
    else:
        file_list = []

    pair_dict = {}
    for filename in file_list:
        tags = filename.split(".")
        prefix = tags[0].split("_")[0]
        if prefix in pair_dict:
            pair_dict[prefix].append(join(file_directory, filename))
        else:
            pair_dict[prefix] = [join(file_directory, filename)]

    return pair_dict



#ToTest: RampelliS_2015

# where to find files? //gfs/data/curated_metagenomes_preprocessed/KarlssonFH_2013
# input_file_pairs: Dictionary where the keys are the ids and the values are the filenames of read pairs.


