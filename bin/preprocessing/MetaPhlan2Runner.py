import os

import time
import datetime
from multiprocessing import Pool


from os.path import *
from os import listdir
import PreprocessingConfig as config

def run_metaphlan2_from_bowtie(metaphlan_savedir, file_dict, threads=12, other_flags='',
                        report_postfix='.profile.txt'):
    metphlan_report_files = {}
    i = 0
    for act_prefix in file_dict:
        if i != 0:
            metphlan_report_files[act_prefix] = join(metaphlan_savedir, act_prefix + report_postfix)
            metaphlan_cmd = 'python3 /gfs/progs/metaphlan2_201910/metaphlan2/metaphlan2.py  \
                      {0} \
                      --input_type bowtie2out \
                      --nproc {1} \
                      {2} \
                      > {3}'.format(file_dict[act_prefix],
                                           threads,
                                           ' '.join(other_flags),
                                           join(metaphlan_savedir, act_prefix + report_postfix),
                                           )

            try:
                os.mkdir(metaphlan_savedir)
                print("Directory " , metaphlan_savedir ,  " created")
            except FileExistsError:
                pass
                #print("Directory " , output_dir ,  " already exists")
            print(metaphlan_cmd)

            os.system("echo " + metaphlan_cmd)
            os.system(metaphlan_cmd)
        else:
            pass
        i = i+1
    return metphlan_report_files

#python3 /gfs/progs/metaphlan2_201910/metaphlan2/metaphlan2.py
#    /gfs/data/curated_metagenomes_metaphlan2/ZellerG_2014/ERR479598.bowtie2out.txt --input_type bowtie2out

# python3 /gfs/progs/metaphlan2_201910/metaphlan2/metaphlan2.py
# ERR479032.bowtie2out.txt --input_type bowtie2out --nproc 24

def get_bowtie_files(dataset_name):
    file_directory = os.path.join("/gfs/data/curated_metagenomes_metaphlan2", dataset_name)

    if os.path.exists(file_directory):
        file_list = [f for f in listdir(file_directory)
                        if isfile(join(file_directory, f)) and f.endswith("bowtie2out.txt")]
    else:
        file_list = []

    file_dict = {}
    for filename in file_list:
        tags = filename.split(".")
        prefix = tags[0].split("_")[0]
        file_dict[prefix] = os.path.join(file_directory, filename)

    return file_dict


def run_metaphlan2_program(metaphlan_savedir, input_file_pairs, threads=12, other_flags='',
                        report_postfix='.profile.txt', bowtie2out_postfix='.bowtie2out.txt'):

    metphlan_report_files = {}
    i = 0
    for act_prefix in input_file_pairs:
        if i != 0:
            metphlan_report_files[act_prefix] = join(metaphlan_savedir, act_prefix + report_postfix)
            metaphlan_cmd = 'python3 /gfs/progs/metaphlan2_201910/metaphlan2/metaphlan2.py  \
                      {0} \
                      --input_type fastq \
                      --nproc {1} \
                      --bowtie2out {2}\
                      {3} \
                      > {4}'.format(','.join(input_file_pairs[act_prefix]),
                                           threads,
                                           join(metaphlan_savedir, act_prefix + bowtie2out_postfix),
                                           ' '.join(other_flags),
                                           join(metaphlan_savedir, act_prefix + report_postfix),
                                           )

            try:
                os.mkdir(metaphlan_savedir)
                print("Directory " , metaphlan_savedir ,  " created")
            except FileExistsError:
                pass
                #print("Directory " , output_dir ,  " already exists")
            print(metaphlan_cmd)

            os.system("echo " + metaphlan_cmd)
            os.system(metaphlan_cmd)
        else:
            pass
        i = i+1
    return metphlan_report_files




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


