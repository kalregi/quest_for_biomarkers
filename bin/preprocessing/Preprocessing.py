
import os

import time
import datetime
from multiprocessing import Pool


from os.path import *
from os import listdir
import TrimmomaticRunner as TrimmomaticRunner
import PreprocessingConfig as config

def _timestamp():
    return datetime.datetime.fromtimestamp(time.time()).strftime('%Y-%m-%d_%H:%M:%S')

def run_preprocessing(dataset_name):
    generate_quality_control_profile(dataset_name)
    run_trimmomatic(dataset_name)
    generate_quality_control_profile(dataset_name, True)


def generate_quality_control_profile(dataset_name, trimmed=False, threads=22):
    fastqc(dataset_name, trimmed, threads)
    multiqc(dataset_name, trimmed)

def get_fastqc_command(dataset_name, trimmed=False, threads=22):
    command = None
    output_dir = config.get_qcreport_destination_path(dataset_name, trimmed)
    filedirectory = config.get_qcreport_source_path(dataset_name, trimmed)
    file_list = list_files(dataset_name, trimmed)
    commands = []
    if len(file_list) != 0:
        try:
            os.mkdir(output_dir)
            #print("Directory " , output_dir ,  " created")
        except FileExistsError:
            pass
            #print("Directory " , output_dir ,  " already exists")
        file_list_size = 500
        for part in range(0, round(len(file_list)/file_list_size)+1):
            part_start = part*file_list_size
            part_end = min([(part+1)*file_list_size, len(file_list)])
            fastqc_command = ["fastqc", " ".join((join(filedirectory, f) for f in file_list[part_start:part_end])), "-q", "-o", output_dir,
                          "-t", str(threads)]

            #print(' '.join(fastqc_command))
            commands.append(' '.join(fastqc_command))
    else:
        print(dataset_name, "has no files")
    return commands

def fastqc(dataset_name, trimmed=False, threads=22):
    print("fastqc on", dataset_name)
    commands = get_fastqc_command(dataset_name=dataset_name, trimmed=trimmed, threads=threads)
    if commands is not None:
        i = 1
        for command in commands:
            if i == 1:
                os.system(command)
            i = i+1


def get_multiqc_command(dataset_name, trimmed=False):
    command = None
    input_dir = config.get_qcreport_destination_path(dataset_name, trimmed)
    if os.path.exists(input_dir):
        multiqc_filename = config.get_multiqc_output_destination_file_name(dataset_name, trimmed)
        multiqc_output_directory = config.get_multiqc_destination_path(dataset_name, trimmed)

        try:
            os.mkdir(multiqc_output_directory)
            #print("Directory " , multiqc_output_directory ,  " created")
        except FileExistsError:
            pass
            #print("Directory " , multiqc_output_directory ,  " already exists")

        multiqc_command = ["multiqc -f -q", input_dir, "-o", multiqc_output_directory, "-n", multiqc_filename]
        #print(' '.join(multiqc_command))
        command = ' '.join(multiqc_command)
    else:
        print(input_dir, "does not exist")
    return command

def multiqc(dataset_name, trimmed=False):
    print("multiqc on", dataset_name)
    command = get_multiqc_command(dataset_name=dataset_name, trimmed=trimmed)
    if command is not None:
        print("multiqc command running...")
        os.system(command)

def list_files(dataset_name, trimmed=False):

        file_directory = config.get_qcreport_source_path(dataset_name, trimmed)
        if os.path.exists(file_directory):
            file_list = [f for f in listdir(file_directory) if isfile(join(file_directory, f))
                       and check_if_file_processable(f, trimmed)]
        else:
            file_list = []
        return file_list


def check_if_file_processable(filename, trimmed=False):
        trimmed_tag_1 = config.get_trimmomatic_output_extensions()["pair_1_ext"]
        trimmed_tag_2 = config.get_trimmomatic_output_extensions()["pair_2_ext"]
        junk_tag_1 = config.get_trimmomatic_output_extensions()["junk_pair_1_ext"]
        junk_tag_2 = config.get_trimmomatic_output_extensions()["junk_pair_2_ext"]
        if junk_tag_1 in filename or junk_tag_2 in filename:
            return False
        if trimmed:
            return filename.endswith(".fastq"+trimmed_tag_1+".gz") \
                   or filename.endswith(trimmed_tag_1+".fastq") \
                   or filename.endswith(".fastq"+trimmed_tag_2+".gz") \
                   or filename.endswith(trimmed_tag_2+".fastq")

        else:
            return trimmed_tag_1 not in filename and trimmed_tag_1 not in filename \
                   and (filename.endswith("fastq.gz") or filename.endswith(".fastq"))


def order_files_to_paires(file_list):
    pair_list = []
    for f in file_list:
        splitted = f.split('.')
        name = splitted[0]
        ext = ".".join(splitted[1:])
        if name[-2:] == "_1":
            other_pair = str(name[:-2] ) +"_2." + ext
            if (None, other_pair) in pair_list:
                pair_list[pair_list.index((None, other_pair))] = (f, other_pair)
            else:
                pair_list.append((f, None))
        if name[-2:] == "_2":
            other_pair = str(name[:-2] ) +"_1." + ext
            if (other_pair, None) in pair_list:
                pair_list[pair_list.index((other_pair, None))] = (other_pair ,f)
            else:
                pair_list.append((None ,f))
    for (pair_1, pair_2) in pair_list:
        if pair_1 is None:
            #raise ValueError("File_list argument does not contain pair for '" + pair_2 + "'")
            pass
        if pair_2 is None :
            #raise ValueError("File_list argument does not contain pair for '" + pair_1 + "'")
            pass
    return pair_list

def _pool_helper(paramdict):
    TrimmomaticRunner.run_trimmomatic_pair_end(
        path_to_fastq_files=paramdict["path_to_fastq_files"],
        output_dir=paramdict["output_dir"],
        illuminaclip=paramdict["illuminaclip"],
        leading=paramdict["leading"],
        trailing=paramdict["trailing"],
        slidingwindow=paramdict["slidingwindow"],
        minlen=paramdict["minlen"]
    )

def run_trimmomatic(dataset_name):
        if config.is_trimmable(dataset_name) and config.is_pair(dataset_name):
            pair_list = order_files_to_paires(list_files(dataset_name, False))
            directory_path = config.get_source_dir(dataset_name)
            trimmomatic_parameters = config.get_trimommatic_parameters(dataset_name)
            print("Run trimmomatic on:", pair_list)
            parameter_dict_list = []
            for files in pair_list:
                parameter_dict_list.append({
                    "path_to_fastq_files": [join(directory_path, f) for f in files],
                    "output_dir": config.get_trimmomatic_destination_path(dataset_name),
                    "illuminaclip": trimmomatic_parameters["illuminaclip"],
                    "leading": trimmomatic_parameters["leading"],
                    "trailing": trimmomatic_parameters["leading"],
                    "slidingwindow": trimmomatic_parameters["slidingwindow"],
                    "minlen": trimmomatic_parameters["minlen"]
                })
            with Pool(5) as p:
                p.map(_pool_helper, parameter_dict_list)
