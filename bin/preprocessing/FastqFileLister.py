
import PreprocessingConfig as config
from os.path import *
from os import listdir

def print_file_list(dataset_name, output_file_name, output_file_path=None, delimiter=";", with_header=True, ):
    file_dict = order_files_to_dict(dataset_name)

    if output_file_path is None:
        output_file = open(output_file_name, "w")
    else:
        output_file = open(join(output_file_path, output_file_name),"w")
    content = []
    if with_header:
        header = delimiter.join(["KEY","DATASET","PAIR_1","PAIR_2","SINGLE","IS_PAIR_END","PATH"])
        content.append(header)
    for k, v in file_dict.items():
        row = delimiter.join([v["key"], v["dataset"], v["pair_1"], v["pair_2"], v["single"], str(v["is_pair_end"]), v["path"]])
        content.append(row)
    output_file.write("\n".join(content))

def order_files_to_dict(dataset_name, file_list=None):
    if file_list is None:
        file_list = list_files(dataset_name)
    file_dict = {}
    for filename in file_list:
        key = get_key(filename)
        if key not in file_dict:
            file_dict[key] = {
                "key": key,
                "dataset": dataset_name,
                "pair_1" : "",
                "pair_2" : "",
                "single" : "",
                "is_pair_end" : "",
                "path" : config.get_source_dir(dataset_name)
            }
        pair_end = check_if_file_pair_end(filename)
        file_dict[key]["is_pair_end"] = pair_end
        if pair_end:
            if is_first_pair(filename):
                file_dict[key]["pair_1"] =  filename
            else:
                file_dict[key]["pair_2"] =  filename
        if check_if_file_single_end(filename):
            file_dict[key]["single"] = filename
    return file_dict


def get_key(filename):
    base_filename = filename.split('.')[0]
    return base_filename.split('_')[0]


def is_first_pair(filename):
    base_filename = filename.split('.')[0]
    return base_filename.split('_')[-1] == '1'


def is_second_pair(filename):
    base_filename = filename.split('.')[0]
    return base_filename.split('_')[-1] == '2'


def list_files(dataset_name):
    file_directory = config.get_qcreport_source_path(dataset_name)

    file_list = [f for f in listdir(file_directory) if isfile(join(file_directory, f))
                   and check_if_file_processable(f)]
    return file_list


def check_if_file_pair_end(filename):
    return ('_1.' in filename or '_2.' in filename) and \
           (filename.endswith(".fastq.gz") or filename.endswith(".fastq"))


def check_if_file_single_end(filename):
    single_tag = config.get_trimmomatic_output_extensions()["single_ext"]
    return "_1" not in filename and "_2" not in filename\
           and (filename.endswith(".fastq.gz") or filename.endswith(".fastq"))


def check_if_file_processable(filename):
    return filename.endswith(".fastq.gz") or filename.endswith(".fastq")