import yaml
from os.path import join

CONFIG_FILE = yaml.load(open("config_preprocessing.yml"))

def set_config_file(path_to_config):
    global CONFIG_FILE
    CONFIG_FILE = yaml.load(open("config_preprocessing.yml"))


# SOURCE DIR

def get_source_dir(dataset_name: str):
    return join(CONFIG_FILE["source_dir"],
                CONFIG_FILE["trimmomatic"]["dataset"][dataset_name]["name"])

# DATASET NAMES

def get_dataset_names():
    dataset_names = CONFIG_FILE["trimmomatic"]["dataset"]
    return [k for k in dataset_names.keys() if k != "default"]

# QC REPORT

def get_qcreport_fastq_thread_count():
    return CONFIG_FILE["qc_report"]["threads"]


def get_qcreport_source_path(dataset_name: str, trimmed: bool = False):
    if trimmed:
        return get_trimmomatic_destination_path(dataset_name)
    else:
        return get_source_dir(dataset_name)


def get_qcreport_destination_path(dataset_name: str, trimmed: bool = False):
    if trimmed:
        return join(CONFIG_FILE["qc_report"]["destination_path"],
                    CONFIG_FILE["qc_report"]["outputfolder"]["trimmed"],
                    CONFIG_FILE["trimmomatic"]["dataset"][dataset_name]["name"])
    else:
        return join(CONFIG_FILE["qc_report"]["destination_path"],
                    CONFIG_FILE["qc_report"]["outputfolder"]["original"],
                    CONFIG_FILE["trimmomatic"]["dataset"][dataset_name]["name"])


def get_multiqc_destination_path(dataset_name: str, trimmed=False):
    return get_qcreport_destination_path(dataset_name, trimmed)


def get_multiqc_output_destination_file_name(dataset_name: str, trimmed=False):
    if trimmed:
        return join(CONFIG_FILE["qc_report"]["destination_path"],
             CONFIG_FILE["trimmomatic"]["dataset"][dataset_name]["name"] + "_" +
             CONFIG_FILE["qc_report"]["outputfolder"]["trimmed"] +
             ".index.html")
    else:
        return join(CONFIG_FILE["qc_report"]["destination_path"],
             CONFIG_FILE["trimmomatic"]["dataset"][dataset_name]["name"] + "_" +
             CONFIG_FILE["qc_report"]["outputfolder"]["original"] +
             ".index.html")


# TRIMMOMATIC


def get_trimommatic_parameters(dataset_name: str):
    config = CONFIG_FILE["trimmomatic"]
    # addig default values for non-existing parameters
    for def_param_name, def_param_value in config["dataset"]["default"].items():
        if def_param_name not in config["dataset"][dataset_name]:
            config["dataset"][dataset_name][def_param_name] = def_param_value
    config = config["dataset"][dataset_name]
    #del config["name"]
    return config

def is_pair(dataset_name: str):
    return get_trimommatic_parameters(dataset_name)["pair"]

def is_trimmable(dataset_name: str):
    return get_trimommatic_parameters(dataset_name)["process"]

def get_trimmomatic_path():
    return CONFIG_FILE["trimmomatic"]["trimmomatic_path"]


def get_trimmomatic_destination_path(dataset_name: str):
    return join(CONFIG_FILE["trimmomatic"]["destination_path"],
                CONFIG_FILE["trimmomatic"]["dataset"][dataset_name]["name"],
                CONFIG_FILE["trimmomatic"]["outputfolder"])


def get_trimmomatic_thread_count():
    return CONFIG_FILE["trimmomatic"]["threads"]


def get_trimmomatic_output_extensions():
    return CONFIG_FILE["trimmomatic"]["output_extension"]