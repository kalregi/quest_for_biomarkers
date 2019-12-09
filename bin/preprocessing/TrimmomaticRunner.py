
import os
from os.path import join
import PreprocessingConfig as config

TRIMMOMATIC_PATH = config.get_trimmomatic_path()

def run_trimmomatic_pair_end(path_to_fastq_files, output_dir, n_process=None,
                             illuminaclip:tuple=("TruSeq3-PE.fa", 2, 15, 10), leading=9,
                             trailing=3, slidingwindow:tuple=(6, 20), minlen=36, path_to_config=None):
    """
    Trimmomatic processing of pair-end sequences

    :param path_to_fastq_files: Path to file pair to trim.
    :param output_dir: The desired output directory where trimmed files go.
    :param n_process: Number of processes.
    :param illuminaclip: Finding full and part adapters or contaminants. Usage: (<fasteWithAdaptersEtc>,
        <seed mismatches>, <palindrom clip threshold>, <simple clip threshold>, <minAdapterLength>, <keepBothReads>)
    :param leading: Cutting those bases from the beginning that are lower quality then it is given.
    :param trailing: Same as LEADING but from the end.
    :param slidingwindow: We can set the window size, that determine the length of a window in the average quality
        is measured. If the average quality is lower than the threshold, that part of the sequence will be cut out.
        Usage: (<windowSize>, <requiedQuality>).
    :param minlen: Removing those sequences that are shorter than the given threshold.
    """

    command = _get_trimmomatic_pair_end_command(path_to_fastq_files, output_dir, n_process, illuminaclip, leading,
                                                trailing, slidingwindow, minlen, path_to_config)
    os.system(command)


def _get_trimmomatic_pair_end_command(path_to_fastq_files, output_dir, n_process=None,
                                     illuminaclip:tuple=("TruSeq3-PE.fa", 2, 15, 10), leading=9,
                                     trailing=3, slidingwindow:tuple=(6, 20), minlen=36, path_to_config=None):

    if path_to_config is not None:
        config.set_config_file(path_to_config)
        global TRIMMOMATIC_PATH
        TRIMMOMATIC_PATH = config.get_trimmomatic_path()

    if n_process is None:
        n_process = config.get_trimmomatic_thread_count()

    (output_fastq_file_pair1, output_fastq_junk_file_pair1, output_fastq_file_pair2, output_fastq_junk_file_pair2) =\
            _get_output_files(path_to_fastq_files, output_dir)


    parameters = _get_parameters(illuminaclip, leading, trailing, slidingwindow, minlen)

    trimmomatic_command = ["java -jar", TRIMMOMATIC_PATH, "PE -threads", str(n_process),
                           "-phred33",
                           path_to_fastq_files[0], path_to_fastq_files[1],
                           output_fastq_file_pair1, output_fastq_junk_file_pair1,
                           output_fastq_file_pair2, output_fastq_junk_file_pair2,
                           parameters]
    return ' '.join(trimmomatic_command)

def _get_parameters(illuminaclip, leading, trailing, slidingwindow, minlen):

    def _join_parameters(parameters, name, values):
        if values is not None:
            try:
                parameters += name + ":" + ":".join(str(v) for v in values) + " "
            except TypeError:
                parameters += name + ":" + str(values) + " "
        return parameters

    if os.path.split(illuminaclip[0])[0] == "":
        trimmomatic_path_dir = os.path.split(TRIMMOMATIC_PATH)
        illuminaclip = list(illuminaclip)
        illuminaclip[0] = join(trimmomatic_path_dir[0], 'adapters', illuminaclip[0])

    parameters = ""
    parameters = _join_parameters(parameters, "ILLUMINACLIP", illuminaclip)
    parameters = _join_parameters(parameters, "LEADING", leading)
    parameters = _join_parameters(parameters, "TRAILING", trailing)
    parameters = _join_parameters(parameters, "SLIDINGWINDOW", slidingwindow)
    parameters = _join_parameters(parameters, "MINLEN", minlen)

    return parameters

def _get_output_files(path_to_fastq_files, output_dir):

    input_filename_pair_1 = os.path.basename(path_to_fastq_files[0])
    input_filename_pair_2 = os.path.basename(path_to_fastq_files[1])
    (base_filename_pair_1, ext_pair_1) = os.path.splitext(input_filename_pair_1)
    (base_filename_pair_2, ext_pair_2) = os.path.splitext(input_filename_pair_2)

    output_extensions = config.get_trimmomatic_output_extensions()

    os.makedirs(output_dir, exist_ok=True)
    """
    try:
        os.mkdir(output_dir)
        print("Directory " , output_dir ,  " created")
    except FileExistsError:
        print("Directory " , output_dir ,  " already exists")
    """

    output_fastq_file_pair1 = \
        join(output_dir,
             base_filename_pair_1 + output_extensions['pair_1_ext'] + ext_pair_1)

    output_fastq_file_pair2 = \
        join(output_dir,
             base_filename_pair_2 + output_extensions['pair_2_ext'] + ext_pair_2)

    output_fastq_junk_file_pair1 =\
        join(output_dir,
             base_filename_pair_1 + output_extensions['junk_pair_1_ext'] + ext_pair_1)

    output_fastq_junk_file_pair2 = \
        join(output_dir,
             base_filename_pair_2 + output_extensions['junk_pair_2_ext'] + ext_pair_2)

    return (output_fastq_file_pair1, output_fastq_junk_file_pair1,
            output_fastq_file_pair2, output_fastq_junk_file_pair2)