import os
import re
from os.path import join
from IPython.display import HTML
import pandas as pd

# functions
def run_qc(fastq_input_folder, quality_control_output, quality_control_summary_output, default_fastq_ext='fastq', threads=20):
    """
    Runs FastQC program on all sequencing file in the given folder and generates the summarized reports of QC.
    It is assumed that FastQC and MultiQC is installed and available from the environment.

    :param fastq_input_folder:              Input folder of short read sequences. It is assumed that the extension of the files are 'default_fastq_ext' ( i.e. .fastq, .fq, etc).
    :param quality_control_output:          Output folder for individual fastq files
    :param quality_control_summary_output:  Output folder for the summarized output (MultiQC)
    :param default_fastq_ext:               Default extension of sequencing files ( i.e. .fastq, .fq, etc).
    :param threads:                         Number of CPU cores to be used during the analysis
    :return:
    """
    raw_fastqc_cmd = 'fastqc --outdir {0} \
     -t {1} \
     {2} \
    '.format(quality_control_output, threads, join(fastq_input_folder, '*' + default_fastq_ext))
    print('Running raw fastqc: %s' % raw_fastqc_cmd)
    os.system(raw_fastqc_cmd)

    # doing multiQC
    raw_multiqc_cmd = 'multiqc -f --outdir {0} \
    {1} \
    '.format(quality_control_summary_output, join(quality_control_output, '*'))
    print('Running raw multiqc: %s' % raw_multiqc_cmd)
    os.system(raw_multiqc_cmd)


def get_trimmomatic_cmd_default(trimmed_output, fastq_files, trimmomatic_path, trimmomatic_parameters):
    """
    Preprocesses a pair-end sequencing data given the parameters and writes the preprocessed reads into a new folder.
    It is assumed that the read pair has the same prefix and the forward one contains the substring 'R1' while the reverse one contains 'R2'.

    :param trimmed_output:          Output folder for the preprocessed reads. Preprocessed fastq files will contain the 'trimmed' substring
    :param fastq_files:             List of input fastq files (one forward and one reverse sequencing short read)
    :param trimmomatic_path:        Absolute path to the trimmomatic program
    :param trimmomatic_parameters:  Parameters for the trimmomatic
    :return:                        A runnable command
    """
    input_filename_forward = os.path.basename(fastq_files[0])
    input_filename_reverse = os.path.basename(fastq_files[1])

    (base_filename_forward, ext_forward) = os.path.splitext(input_filename_forward)
    (base_filename_reverse, ext_reverse) = os.path.splitext(input_filename_reverse)

    output_fastq_file_forward = join(trimmed_output, base_filename_forward + '.trimmed' + ext_forward)
    output_fastq_file_reverse = join(trimmed_output, base_filename_reverse + '.trimmed' + ext_reverse)

    output_fastq_junk_file_forward = join(trimmed_output, base_filename_forward + '.junk.trimmed' + ext_forward)
    output_fastq_junk_file_reverse = join(trimmed_output, base_filename_reverse + '.junk.trimmed' + ext_reverse)

    trimmomatic_cmd = 'java -jar {0} PE -phred33 \
                   {1} {2}                   \
                   {3} {4}                   \
                   {5} {6}                  \
                   {7} \
    '.format(trimmomatic_path, fastq_files[0], fastq_files[1],
             output_fastq_file_forward, output_fastq_junk_file_forward,
             output_fastq_file_reverse, output_fastq_junk_file_reverse,
             trimmomatic_parameters)
    return trimmomatic_cmd


def get_prefix_from_filename(filename):
    """
    Gets the prefix (i.e. id) from a filename of a fastq file.
    :param filename:  filename of a fastq file
    :return:
    """
    # longest prefix before R1 or R2
    pair_end_match = re.compile("(\_R[1|2])|(\-[1|2])\.")
    match = pair_end_match.search(filename)
    try:
        prefix = filename[0:match.start()]
    except AttributeError:
        prefix = None
    return prefix

def get_file_pairs_orientations(filename):
    if '_R1_' in filename:
        act_orientation = 'forward'
    if '_R2_' in filename:
        act_orientation = 'reverse'
        
    return act_orientation
    


def get_file_pairs(save_path, suffix_filter_to = '_fastqc.html'):
    """
    Tries to assign the pair-end read related files into a common id.    
    :param save_path:         Folder to scan
    :param suffix_filter_to:  Consider only files ending with string
    :return: 
    """
    filtered_filenames = sorted([f for path, dfd, filenames in os.walk(save_path) for f in filenames if f.endswith(suffix_filter_to)])
    prefix_to_pair = {}
    
    print('Bllalalaldsgffdfd')
    
    for filename in filtered_filenames:
        act_orientation = None
        act_prefix = get_prefix_from_filename(filename)
        print(act_prefix)
        act_orientation = get_file_pairs_orientations(filename)
        print(act_orientation)

        if act_prefix in prefix_to_pair:
            prefix_to_pair[act_prefix][act_orientation] = join( filename)
        else:
            prefix_to_pair[act_prefix] = {'forward' : '', 'reverse' : ''}
            prefix_to_pair[act_prefix][act_orientation] = join( filename)
    return prefix_to_pair


def get_qc_html_table(output_dir):
    """
    Generates a HTML linkable output dataframe
    :param output_dir: 
    :return: 
    """
    pd.set_option('display.max_colwidth', -1)
    seq2output = get_file_pairs(output_dir, suffix_filter_to = '_fastqc.html')
    act_columns = ['ID' , 'Forward Reads', 'Reverse Reads']
    qc_results_table = pd.DataFrame.from_dict(seq2output, orient='index').reset_index()
    qc_results_table.columns = act_columns
    qc_results_table['Forward Reads'] = qc_results_table['Forward Reads'].apply(lambda x:
                                                        '<a href="https://localhost:8888/files/{0}">{1}</a>'.format(join(output_dir, x), x))
    qc_results_table['Reverse Reads'] = qc_results_table['Reverse Reads'].apply(lambda x:
                                                        '<a href="https://localhost:8888/files/{0}">{1}</a>'.format(join(output_dir, x), x))
    return qc_results_table






def get_illumina_pairs(save_path, filtering_string='junk'):
    """
    Assigns the sequencing pairs to an ID (sequencing)
    :param save_path:         Input folder finding the pairs in
    :param filtering_string:  Ignore those files that contain that substring. Typically 'junk' parts should be ignored
    :return:                  Dictionary where the keys are the ids and the values are the filenames of read pairs.
    """
    raw_files = sorted([f for path, dfd, filenames in os.walk(save_path) for f in filenames if os.path.splitext(f)[1] == '.fastq'
                        and filtering_string not in f])
    raw_prefix = list(set([get_prefix_from_filename(filename) for filename in raw_files]))
    prefix_to_pair = {} 
    for filename in raw_files:
        act_prefix = get_prefix_from_filename(filename)      
        
        if act_prefix in prefix_to_pair:
            prefix_to_pair[act_prefix].append(join(save_path, filename))
        else:
            prefix_to_pair[act_prefix] = [join(save_path, filename)]
    return prefix_to_pair


def run_kraken2_program(kraken_savedir, input_file_pairs, threads=20, db_path='', other_flags='', report_postfix='.report.out', output_postfix='.output.out'):
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
    for act_prefix in input_file_pairs:
        kraken_report_files[act_prefix] = join(kraken_savedir, act_prefix + report_postfix)
        kraken_cmd = '/gfs/progs/kraken2/kraken2  \
                  --report {0} \
                  --output {1} \
                  --db {2} \
                  --threads {3} \
                  {4} \
                  --paired {5}'.format(join(kraken_savedir, act_prefix + report_postfix),
                                       join(kraken_savedir, act_prefix + output_postfix),
                                       db_path,
                                       threads,
                                       ' '.join(other_flags),
                                       ' '.join(input_file_pairs[act_prefix]))
        print(kraken_cmd)
        os.system(kraken_cmd)
    return kraken_report_files

def get_bracken_commands(kraken_report_files, bracken_output_folder, db_path, report_levels, bracken_db_length, bracken_threshold):
    """
    Generates runnable bracken commands in order to estimate the abundance of species in the sample.

    :param kraken_report_files:
    :param bracken_output_folder:  Output folder for storing the bracken output
    :param db_path:                Path to the indexed database
    :param report_levels:          List of taxonomic levels to be reported.
    :param bracken_db_length:      Length of k-mers to be used for the analysis
    :param bracken_threshold:      Threshold for
    :return:
    """

    print('Running bracken!')
    bracken_cmds = []
    for act_sample_id, act_report in kraken_report_files.items():
        for act_level in report_levels:
            act_bracken_output = join(bracken_output_folder, act_sample_id + '_' + act_level + '.bracken_report.txt')
            bracken_cmd = 'bracken \
              -d {0} \
              -i {1} \
              -o {2} \
              -l {3} \
              -r {4} \
              -t {5} \
              '.format(db_path, act_report, act_bracken_output, act_level, bracken_db_length, bracken_threshold)

            bracken_cmds.append(bracken_cmd)
    return bracken_cmds
