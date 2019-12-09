# -*- coding: UTF-8 -*-
# Copyright 2016 by Balázs Ligeti (obalasz@gmail.com), xxxx xxx. All rights reserved.
import os
import copy
from os.path import join, splitext
import re
from MicrobiomeDataManagement import MicrobiomeData
from DBCreator import BackgroundDatabase
from BowtieAligner import BowtieSequenceAligner


class PhageAnalysis(object):
    """ Class for analizing the phage profiles in microbiome data """

    _metagenome_id_recognizer = re.compile(r"(([E|S]R[R|S])(\d+))|((MGM)(\d+))", flags=re.IGNORECASE)

    def __init__(self, microbiome_data=None, background_database=None, work_directory=None, trimmomatic_path=None, metaphlan2_path=None,
                 python2_path="/gfs/progs/anaconda2/bin/python2"):
        """
        Initialize phage analysis class.
        :param background_database: A background database, that represents the phage sequences and its related data
        #python2 path is needed for the metaphlan
        :param work_directory: Path to the directory, where the results are stored
        """
        print("Creating PhageAnalysis class!")
        self.microbiome_data = microbiome_data  # type: MicrobiomeData
        self.background_database = background_database  # type: BackgroundDatabase
        self.trimmomatic_path = trimmomatic_path
        self.metaphlan2_path = metaphlan2_path
        self.base_name_taxon = "aggregated_results_taxon_id.tsv"
        self.base_name_genome = "aggregated_results_genome_id.tsv"
        self.python2_path = python2_path

        if work_directory is not None:
            self.work_directory = work_directory
        else:
            self.work_directory = background_database.savedir
        self.temp_directory = join(self.work_directory, "temp")
        if not os.path.exists(self.temp_directory):
            print('MicrobiomeDataProcessor: Creating directory: %s' % self.temp_directory)
            os.makedirs(self.temp_directory)

    def _get_list_of_analized_samples(self, alignment_save_directory):
        """ Finds the samples that have been processed. It is assumed, that the name contains the
        the sample ID
        :param alignment_save_directory:  The directory where the final results are stored
        :return: List of sample files and the mapping between the id and the filename that contains the alignment information
        """
        print("Collecting processed samples!")
        all_files = [f for f in os.listdir(alignment_save_directory) if os.path.isfile(os.path.join(
            alignment_save_directory, f))]
        sample_ids = list(self.microbiome_data.data["sample_id"])
        sample_id2file = {}
        processed_sample_ids = []
        for filename in all_files:
            for sample_id in sample_ids:
                if sample_id in filename and ".sorted" in filename and os.path.getsize(os.path.join(
                        alignment_save_directory, filename)) > 0:
                    processed_sample_ids.append(sample_id)
                    if sample_id in sample_id2file:
                        print("Id has already contains in the dictionary")
                        print("current files: ", sample_id2file[sample_id])
                        a = None.strip()
                    else:
                        sample_id2file[sample_id] = filename
                    break
        print("Number of processed samples: %d" % len(processed_sample_ids))
        return processed_sample_ids, sample_id2file

    def _copy_raw_sequence_file_into_tempory_workdir(self, sample_id, is_soft_copy=True):
        """
        It copies the raw (probably compressed) sequence data file into to temporary directory for further processing
        :param sample_id: the identifier of the sample to work with
        :return: path the to the copied raw sequence data
        """
        sample_temp_directory = self._get_temporary_directory_for_sample(sample_id)
        sample_paths = self.microbiome_data.get_raw_sequence_data_file_name(sample_id)
        paths_to_sample = []
        for sample_path in sample_paths:
            path_to_sample = join(sample_temp_directory, self.microbiome_data.data.loc[sample_id, MicrobiomeData.add_expected_source_file_key])
            if is_soft_copy:
                # using soft links to avoid copy (using NFS is it not neccessery)
                copy_command = "ln -s " + sample_path + " " + join(sample_temp_directory, os.path.split(sample_path)[1])
            else:
                copy_command = "cp " + sample_path + " " + sample_temp_directory
            paths_to_sample.append(path_to_sample)
            # print("Copy the sample %s from \n\t'%s'\n to \n\t'%s'" % (sample_id, sample_path, path_to_sample))
            print(copy_command)
            os.system(copy_command)
        return paths_to_sample

    def _get_temporary_directory_for_sample(self, sample_id):
        """
        Generate a temporary directory for a sample. If the directory is not exist, then it will be created.
        :param sample_id: identifier for the sample
        :return: path to the temporary directory
        """
        sample_temp_directory = join(self.temp_directory, sample_id)
        if not os.path.exists(sample_temp_directory):
            os.mkdir(sample_temp_directory)
        return sample_temp_directory

    def _guess_fastq_types_files_guess_run_id_mappings(self, fastq_files):
        """
        Guessing the sample ids -> filename mappings.
        :param fastq_files:  fastq_files match the ids with
        :return:
        """
        sequence_run_to_id = {}
        for fastq_file in fastq_files:
            matches = re.findall(self._metagenome_id_recognizer, fastq_file)
            if matches:
                sra_id = matches[0][0]
                mgrast_id = matches[0][3]
                if sra_id and sra_id in sequence_run_to_id:
                    sequence_run_to_id[sra_id].append(fastq_file)
                elif sra_id and sra_id not in sequence_run_to_id:
                    sequence_run_to_id[sra_id] = [fastq_file]
                elif mgrast_id and mgrast_id in sequence_run_to_id:
                    sequence_run_to_id[mgrast_id].append(fastq_file)
                elif mgrast_id and mgrast_id not in sequence_run_to_id:
                    sequence_run_to_id[mgrast_id] = [fastq_file]

        return sequence_run_to_id

    @staticmethod
    def _get_fastq_files_in_folder(path, gzip_streaming=False):
        """
        Get the fastq files in a folder
        :return:
        """
        if gzip_streaming:
            fastq_files = [file for file in os.listdir(path) if
                           os.path.isfile(join(path, file))
                           and (splitext(join(path, file))[1] == '.fastq' or splitext(join(path, file))[1] == '.fq'
                                or splitext(join(path, file))[1] == '.gz' or splitext(join(path, file))[1] == '.gz'
                                )]
        else:
            fastq_files = [file for file in os.listdir(path) if
                           os.path.isfile(join(path, file)) and (splitext(join(path, file))[1] == '.fastq' or splitext(join(path, file))[1] == '.fq')]
        return fastq_files

    @staticmethod
    def _guess_fastq_types_files__pair_end_vs_single_decision(run_id_to_filename_mapping):
        """
        Guessing whether the list contains pair-end sequences.
        If it contains, then it returns with the filename of pairs as tuple as well as with the "single" sequences as a tuple
        :param run_id_to_filename_mapping:
        :return:
        """
        potential_pair_end_pairs = [run_id for run_id, files in run_id_to_filename_mapping.items() if len(files) > 1]
        all_pair_end_filenames = []
        for run_id in potential_pair_end_pairs:
            pair_end_filenames = run_id_to_filename_mapping[run_id]
            # direct finding assuming the naming convention __.fastq :
            paired_1_file = [filename for filename in pair_end_filenames if '_1.fastq' in filename or '_1.fq' in filename]
            paired_2_file = [filename for filename in pair_end_filenames if '_2.fastq' in filename or '_2.fq' in filename]
            direct_paired_files = paired_1_file + paired_2_file
            if len(direct_paired_files) == 2:
                pair_end_filenames = direct_paired_files

            if len(pair_end_filenames) > 2:
                pair_end_filenames = [filename for filename in pair_end_filenames if 'singleton' not in filename]

            if len(pair_end_filenames) == 2:
                all_pair_end_filenames.append(pair_end_filenames)
            else:
                print("Filtering was unsuccessful")
        if len(all_pair_end_filenames) == 1:
            all_pair_end_filenames = all_pair_end_filenames[0]
        single_ids = set(list(run_id_to_filename_mapping.keys())) - set(potential_pair_end_pairs)
        single_end_filenames = [run_id_to_filename_mapping[single_id][0] for single_id in single_ids]
        return [all_pair_end_filenames, single_end_filenames]

    def _guess_fastq_types_files(self, sample_id, pair_end_policy="preferred", single_end_policy="concatenate", gzip_streaming=False):
        """
        Guessing the type of fastq file. (Is pair-end or single-end)
        :param sample_id:
        :param pair_end_policy:
        :param single_end_policy:
        :return:
        """
        sample_specific_temp_directory = self._get_temporary_directory_for_sample(sample_id)
        print("Guessing fastq file format!")
        original_fastq_files = self._get_fastq_files_in_folder(sample_specific_temp_directory,gzip_streaming)
        is_pair_end = None
        fastq_files = None
        if not original_fastq_files:
            print("NO fastq file have been found! Check the input!")
        elif len(original_fastq_files) == 1:
            print("Single file has been found!")
            is_pair_end = False
            fastq_files = [original_fastq_files[0]]
        elif len(original_fastq_files) > 1:
            run_id_to_filename_mapping = self._guess_fastq_types_files_guess_run_id_mappings(original_fastq_files)
            [all_pair_end_filenames, dummy] = self._guess_fastq_types_files__pair_end_vs_single_decision(run_id_to_filename_mapping)
            if 'preferred' in pair_end_policy and all_pair_end_filenames:
                print("Selecting pair-end sequences (if it is available)!")
                is_pair_end = True
                fastq_files = all_pair_end_filenames
            elif "concatenate" in single_end_policy:
                import random
                import string
                new_filename = sample_id + ''.join(random.sample(string.ascii_lowercase, 5)) + ".fastq"
                print("concatenating all the available (single) fastq files into one! The new filename is % s" % new_filename)
                with open(join(sample_specific_temp_directory, new_filename), "w") as new_fastq_file:
                    for fastq_files in original_fastq_files:
                        with open(join(sample_specific_temp_directory, fastq_files)) as fastq_file_input:
                            for line in fastq_file_input:
                                new_fastq_file.write(line)
                is_pair_end = False
                fastq_files = [new_filename]
            else:
                print("No available data for further evaluation!")
        else:
            is_pair_end = False
            fastq_files = [original_fastq_files[0]]
        if fastq_files is None:
            path_to_fastq_files = None
        else:
            path_to_fastq_files = [join(sample_specific_temp_directory, fq_file) for fq_file in fastq_files]

        print('[is_pair_end, fastq_files, path_to_fastq_files, original_fastq_files]:', [is_pair_end, fastq_files, path_to_fastq_files, original_fastq_files])
        return [is_pair_end, fastq_files, path_to_fastq_files, original_fastq_files]

    def _screen_for_phage_sequences_compute_bowtie_alignment(self, bsa, is_pair_end, fastq_file, alignment_save_directory, is_delete_sam=False, force_single_sam_processing=False):
        """
        Runs the bowtie with the appproriate paramters
        :param bsa: bwotie sequence aligner
        :param is_pair_end: whether the sequences are pair-end or not
        :param fastq_file:  list of paths to fastq files
        :param alignment_save_directory: where the alignments are deposited
        :return:
        """
        if is_pair_end:
            # pair-end pipeline
            print("Bowtie . using pair-end pipeline")
            bsa.compute_paired_alignment(self.background_database.bowtie_databases,
                                         self.background_database.bowtie_dir,
                                         alignment_save_directory,
                                         fastq_file[0],
                                         fastq_file[1],
                                         is_delete=is_delete_sam,
                                         force_single_sam_processing=force_single_sam_processing)
        else:
            print("Bowtie Using single-end pipeline")
            bsa.compute_alignment(self.background_database.bowtie_databases,
                                  self.background_database.bowtie_dir,
                                  alignment_save_directory,
                                  fastq_file[0],
                                  is_delete=is_delete_sam)
        print(self.background_database.bowtie_databases)

    def _preprocess_fastq_files_trimmomatic_pair_end(self, path_to_fastq_files, n_process=22, gzip_streaming=False):
        """
        Trimmomatic processing of pair-end sequences
        :param path_to_fastq_files:
        :param n_process:
        :return:
        """
        splitted_path_pair1 = os.path.split(path_to_fastq_files[0])
        splitted_path_pair2 = os.path.split(path_to_fastq_files[1])

        default_extension='fastq'
        extension_1 = os.path.splitext(splitted_path_pair1[1])[1]
        extension_2 = os.path.splitext(splitted_path_pair2[1])[1]
        if gzip_streaming and extension_1 == '.gz' and extension_2 == '.gz':
            default_extension='.gz'

        output_fastq_file_pair1 = join(splitted_path_pair1[0], os.path.basename(splitted_path_pair1[1]) + ".trimmed." + default_extension)
        output_fastq_file_pair2 = join(splitted_path_pair2[0], os.path.basename(splitted_path_pair2[1]) + ".trimmed." + default_extension)

        output_fastq_junk_file_pair1 = join(splitted_path_pair1[0], os.path.basename(splitted_path_pair1[1]) + ".junk.trimmed." + default_extension)
        output_fastq_junk_file_pair2 = join(splitted_path_pair2[0], os.path.basename(splitted_path_pair2[1]) + ".junk.trimmed." + default_extension)
        trimmomatic_path_dir = os.path.split(self.trimmomatic_path)
        os.system('cd ' + trimmomatic_path_dir[0])
        print('cd ' + trimmomatic_path_dir[0])
        os.chdir(join(trimmomatic_path_dir[0], 'adapters'))
        trimmomatic_command = ["java -jar", join('..', trimmomatic_path_dir[1]), "PE -threads", str(n_process), "-phred33",
                               path_to_fastq_files[0], path_to_fastq_files[1],
                               output_fastq_file_pair1, output_fastq_junk_file_pair1, output_fastq_file_pair2, output_fastq_junk_file_pair2,
                               "ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:6:20 MINLEN:36"]
        print(' '.join(trimmomatic_command))
        os.system(' '.join(trimmomatic_command))
        print([output_fastq_file_pair1, output_fastq_file_pair2])
        return [output_fastq_file_pair1, output_fastq_file_pair2]

    def _preprocess_fastq_files_trimmomatic_single(self, path_to_fastq_files, n_processs=22, gzip_streaming=False):
        print("preprocessing fastq files: ", path_to_fastq_files)
        print(path_to_fastq_files)
        splited_path = os.path.split(path_to_fastq_files[0])

        default_extension = 'fastq'
        extension_1 = os.path.splitext(splited_path[1])[1]
        if gzip_streaming and extension_1 == '.gz':
            default_extension = '.gz'

        output_fastq_file = join(splited_path[0], os.path.basename(splited_path[1]) + ".trimmed." + default_extension)
        print(output_fastq_file)
        trimmomatic_path_dir = os.path.split(self.trimmomatic_path)
        os.system('cd ' + trimmomatic_path_dir[0])
        print(trimmomatic_path_dir)
        print('cd ' + trimmomatic_path_dir[0])
        os.chdir(join(trimmomatic_path_dir[0], 'adapters'))
        trimmomatic_command = ["java -jar", join('..', trimmomatic_path_dir[1]), "SE -threads", str(n_processs), "-phred33", path_to_fastq_files[0], output_fastq_file,
                               "ILLUMINACLIP:TruSeq3-SE.fa:2:30:10 HEADCROP:4 SLIDINGWINDOW:6:10 MINLEN:36"]
        print(' '.join(trimmomatic_command))
        os.system(' '.join(trimmomatic_command))
        return [output_fastq_file]

    def _run_metaphlan(self, metaphlan_save_directory, sample_id, is_pair_end, path_to_fastq_files, n_process=23):

        print("Running metaphlan. path to mph: ", self.metaphlan2_path)
        save_save_name = join(metaphlan_save_directory, sample_id + ".txt")
        bowtie2_parsed_output_path = join(metaphlan_save_directory, "mph_bowtie_output." + sample_id + ".txt")
        if os.path.exists(bowtie2_parsed_output_path):
            os.remove(bowtie2_parsed_output_path)
        if os.path.exists(save_save_name):
            os.remove(save_save_name)

        if self.metaphlan2_path is None:
            print("Please provide a metaphlan path!")
            return
        if is_pair_end:
            print("Pair-end pipeline!")
            run_metaphlan = [self.python2_path, self.metaphlan2_path + 'metaphlan2.py', "--input_type fastq --nproc", str(n_process),
                             '--mpa_pkl', join(self.metaphlan2_path, 'db_v20', 'mpa_v20_m200.pkl'), '--bowtie2db',
                             join(self.metaphlan2_path, 'db_v20', 'mpa_v20_m200'), "--bowtie2out",
                             bowtie2_parsed_output_path, path_to_fastq_files[0] + "," + path_to_fastq_files[1], bowtie2_parsed_output_path]
        else:
            print("Single-end pipeline")
            run_metaphlan = [self.python2_path, self.metaphlan2_path + 'metaphlan2.py', "--input_type fastq --nproc", str(n_process),
                             '--mpa_pkl', join(self.metaphlan2_path, 'db_v20', 'mpa_v20_m200.pkl'), '--bowtie2db',
                             join(self.metaphlan2_path, 'db_v20', 'mpa_v20_m200'),
                             "--bowtie2out",
                             bowtie2_parsed_output_path, path_to_fastq_files[0], bowtie2_parsed_output_path]
        print(" ".join(run_metaphlan))
        os.system(" ".join(run_metaphlan))

    @staticmethod
    def _count_lines_fastq_files(path_to_fastq_files):
        print("Counting the lines in the files!")
        line_counts = 0
        for fastq_file in path_to_fastq_files:
            with open(fastq_file) as fin:
                for _ in fin:
                    line_counts += 1
        return line_counts

    def _screen_for_phage_sequences_rename_fastq_files(self, sample_id, fastq_file, original_fastq_files):
        """
        Adding a prefix in front of all the filenames rendering the later processing the matching with id easier

        :param sample_id:
        :param fastq_file:
        :param path_to_fastq_files:
        :param original_fastq_files:
        :return:
        """
        print('Renaming the files ... assuming unix os and path settings ')

        prefix = '.idnd.'
        sample_specific_temp_directory = self._get_temporary_directory_for_sample(sample_id)

        new_fastq_files = []
        new_path_to_fastq_files = []
        # new_original_fastq_files = []
        # print('sample_id %s' % sample_id)
        # print('fastq_files %s' % str(fastq_file))
        # print('path_to_fastq_files %s' % str(path_to_fastq_files))
        # print('original_fastq_files %s' % str(original_fastq_files))

        new_original_fastq_files = copy.deepcopy(original_fastq_files)

        for act_fastq_file in fastq_file:
            new_fastq_file = sample_id + prefix + act_fastq_file
            create_link_command = ['ln -s', join(sample_specific_temp_directory,act_fastq_file), join(sample_specific_temp_directory,new_fastq_file)]
            print(' '.join(create_link_command))
            os.system(' '.join(create_link_command))
            new_fastq_files.append(new_fastq_file)
        new_path_to_fastq_files = [join(sample_specific_temp_directory, fastq_file) for fastq_file in new_fastq_files]

        # print('new sample_id %s' % sample_id)
        # print('new fastq_file %s' % str(new_fastq_files))
        # print('new path_to_fastq_files %s' % str(new_path_to_fastq_files))
        # print('new original_fastq_files %s' % str(new_original_fastq_files))

        return new_fastq_files, new_path_to_fastq_files, new_original_fastq_files




    def _screen_for_phage_sequences(self, sample_id, bsa, alignment_save_directory, metaphlan_save_directory=None, is_running_metaphlan=True,
                                    is_remove_fastq=False, is_preprocessing=True, sample_statistics=None,
                                    pair_end_policy="preferred", single_end_policy="concatenate", n_process=23, gzip_streaming=True, force_single_sam_processing=False, is_delete_sam=True):
        """
        Finding potentially phage sequences in a specific microbiome data .
        :param sample_id: identifier of the microbiome
        :param is_remove_fastq: It is True if one wants to remove the raw fastq file
        """
        from time import gmtime, strftime
        import shutil
        sample_specific_directory = self._get_temporary_directory_for_sample(sample_id)
        sample_statistics_header = ["sample_id", "read_count", "is_pair_end", "is_concatenated", "is_preprocessed", "analysis_start_time"]
        is_preprocessed = False
        analysis_start_time = strftime("%Y-%m-%d %H:%M:%S", gmtime())
        sample_statistics=False
        print("Processing samnple ", sample_id)

        if os.path.exists(sample_specific_directory):
            shutil.rmtree(sample_specific_directory)
        self._copy_raw_sequence_file_into_tempory_workdir(sample_id)
        if gzip_streaming:
            # check for gzip files
            fastqfiles = self._get_fastq_files_in_folder(sample_specific_directory, gzip_streaming)
        else:
            self.microbiome_data.extract_compressed_raw_sequence_data(sample_specific_directory)
        [is_pair_end, fastq_file, path_to_fastq_files, original_fastq_files] = self._guess_fastq_types_files(sample_id, pair_end_policy, single_end_policy, gzip_streaming)
        # here we should 'rename' the fastq files. The name should begin with the sample_id + '.idns'
        [fastq_file, path_to_fastq_files, original_fastq_files] = self._screen_for_phage_sequences_rename_fastq_files(sample_id, fastq_file, original_fastq_files)

        if is_pair_end and is_preprocessing:
            path_to_fastq_files = self._preprocess_fastq_files_trimmomatic_pair_end(path_to_fastq_files, n_process, gzip_streaming=False)
            is_preprocessed = True
        elif not is_pair_end and is_preprocessing:
            path_to_fastq_files = self._preprocess_fastq_files_trimmomatic_single(path_to_fastq_files, n_process, gzip_streaming=False)
            is_preprocessed = True

        if sample_statistics:
            if not os.path.exists(sample_statistics):
                with open(sample_statistics, "w") as stat_out:
                    stat_out.write("\t".join(sample_statistics_header) + "\n")
            line_counts = self._count_lines_fastq_files(path_to_fastq_files)
            is_concatenated = False
            if not is_pair_end and len(original_fastq_files) > 1:
                is_concatenated = True
            records = [str(sample_id), str(int(line_counts / 4)), str(is_pair_end), str(is_concatenated), str(is_preprocessed), str(analysis_start_time)]
            print("Writing sample statistics!")
            with open(sample_statistics, 'a') as stat_out:
                stat_out.write("\t".join(records) + "\n")
        if is_running_metaphlan and metaphlan_save_directory:
            self._run_metaphlan(metaphlan_save_directory, sample_id, is_pair_end, path_to_fastq_files, n_process=n_process)
        if is_pair_end is not None and fastq_file is not None:
            self._screen_for_phage_sequences_compute_bowtie_alignment(bsa, is_pair_end, path_to_fastq_files, alignment_save_directory,is_delete_sam, force_single_sam_processing)
        print(fastq_file, is_pair_end)
        if is_remove_fastq and original_fastq_files:
            [os.remove(join(sample_specific_directory, fastq_file)) for fastq_file in self._get_fastq_files_in_folder(sample_specific_directory)]
            shutil.rmtree(sample_specific_directory)

    def screen_for_phage_sequences(self, body_site='stool', is_remove_fastq=True, is_resume=False, is_running_metaphlan=False,
                                   is_preprocessing=False, sample_statistics=None, n_process=23,
                                   pair_end_policy="preferred", single_end_policy="concatenate", is_pair_end_as_single=False, surpassing_commands=None, is_delete_sam=True):
        """
        Finding potentially phage sequences in microbiome sequence data.
        :param body_site:
        :param is_remove_fastq:
        :param is_resume:
        :param is_running_metaphlan:
        :param is_preprocessing:
        :param sample_statistics:
        :param n_process:
        :param pair_end_policy:
        :param single_end_policy:
        :return:
        """
        from time import gmtime, strftime
        print('MicrobiomeDataProcessor: Start screeening for all dataset!')
        print('surpassing_commands: passing: ', surpassing_commands)

        if surpassing_commands is None :
            bsa = BowtieSequenceAligner(max_cpus=n_process,bowtie_st_surpassing_cmds=None)
        else:
            #surpassing_commands = ["--no-hd", "--no-unal", "--no-sq", "--very-sensitive"]
            bsa = BowtieSequenceAligner(max_cpus=n_process, bowtie_st_surpassing_cmds=surpassing_commands)
        test_mode = False
        temp_sample_ids = []
        if test_mode:
            temp_sample_ids = ["SRS056695", "SRS055298"]
            temp_sample_ids = ["SRS019976", "SRS056695", "SRS055298"]
            temp_sample_ids = ["SRS049712"]
            is_remove_fastq = True

        if self.background_database is not None:
            alignment_save_directory = os.path.join(self.background_database.savedir, self.microbiome_data.dataset_name + "_results", "alignments")
            metaphlan_save_directory = os.path.join(self.background_database.savedir, self.microbiome_data.dataset_name + "_results", "metaphlan")
        else:
            alignment_save_directory = os.path.join(self.work_directory, self.microbiome_data.dataset_name + "_results", "alignments")
            metaphlan_save_directory = os.path.join(self.work_directory, self.microbiome_data.dataset_name + "_results", "metaphlan")
        if self.background_database is not None:
            error_log_file_name = os.path.join(self.background_database.savedir, self.microbiome_data.dataset_name + "_results", "error_log.txt")
        else:
            error_log_file_name = os.path.join(self.work_directory, self.microbiome_data.dataset_name + "_results", "error_log.txt")
        if sample_statistics is None:
            if self.background_database is not None:
                sample_statistics_file = os.path.join(self.work_directory, self.microbiome_data.dataset_name + "_results", "statistics_log.txt")
            else:
                sample_statistics_file = os.path.join(self.work_directory, self.microbiome_data.dataset_name + "_results", "statistics_log.txt")
        else:
            sample_statistics_file = sample_statistics

        try:
            if not os.path.exists(alignment_save_directory):
                print('MicrobiomeDataProcessor: Creating directory: %s' % alignment_save_directory)
                os.makedirs(alignment_save_directory)
            if not os.path.exists(metaphlan_save_directory):
                print('MicrobiomeDataProcessor: Creating directory: %s' % metaphlan_save_directory)
                os.makedirs(metaphlan_save_directory)

        except FileExistsError as inst:
            print("MicrobiomeDataProcessor: - screen_for_phage_sequences (continue): Error,", inst)

        if not os.path.exists(error_log_file_name):
            with open(error_log_file_name, "w") as fout:
                fout.write("#errors:\n")
        if not os.path.exists(alignment_save_directory):
            os.makedirs(alignment_save_directory)
        downloaded_samples = self.microbiome_data.get_downloaded_body_type_samples(body_site)
        print("MicrobiomeDataProcessor: Number of samples: %d" % len(downloaded_samples))
        if is_resume and not test_mode:
            [processed_elements, sample_id2file] = self._get_list_of_analized_samples(alignment_save_directory)
            print(processed_elements)
            samples = list(set(downloaded_samples) - set(processed_elements))
            print("MicrobiomeDataProcessor: Size of dataset: %d" % len(samples))
        elif test_mode:
            samples = temp_sample_ids
        else:
            samples = downloaded_samples
        print("MicrobiomeDataProcessor: Size of dataset: %d" % len(samples))
        for sample_id in samples:
            print("MicrobiomeDataProcessor: Working with the sample: %s" % sample_id)
            try:
                self._screen_for_phage_sequences(sample_id, bsa, alignment_save_directory, metaphlan_save_directory, is_running_metaphlan, is_remove_fastq,
                                                 is_preprocessing, sample_statistics_file, pair_end_policy, single_end_policy, n_process,
                                                 force_single_sam_processing=is_pair_end_as_single, is_delete_sam=is_delete_sam)
            except Exception as se:
                print("MicrobiomeDataProcessor: Error during with the sample: " + str(sample_id) + "; analysis error:", se)
                with open(error_log_file_name, 'a') as fout:
                    fout.write(str(sample_id) + "\t" + strftime("%Y-%m-%d %H:%M:%S", gmtime()) + "\terror: " + str(se) + "\n")
                    raise se

    @staticmethod
    def _add_count_to_dictionary(dictionary, key):
        """ Adding data to a dictionary"""
        if key in dictionary:
            dictionary[key] += 1
        elif key is not None:
            dictionary[key] = 1

    @staticmethod
    def _process_alignment_file_guessing_format(filename):
        """
        Guessing the format of the file by reading its first 10 rows.

        :param filename: input filename
        :return: True if the rows contain 14 record, false otherwise
        """
        data_lengths = []
        with open(filename) as fin:
            for i, line in enumerate(fin):
                data = line.strip().split()
                data_lengths.append(len(data))
                if i > 10:
                    break
        avg_lengths = sum(data_lengths) / len(data_lengths)

        if avg_lengths == 14:
            is_paired = True
        else:
            is_paired = False
        return is_paired

    @staticmethod
    def _process_alignment_file_process_record(data):
        """

        :return:
        """
        read_id = data[0]
        hit_pair1 = data[1]
        try:
            [genome, taxon] = hit_pair1.split(";")
            genome = int(genome)
            taxon = int(taxon)
        except ValueError:
            print("Value error! Invalid distance values: ", data)
            genome = None
            taxon = None
        try:
            coord = int(data[2])
        except ValueError:
            print("Value error! Invalid distance values: ", data)
            coord = -1
        try:
            hit_quality = float(data[5])
        except ValueError:
            print("Value error! Invalid distance values: ", data)
            hit_quality = -1
        return [read_id, genome, taxon, coord, hit_quality]

    @staticmethod
    def _process_alignment_file_lca_single(filename, taxonomy, score_tsh, number_of_hits):
        """
        Processing the alignment file of single, simple alignmetns
        :param filename: path to the file to be processed
        :param taxonomy: taxonomy, that contains the record
        :param score_tsh: alignment having a score smnaller then this tresholds are ignored
        :param number_of_hits: number of hits should be consired i.e. filtering to unique hits only
        :return: aggregated counts
        """
        act_read_id = None
        act_taxa = []
        taxon2count = {}
        with open(filename) as fin:
            for line in fin:
                data = line.strip().split()
                [read_id, genome, taxon, coord, hit_quality] = PhageAnalysis._process_alignment_file_process_record(data)
                if read_id != act_read_id:
                    if 0 < len(act_taxa) <= number_of_hits:
                        lca = taxonomy.find_common_ancestor(act_taxa)
                        PhageAnalysis._add_count_to_dictionary(taxon2count, lca)
                    act_read_id = read_id
                    act_taxa = []
                if taxon is not None and hit_quality > score_tsh:
                    act_taxa.append(taxon)
        return taxon2count

    @staticmethod
    def _process_alignment_file_lca_paired(filename, taxonomy, score_tsh, number_of_hits, is_concordant=True):
        """
        Processing the alignment file of paired alignments
        :param filename: path to the file to be processed
        :param taxonomy: taxonomy, that contains the record
        :param score_tsh: alignment having a score smnaller then this tresholds are ignored
        :param number_of_hits: number of hits should be consired i.e. filtering to unique hits only
        :return: aggregated counts
        """
        act_read_id = None
        act_taxa = []
        taxon2count = {}
        with open(filename) as fin:
            for line in fin:
                data = line.strip().split()
                [read_id1, genome1, taxon1, coord1, hit_quality1] = PhageAnalysis._process_alignment_file_process_record(data[0:7])
                [_, genome2, taxon2, coord2, hit_quality2] = PhageAnalysis._process_alignment_file_process_record(data[7:])
                if read_id1 != act_read_id:
                    if 0 < len(act_taxa) <= number_of_hits:
                        lca = taxonomy.find_common_ancestor(act_taxa)
                        PhageAnalysis._add_count_to_dictionary(taxon2count, lca)
                    act_read_id = read_id1
                    act_taxa = []
                if is_concordant and genome2 == genome1 and abs(coord2 - coord1) < 1500 \
                        and hit_quality1 > score_tsh and hit_quality2 > score_tsh:
                    act_taxa.append(taxon1)
                elif not is_concordant:
                    if genome1 and hit_quality1 > score_tsh:
                        act_taxa.append(taxon1)
                    if genome2 and hit_quality2 > score_tsh:
                        act_taxa.append(taxon2)
        return taxon2count

    @staticmethod
    def _process_alignment_file_lca(taxonomy, sample_id, filename, output_filename_taxon, is_concordant=True, score_tsh=0, number_of_hits=500):
        """
        Generating LCA results for the file
        :param taxonomy:
        :return:
        """
        print("Working with file: ", filename)
        if os.path.getsize(filename) < 1:
            print("Empty file!")
            return
        is_paired = PhageAnalysis._process_alignment_file_guessing_format(filename)
        if not is_paired:
            taxon2count = PhageAnalysis._process_alignment_file_lca_single(filename, taxonomy, score_tsh, number_of_hits)
        else:
            print("Using paired strategy!")
            taxon2count = PhageAnalysis._process_alignment_file_lca_paired(filename, taxonomy, score_tsh, number_of_hits, is_concordant)
        with open(output_filename_taxon, "w") as fout_taxon:
            for taxon, count in taxon2count.items():
                fout_taxon.write(str(taxon) + "\t" + str(count) + "\t" + sample_id + "\n")

    @staticmethod
    def _process_alignment_file_raw_single(filename, score_tsh, number_of_hits):
        """
        Processing the alignment file of single, simple alignmetns
        :param filename: path to the file to be processed
        :param score_tsh: alignment having a score smnaller then this tresholds are ignored
        :param number_of_hits: number of hits should be consired i.e. filtering to unique hits only
        :return: aggregated counts
        """
        act_read_id = None
        act_taxa = []
        act_genomes = []
        taxon2count = {}
        genome2count = {}
        with open(filename) as fin:
            for line in fin:
                data = line.strip().split()
                [read_id, genome, taxon, _, hit_quality] = PhageAnalysis._process_alignment_file_process_record(data)
                if read_id != act_read_id:
                    if 0 < len(act_taxa) <= number_of_hits:
                        for act_taxon in act_taxa:
                            PhageAnalysis._add_count_to_dictionary(taxon2count, act_taxon)
                        for act_genome in act_genomes:
                            PhageAnalysis._add_count_to_dictionary(genome2count, act_genome)
                    act_read_id = read_id
                    act_taxa = []
                    act_genomes = []
                if taxon is not None and hit_quality > score_tsh:
                    act_taxa.append(taxon)
                if genome is not None and hit_quality > score_tsh:
                    act_genomes.append(genome)
        return [genome2count, taxon2count]

    @staticmethod
    def _process_alignment_file_raw_paired(filename, score_tsh, number_of_hits, is_concordant=True):
        """
        Processing the alignment file of paired alignments
        :param filename: path to the file to be processed
        :param score_tsh: alignment having a score smnaller then this tresholds are ignored
        :param number_of_hits: number of hits should be consired i.e. filtering to unique hits only
        :return: aggregated counts
        """
        act_read_id = None
        act_taxa = []
        act_genomes = []
        taxon2count = {}
        genome2count = {}
        with open(filename) as fin:
            for line in fin:
                data = line.strip().split()
                [read_id1, genome1, taxon1, coord1, hit_quality1] = PhageAnalysis._process_alignment_file_process_record(data[0:7])
                [_, genome2, taxon2, coord2, hit_quality2] = PhageAnalysis._process_alignment_file_process_record(data[7:])
                if read_id1 != act_read_id:
                    if 0 < len(act_taxa) <= number_of_hits:
                        for act_taxon in act_taxa:
                            PhageAnalysis._add_count_to_dictionary(taxon2count, act_taxon)
                        for act_genome in act_genomes:
                            PhageAnalysis._add_count_to_dictionary(genome2count, act_genome)
                        act_read_id = read_id1
                        act_taxa = []
                if is_concordant and genome2 == genome1 and abs(coord2 - coord1) < 1500 \
                        and hit_quality1 > score_tsh and hit_quality2 > score_tsh:
                    act_taxa.append(taxon1)
                    act_genomes.append(genome1)
                elif not is_concordant:
                    if genome1 and hit_quality1 > score_tsh:
                        act_taxa.append(taxon1)
                        act_genomes.append(genome1)
                    if genome2 and hit_quality2 > score_tsh:
                        act_taxa.append(taxon2)
                        act_genomes.append(genome2)
        return [genome2count, taxon2count]

    @staticmethod
    def _process_alignment_file(sample_id, filename, output_filename_genome, output_filename_taxon, is_concordant=True, score_tsh=0, number_of_hits=500):
        """
        Generating results to the file!
        :param taxonomy:
        :return:
        """
        print("Processing a file! %s" % filename)
        if os.path.getsize(filename) < 1:
            print("Empty file!")
            return
        is_paired = PhageAnalysis._process_alignment_file_guessing_format(filename)
        if not is_paired:
            [genome2count, taxon2count] = PhageAnalysis._process_alignment_file_raw_single(filename, score_tsh, number_of_hits)
        else:
            [genome2count, taxon2count] = PhageAnalysis._process_alignment_file_raw_paired(filename, score_tsh, number_of_hits, is_concordant)
        with open(output_filename_genome, "w") as fout_genome, open(output_filename_taxon, "w") as fout_taxon:
            for taxon, count in taxon2count.items():
                fout_taxon.write(str(taxon) + "\t" + str(count) + "\t" + sample_id + "\n")
            for genome, count in genome2count.items():
                fout_genome.write(str(genome) + "\t" + str(count) + "\t" + sample_id + "\n")

    @staticmethod
    def _pivot_separate_tables_from_files(filenames):
        """
        Generating a summary table for the counts
        :param filenames:
        :return:
        """
        print("Number of file names: %d" % len(filenames))
        id2sample_id2counts = {}
        for filename in filenames:
            for line in open(filename):
                [organism_id, count, sample_id] = line.strip().split()
                if organism_id in id2sample_id2counts:
                    id2sample_id2counts[organism_id][sample_id] = int(count)
                else:
                    id2sample_id2counts[organism_id] = {sample_id: int(count)}
        return id2sample_id2counts

    @staticmethod
    def _writing_nested_records_into_file(id_type_name, filename, dictionary):
        """
        Writing nested (dictionary in dictionary) records into a file.
        :param id_type_name: type of the index (genome or taxon id)
        :param filename: filename of the output file
        :param dictionary: data in the records
        :return:
        """
        print("Collecting column ids")
        sample_ids = [sample_id for organism_id in dictionary for sample_id in dictionary[organism_id]]
        sample_ids = list(set(sample_ids))
        print("Number of samples: %d" % len(sample_ids))
        header = [id_type_name] + sample_ids
        with open(filename, "w") as out_file:
            out_file.write("\t".join(header) + "\n")
            for organism_id in dictionary:
                record = [str(organism_id)]
                act_sample2count = dictionary[organism_id]
                for sample_id in sample_ids:
                    if sample_id in act_sample2count:
                        record.append(act_sample2count[sample_id])
                    else:
                        record.append(0)
                out_file.write("\t".join([str(attr) for attr in record]) + "\n")

    @staticmethod
    def _generate_parameters_for_aggregations(taxonomy, sample_ids, sample_id2file, alignment_folder, output_folder, base_name_genome, base_name_taxon, score_tsh, number_of_hits,
                                              is_concordant=True):
        """
        Generates the necessary parameters for the function calls (it is necessary for the parallel processing)
        :param taxonomy:  a taxonomy object for calculating the LCA value
        :param sample_ids: sample ids to be processed
        :param sample_id2file:  mapping between the file names the sample ids
        :param alignment_folder: path to the folder, where the alignments are deposited (all file will be processed)
        :param output_folder:  path to the folder where the temporary and the final results are deposited
        :param base_name_genome: base name substring that the temporary results file will begin
        :param base_name_taxon: base name substring that the temporary results file will begin
        :return:
        """
        params = []
        params_lca = []
        for sample_id in sample_ids:
            filename = sample_id2file[sample_id]
            output_filename_genome = join(output_folder, sample_id + base_name_genome)
            output_filename_taxon = join(output_folder, sample_id + base_name_taxon)
            act_params = (sample_id, join(alignment_folder, filename), output_filename_genome, output_filename_taxon, is_concordant, score_tsh, number_of_hits)
            param_lca = (taxonomy, sample_id, join(alignment_folder, filename), output_filename_taxon, is_concordant, score_tsh, number_of_hits)
            params.append(act_params)
            params_lca.append(param_lca)
        return params, params_lca

    def _process_alignments_lca(self, sample_ids, sample_id2file, output_folder, alignment_folder, k_pool, taxonomy=None, score_tsh=0.9, number_of_hits=500):
        """
        Counts the assigned reads to the taxa and genomes using Lowest Common Ancestor principle
        :param sample_ids: sample ids to be processed
        :param sample_id2file: mapping between the filenames the sample ids
        :param output_folder: path to the folder where the temporary and the final results are deposited
        :param alignment_folder: path to the folder, where the alignments are deposited (all file will be processed)
        :param k_pool: number of threads for parallel run
        :param taxonomy: a taxonomy object for calculating the LCA value
        :return:
        """
        from multiprocessing import Pool
        base_name_genome = ".lca_" + self.base_name_genome + "_alignment_tsh." + str(score_tsh) + "_max_hit_per_read" + str(number_of_hits)
        base_name_taxon = ".lca_" + self.base_name_taxon + "_alignment_tsh." + str(score_tsh) + "_max_hit_per_read" + str(number_of_hits)
        print(base_name_taxon)
        print("output folder:", output_folder)
        print(join(output_folder, base_name_taxon))
        output_summary_filename_taxon = join(output_folder, base_name_taxon)
        [params, params_lca] = PhageAnalysis._generate_parameters_for_aggregations(taxonomy, sample_ids, sample_id2file, alignment_folder, output_folder, base_name_genome,
                                                                                   base_name_taxon, score_tsh, number_of_hits)
        pool = Pool(k_pool)
        pool.starmap(PhageAnalysis._process_alignment_file_lca, params_lca)
        output_filenames_taxa = [record[3] for record in params]
        taxa2samples2count = PhageAnalysis._pivot_separate_tables_from_files(output_filenames_taxa)
        print("Start writing!")
        print("OUtput summary file: ", output_summary_filename_taxon[1:])
        PhageAnalysis._writing_nested_records_into_file("taxon_id", output_summary_filename_taxon, taxa2samples2count)
        try:
            # [os.remove(filename) for filename in output_filenames_taxa]
            pass
        except FileNotFoundError:
            print("File not found!")
        return output_summary_filename_taxon

    def _process_alignments_raw(self, sample_ids, sample_id2file, output_folder, alignment_folder, k_pool, score_tsh=0.9, number_of_hits=500):
        """
        Counts the assigned reads to the taxa and genomes.
        :param sample_ids:  sample ids to be processed
        :param sample_id2file: mapping between the file names the sample ids
        :param output_folder:  path to the folder where the temporary and the final results are deposited
        :param alignment_folder: path to the folder, where the alignments are deposited (all file will be processed)
        :param k_pool: number of threads for parallel run
        :return:
        """
        from multiprocessing import Pool

        base_name_genome = self.base_name_genome + "_alignment_tsh." + str(score_tsh) + "_max_hit_per_read" + str(number_of_hits)
        base_name_taxon = self.base_name_taxon + "_alignment_tsh." + str(score_tsh) + "_max_hit_per_read" + str(number_of_hits)
        output_summary_filename_genome = join(output_folder, base_name_genome)
        output_summary_filename_taxon = join(output_folder, base_name_taxon)

        pool = Pool(k_pool)
        [params, _] = PhageAnalysis._generate_parameters_for_aggregations(None, sample_ids, sample_id2file, alignment_folder, output_folder, self.base_name_genome,
                                                                          self.base_name_taxon, score_tsh, number_of_hits)
        if k_pool == 1:
            print("Processing data directly")
            for i in range(len(params)):
                print(len(params[0]))
                print(params[i][0], params[i][1], params[i][2], params[i][3], params[i][4], params[i][5], params[i][6])
                self._process_alignment_file(params[i][0], params[i][1], params[i][2], params[i][3], params[i][4], params[i][5], params[i][6])
        else:
            pool.starmap(PhageAnalysis._process_alignment_file, params)
        output_filenames_taxa = [record[3] for record in params]
        output_filenames_genome = [record[2] for record in params]
        taxa2samples2count = PhageAnalysis._pivot_separate_tables_from_files(output_filenames_taxa)
        genome2samples2count = PhageAnalysis._pivot_separate_tables_from_files(output_filenames_genome)
        PhageAnalysis._writing_nested_records_into_file("taxon_id", output_summary_filename_taxon, taxa2samples2count)
        PhageAnalysis._writing_nested_records_into_file("genome_id", output_summary_filename_genome, genome2samples2count)
        try:
            pass
            # [os.remove(filename) for filename in output_filenames_taxa]
            # [os.remove(filename) for filename in output_filenames_genome]
        except FileNotFoundError:
            print("File not found!")
        return output_summary_filename_genome, output_summary_filename_taxon

    def process_alignments(self, alignment_folder, output_folder, k_pool=10, is_lca=False, taxonomy=None, score_tsh=0.9, number_of_hits=500):
        """
        Extract the information from processed alignment files
        :param alignment_folder: path to the folder, where the alignments are deposited (all file will be processed)
        :param output_folder: path to the folder where the temporary and the final results are deposited
        :param k_pool: number of threads for parallel run
        :param is_lca: True if the calculation of LCA value  is required
        :param taxonomy: a taxonomy object for calculating the LCA value
        :param score_tsh:
        :param number_of_hits:
        :return: path to the output summary file
        """

        sample_ids, sample_id2file = self._get_list_of_analized_samples(alignment_folder)
        if not os.path.exists(output_folder):
            os.mkdir(output_folder)
        if is_lca:
            output_summary_filename_taxon = self._process_alignments_lca(sample_ids, sample_id2file, output_folder, alignment_folder, k_pool, taxonomy, score_tsh, number_of_hits)
            return output_summary_filename_taxon
        else:
            output_summary_filename_genome, output_summary_filename_taxon = self._process_alignments_raw(sample_ids, sample_id2file, output_folder, alignment_folder, k_pool,
                                                                                                         score_tsh, number_of_hits)
            return output_summary_filename_genome, output_summary_filename_taxon
