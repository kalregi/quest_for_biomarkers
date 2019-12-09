#!/usr/bin/env python
# -*-coding: utf8 -*-

import os
from os.path import join, isfile
import yaml



def yaml_parser_preprocess_pipeline_input_parameter_get_description_df(config):
    """
    Creates a very basic description table for the proteins in the table
    :param config:  Pipeline info
    :return:
    """
    df= pandas.DataFrame.from_dict(config['query_system']['profile_description'], orient = 'index')
    df['query_protein_name'] = df.index
    return df


def yaml_parser_preprocess_pipeline_input_parameter(config):
    """
    Preprocessing the pipeline parameters (extending the prefixes and stuff)
    :param config:  YAML configuration
    :return:
    """
    description_df = yaml_parser_preprocess_pipeline_input_parameter_get_description_df(config)
    config['query_system']['description'] = description_df


def yaml_parser_preprocess_pipeline_outpout_parameter(config):
    """
    Preprocessing the pipeline parameters (extending the prefixes and stuff)
    :param config:  YAML configuration
    :return:
    """
    qs_system_name = config['query_system']['system_short_name']
    seq_db_name = config['refseq']['refseq_data']
    save_dir = config['output']['save_dir']

    config['output']['system_save_dir'] = join(save_dir, qs_system_name)
    config['output']['hmm_hits_raw_output_dir'] = join(save_dir, qs_system_name, config['output']['hmm_hits_raw_output_prefix_dir'])
    config['output']['hmm_hits_dir'] = join(save_dir, qs_system_name, config['output']['hmm_hits_prefix_dir'])
    config['output']['filtered_hmm_hits_file'] = join(save_dir,qs_system_name, seq_db_name + '.' + qs_system_name + '.' + config['output']['filtered_hmm_hits_postfix'])
    config['output']['system_specific_annotations_dir'] = join(save_dir, qs_system_name, config['output']['system_specific_annotations_dirname_postfix'])
    config['output']['system_neighborhood_specific_annotations_dir'] = join(save_dir, qs_system_name, config['output']['system_neighborhood_specific_annotations_dirname_postfix'])

    config['output']['systems_file'] = join(save_dir, qs_system_name, seq_db_name + '.' + qs_system_name + '.' + config['output']['act_qs_system_patterns_summary_postfix'])
    config['output']['systems_components_file'] = join(save_dir,qs_system_name, seq_db_name + '.' + qs_system_name + '.' + config['output']['act_qs_system_patterns_postfix'])

    #database:

    config['database']['password'] = ''.join([line.strip() for line in open(config['database']['password_file'])])
