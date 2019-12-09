import pandas as pd
import os
import seaborn as sns
from skbio.diversity import alpha_diversity
from skbio.diversity import beta_diversity
from skbio import TreeNode
from io import StringIO
from skbio.stats.ordination import pcoa

import cx_Oracle

SERVER = 'dboracle.itk.ppke.hu'
PORT = 1521
SERVICE = 'rsc.itk.ppke.hu'
USER = 'metagenome'
PASSWORD = '.'

class Connection():

    def __init__(self):
        dsn_tns = cx_Oracle.makedsn(SERVER, PORT, service_name=SERVICE)
        self.connection = cx_Oracle.connect(USER, PASSWORD, dsn_tns)

    def get_connection(self):
        return self.connection

    def select(self, sql):
        return pd.read_sql(sql, con=self.connection)

def get_data_for_sample(path_to_file):
    dataset = path_to_file.split('/')[-2]
    df = pd.read_csv(path_to_file, delimiter = '\t',
                     names=["Abundance", "Root_count", "Taxon_count", "Rank", "NCBI_ID", "Name"])
    sample_id = path_to_file.split('/')[-1].split('.')[0]
    df['Sample_ID'] = sample_id
    df['Dataset'] = dataset
    not_found_files = set()
    return df

def join_metadata(df, dataset_name):

    SQL_METADATA = """
        select
            DATASET_NAME,
            STUDY_CONDITION,
            BODY_SITE,
            SUBJECTID,
            NCBI_ACCESSION,
            COUNTRY,
            DISEASE,
            GENDER,
            AGE,
            AGE_CATEGORY
        from METADATA
        where dataset_name = '{dataset_name}'
        """.replace("{dataset_name}", dataset_name)
    df_metadata = con.select(SQL_METADATA)
    df_metadata = df_metadata.loc[[x is not None for x in df_metadata.NCBI_ACCESSION]]

    tmp = pd.DataFrame(df_metadata.NCBI_ACCESSION.str.split(';').tolist(),
                                    index = [df_metadata.DATASET_NAME,
                                             df_metadata.STUDY_CONDITION, df_metadata.BODY_SITE,
                                            df_metadata.SUBJECTID, df_metadata.COUNTRY,
                                            df_metadata.DISEASE, df_metadata.GENDER,
                                            df_metadata.AGE, df_metadata.AGE_CATEGORY])
    df_metadata_distinct = tmp.stack()

    new_metadata_df = df_metadata_distinct.reset_index([0, 'DATASET_NAME', 'STUDY_CONDITION',
        'BODY_SITE', 'SUBJECTID', 'COUNTRY', 'DISEASE', 'GENDER', 'AGE', 'AGE_CATEGORY'])
    new_metadata_df['NCBI_ACCESSION'] = new_metadata_df[0]
    new_metadata_df = new_metadata_df.drop(columns=0)

    df["DATASET_NAME"] = df["Dataset"]
    df["NCBI_ACCESSION"] = df["Sample_ID"]

    joined = df.merge(new_metadata_df, on='NCBI_ACCESSION')

    """
    col_to_idx = {}
    for i in range(len(merged.columns)):
        col_to_idx[merged.columns[i]] = i

    joined = merged.loc[[row[1][col_to_idx["Sample_ID"]] in row[1][col_to_idx["NCBI_ACCESSION"]]
                        for row in merged.iterrows()]]
    """

    output = joined[["Sample_ID", "Dataset", "Rank", "Name", "Abundance", "STUDY_CONDITION", "BODY_SITE", "SUBJECTID",
                 "COUNTRY", "DISEASE", "AGE", "AGE_CATEGORY"]]
    output = output.reset_index()
    output = output.drop(columns="index")
    return output

con = Connection()

path = '//gfs/data/curated_metagenomes_kraken2/'


files = {}
# r=root, d=directories, f = files
for r, d, f in os.walk(path):
    for directory in d:
        for r1, d1, f1 in os.walk(os.path.join(path,directory)):
            for file in f1:
                if '.report.out' in file:
                    if directory in files:
                        files[directory].append(os.path.join(r1, file))
                    else:
                        files[directory] = [(os.path.join(r1, file))]

dataset_dfs = {}
for dataset_name, filelist in files.items():
    df_list = []
    for file in filelist:
        df = get_data_for_sample(file)
        if df is not None:
            df_list.append(get_data_for_sample(file))
    df_dataset = pd.concat(df_list)

    joined = join_metadata(df_dataset, dataset_name)

    path_to_output = os.path.join(path, dataset_name + ".csv")
    with open(path_to_output, "w") as f:
        f.write(joined.to_csv())
    print(dataset_name, "done.")