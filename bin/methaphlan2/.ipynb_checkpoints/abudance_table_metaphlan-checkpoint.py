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
PASSWORD = 'LrJPRUS73r84'

class Connection():

    def __init__(self):
        dsn_tns = cx_Oracle.makedsn(SERVER, PORT, service_name=SERVICE)
        self.connection = cx_Oracle.connect(USER, PASSWORD, dsn_tns)

    def get_connection(self):
        return self.connection

    def select(self, sql):
        return pd.read_sql(sql, con=self.connection)

con = Connection()

path = '//gfs/data/curated_metagenomes_metaphlan2/'

files = {}
# r=root, d=directories, f = files
for r, d, f in os.walk(path):
    for directory in d:
        for r1, d1, f1 in os.walk(os.path.join(path,directory)):
            for file in f1:
                if '.profile.txt' in file:
                    if directory in files:
                        files[directory].append(os.path.join(r1, file))
                    else:
                        files[directory] = [(os.path.join(r1, file))]

def get_data_for_sample(path_to_file):
    dataset = path_to_file.split('/')[-2]
    sample_id = path_to_file.split('/')[-1].split('.')[0]

    df = pd.read_csv(path_to_file, delimiter = '\t')

    profile_df_2 = df
    profile_df_2['Rank'] = None
    profile_df_2['K'] = None
    profile_df_2['P'] = None
    profile_df_2['C'] = None
    profile_df_2['O'] = None
    profile_df_2['F'] = None
    profile_df_2['G'] = None
    profile_df_2['S'] = None
    profile_df_2['T'] = None
    profile_df_2['Abundance'] = profile_df_2["Metaphlan2_Analysis"]

    for p in profile_df_2["#SampleID"]:
        row_level = None
        tags = p.split('|')
        for tag in tags:
            tmp = tag.split('__')
            try:
                level = tmp[0]
                name = tmp[1]
                profile_df_2
                row_level = level

                profile_df_2.loc[profile_df_2['#SampleID'] == p, level] = name
            except:
                print("Can not split:", tmp)
        profile_df_2.loc[profile_df_2['#SampleID'] == p, 'Rank'] = row_level

    def get_name_for_row(row):
        try:
            return row[row.Rank]
        except:
            return '???'

    profile_df_2["Name"] = profile_df_2.apply(get_name_for_row, axis=1)

    selected = profile_df_2[["Rank", "Name", "Abundance"]]

    selected['Sample_ID'] = sample_id
    selected['Dataset'] = dataset


    return selected[["Dataset", "Sample_ID", "Rank", "Name", "Abundance"]]

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

    df["DATASET_NAME"] = df["Dataset"]
    merged = df.merge(df_metadata, on='DATASET_NAME')

    col_to_idx = {}
    for i in range(len(merged.columns)):
        col_to_idx[merged.columns[i]] = i

    joined = merged.loc[[row[1][col_to_idx["Sample_ID"]] in row[1][col_to_idx["NCBI_ACCESSION"]]
                        for row in merged.iterrows()]]

    output = joined[["Sample_ID", "Dataset", "Rank", "Name", "Abundance",
                     "STUDY_CONDITION", "BODY_SITE", "SUBJECTID",
                     "COUNTRY", "DISEASE", "AGE", "AGE_CATEGORY"]]
    output = output.reset_index()
    output = output.drop(columns="index")
    return output

dataset_dfs = {}
i = 0;
for dataset_name, filelist in files.items():
    if i == 0:
        df_list = []
        for file in filelist:
            print(file)
            df = get_data_for_sample(file)
            if df is not None:
                df_list.append(get_data_for_sample(file))
        df_dataset = pd.concat(df_list)
        output = join_metadata(df_dataset, dataset_name)

        path_to_output = os.path.join(path, dataset_name + ".csv")
        with open(path_to_output, "w") as f:
            f.write(output.to_csv())
    i = i+1
    print(dataset_name, " kész.")

