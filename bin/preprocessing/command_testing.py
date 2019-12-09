from os.path import *

import MetaPhlan2Runner as MetaPhlan2Runner
import PreprocessingConfig as config

datasetnames = config.get_dataset_names()

# ---- For krakren2 ----
"""
db_path = '/gfs/progs/kraken2/Database/minikraken2/minikraken2_v2_8GB_201904_UPDATE'

#db_path = '/gfs/progs/kraken2/Database/minikraken_20171019_8GB'

# For tring:

dataset = 'RampelliS_2015'
input_file_pairs = Kraken2Runner.get_file_pairs(dataset)
kraken_savedir = join('/gfs/data/curated_metagenomes_kraken2', dataset)

Kraken2Runner.run_kraken2_program(kraken_savedir, input_file_pairs, threads=12, db_path=db_path)

#For all:

for dataset in datasetnames:
    input_file_pairs = Kraken2Runner.get_file_pairs(dataset)
    kraken_savedir = join('/gfs/data/curated_metagenomes_kraken2', dataset)

    Kraken2Runner.run_kraken2_program(kraken_savedir, input_file_pairs, threads=12, db_path=db_path)
"""

# ---- For metahplan2 ----

# For try:

"""
dataset = 'RampelliS_2015'
input_file_pairs = MetaPhlan2Runner.get_file_pairs(dataset)
metphlan_savedir = join('/gfs/data/curated_metagenomes_metaphlan2', dataset)

MetaPhlan2Runner.run_metaphlan2_program(metphlan_savedir, input_file_pairs, threads=12)
"""

"""
# For all
for dataset in datasetnames:
    #if dataset not in ["AsnicarF_2017", "ChngKR_2016", "FengQ_2015", "Heitz-BuschartA_2016"]:
    if dataset in ["VogtmannE_2016", "YuJ_2015", "ZellerG_2014"]:
        input_file_pairs = MetaPhlan2Runner.get_file_pairs(dataset)
        metphlan_savedir = join('/gfs/data/curated_metagenomes_metaphlan2', dataset)

        MetaPhlan2Runner.run_metaphlan2_program(metphlan_savedir, input_file_pairs, threads=24)
"""

bowtie_files = MetaPhlan2Runner.get_bowtie_files("ZellerG_2014")
metphlan_savedir = join('/gfs/data/curated_metagenomes_metaphlan2', "ZellerG_2014")

MetaPhlan2Runner.run_metaphlan2_from_bowtie(metphlan_savedir, bowtie_files, threads=24)