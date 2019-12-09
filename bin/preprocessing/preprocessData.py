
import Preprocessing as Preprocessing
import PreprocessingConfig as config


#datasetList = ["AsnicarF_2017", "BritoIL_2016", "Castro-NallarE_2015", "ChngKR_2016", "FengQ_2015",
#               "Heitz-BuschartA_2016", "LeChatelierE_2013", "LiuW_2016", "LomanNJ_2013", "HMP_2012", "KarlssonFH_2013", "NielsenHB_2014",
#               "Obregon-TitoAJ_2015", "OhJ_2014", "OhJ_2016", "QinJ_2012", "QinN_2014", "RampelliS_2015",
#               "RaymondF_2016", "SchirmerM_2016", "TettAJ_2016", "VatanenT_2016", "VincentC_2016", "VogtmannE_2016",
#               "XieH_2016", "YuJ_2015", "ZellerG_2014"]

#for dataset in datasetList:
#    Preprocessing.generate_quality_control_profile(dataset, threads=10)

#Preprocessing.multiqc("KarlssonFH_2013")

#print(datetime.datetime.fromtimestamp(time.time()).strftime('%Y-%m-%d_%H:%M:%S'))

#file_list_path = "/gfs/data/curated_metagenomes_qc/FileLists"
#lister.print_file_list("AsnicarF_2017", "AsnicarF_2017"+"_filelist.csv", output_file_path=file_list_path)

#for dataset in datasetList:
#    lister.print_file_list(dataset, dataset+"_filelist.csv", output_file_path=file_list_path)

# Trimmomatic

#Preprocessing.run_trimmomatic("AsnicarF_2017")

#dataset_names = config.get_dataset_names()
#dataset_names = ["BritoIL_2016", "Castro-NallarE_2015", "ChngKR_2016", "FengQ_2015",
#               "Heitz-BuschartA_2016", "LeChatelierE_2013", "LiuW_2016", "LomanNJ_2013", "HMP_2012", "KarlssonFH_2013", "NielsenHB_2014",
#               "Obregon-TitoAJ_2015", "OhJ_2014", "OhJ_2016", "QinJ_2012", "QinN_2014", "RampelliS_2015",
#               "RaymondF_2016", "SchirmerM_2016", "TettAJ_2016", "VatanenT_2016", "VincentC_2016", "VogtmannE_2016",
#               "XieH_2016", "YuJ_2015", "ZellerG_2014"]
#print(dataset_names)
#for dataset_name in dataset_names:
#    Preprocessing.run_trimmomatic(dataset_name)

#for dataset_name in dataset_names:
#    Preprocessing.generate_quality_control_profile(dataset_name=dataset_name, trimmed = True, threads=10)

missed_dataset_names = ["OhJ_2014"]

#["BritoIL_2016", "NielsenHB_2014", 'OhJ_2014', 'RaymondF_2016']

for dataset in missed_dataset_names:
    Preprocessing.generate_quality_control_profile(dataset, threads=2)