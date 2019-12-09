
from os.path import *
from Preprocessing import Preprocessing

dir_name = "/gfs/data/curated_metagenomes/KarlssonFH_2013"
processor = Preprocessing(join(dirname(__file__), dir_name))
#processor.run_original_qc_report()

#processor.run_trimm_and_qc_report()
processor.list_files(trimmed=True)
processor.generate_quality_control_profile(trimmed=True)