
from TrimmomaticRunner import TrimmomaticRunner
from os.path import *

runner = TrimmomaticRunner("/gfs/progs/Trimmomatic-0.36/trimmomatic-0.36.jar")

folder = join(dirname(__file__),"exampleFiles")

path_to_files = [join(folder,"EXAMPLE_good_1.fastq"), join(folder,"EXAMPLE_good_2.fastq")]

runner.preprocess_fastq_files_trimmomatic_pair_end(path_to_fastq_files=path_to_files)
runner.preprocess_fastq_files_trimmomatic_pair_end(path_to_fastq_files=path_to_files, leading=40)