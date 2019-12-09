

import os
from os.path import *
from os import listdir
from TrimmomaticRunner import TrimmomaticRunner


class Preprocessing:

    # Python script írása, mely a következő műveleteket végzi el:
    # 1. Fájlok listázása: Az adott szekcióhoz kilistázza a fájlokat
    # 2. QC profil: Lefuttatja rajtuk a fastqc-t és a multiqc-t
    # 3. Paraméterek meghatározása: Trimmomatic-hoz meghatározza a paramétereket
    # 4. Trimmomatic: Lefuttatja a trimmomaticot a trimmomatic modul segítségével

    def __init__(self, directory):
        self.directory = directory
        self.file_list = []
        self.trimmomatic_parameters = {}

    def list_files(self):
        fastq_files = [f for f in listdir(self.directory) if isfile(join(self.directory, f)) and self.check_if_file_processable(f)]
        print(fastq_files)

    def check_if_file_processable(self, filename):
        check = filename.split(".")
        print(splitext(filename))
        return (len(check) == 2 and filename[1] == ".fastq") or \
               (len(check) == 3 and filename[1] == ".fastq" and filename[2] == ".gz")

    def generate_quality_control_profile(self):
        pass

    def get_parameters_for_trimmomatic(self):
        # TODO: find appropriate parameters for Trimmomatic
        self.trimmomatic_parameters = {
            "illuminaclip" : ("TruSeq3-PE.fa", 2, 30, 10),
            "leading": 3,
            "trailing":3,
            "slidingwindow" :(6, 20),
            "minlen":36
        }

    def run_trimmomatic(self):
        runner = TrimmomaticRunner("/gfs/progs/Trimmomatic-0.36/trimmomatic-0.36.jar")
        for files in self.file_list:
            path_to_files = [files[0], files[1]]
            runner.preprocess_fastq_files_trimmomatic_pair_end(path_to_fastq_files=path_to_files,
                                                               illuminaclip=self.trimmomatic_parameters["illuminaclip"],
                                                               leading=self.trimmomatic_parameters["leading"],
                                                               trailing=self.trimmomatic_parameters["leading"],
                                                               slidingwindow=self.trimmomatic_parameters["slidingwindow"],
                                                               minlen=self.trimmomatic_parameters["minlen"]
                                                               )
