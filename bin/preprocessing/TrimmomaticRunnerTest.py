import unittest
import TrimmomaticRunner as TrimmomaticRunner

class TrimmomaticRunnerTest(unittest.TestCase):

    def test_get_trimmomatic_pair_end_command(self):
        self.assertEqual(TrimmomaticRunner._get_trimmomatic_pair_end_command(
                                path_to_fastq_files=("gfs/data/curated_metagenomes/KarlssonFH_2013/ERR260132_1.fastq.gz",
                                                     "gfs/data/curated_metagenomes/KarlssonFH_2013/ERR260132_2.fastq.gz"),
                                output_dir="gfs/data/curated_metagenomes_trimmed/KarlssonFH_2013",
                                n_process=3,
                                illuminaclip=("TruSeq3-PE.fa", 2, 30, 10),
                                leading=3,
                                trailing=3,
                                slidingwindow=(6,20),
                                minlen=36
                                ),
                            "java -jar /gfs/progs/Trimmomatic-0.36/trimmomatic-0.36.jar PE -threads 3 -phred33 gfs/data/curated_metagenomes/KarlssonFH_2013/ERR260132_1.fastq.gz gfs/data/curated_metagenomes/KarlssonFH_2013/ERR260132_2.fastq.gz gfs/data/curated_metagenomes_trimmed/KarlssonFH_2013/ERR260132_1.fastq.trimmed.gz gfs/data/curated_metagenomes_trimmed/KarlssonFH_2013/ERR260132_1.fastq.junk.trimmed.gz gfs/data/curated_metagenomes_trimmed/KarlssonFH_2013/ERR260132_2.fastq.trimmed.gz gfs/data/curated_metagenomes_trimmed/KarlssonFH_2013/ERR260132_2.fastq.junk.trimmed.gz ILLUMINACLIP:/gfs/progs/Trimmomatic-0.36/adapters/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:6:20 MINLEN:36 "
        )

    def test_get_parameters(self):
        self.assertEqual(TrimmomaticRunner._get_parameters(
            illuminaclip=("TruSeq3-PE.fa", 2, 30, 10),
            leading=3,
            trailing=3,
            slidingwindow=(6,20),
            minlen=36),
            "ILLUMINACLIP:/gfs/progs/Trimmomatic-0.36/adapters/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:6:20 MINLEN:36 ")

    def test_get_output_files(self):
         self.assertEqual(TrimmomaticRunner._get_output_files(
            ["gfs/data/curated_metagenomes/KarlssonFH_2013/ERR260132_1.fastq.gz",
             "gfs/data/curated_metagenomes/KarlssonFH_2013/ERR260132_2.fastq.gz"],
            "gfs/data/curated_metagenomes_trimmed/KarlssonFH_2013/"),
            ("gfs/data/curated_metagenomes_trimmed/KarlssonFH_2013/ERR260132_1.fastq.trimmed.gz",
             "gfs/data/curated_metagenomes_trimmed/KarlssonFH_2013/ERR260132_1.fastq.junk.trimmed.gz",
             "gfs/data/curated_metagenomes_trimmed/KarlssonFH_2013/ERR260132_2.fastq.trimmed.gz",
             "gfs/data/curated_metagenomes_trimmed/KarlssonFH_2013/ERR260132_2.fastq.junk.trimmed.gz"))

if __name__ == '__main__':
    unittest.main()