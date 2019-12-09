import unittest
import PreprocessingConfig as cf

class PreprocessingConfigTest(unittest.TestCase):

    def test_get_source_dir(self):
        self.assertEqual(cf.get_source_dir("KarlssonFH_2013"),
                         "/gfs/data/curated_metagenomes/KarlssonFH_2013")

    def test_get_qcreport_fastq_thread_count(self):
        self.assertEqual(cf.get_qcreport_fastq_thread_count(),22)

    def test_get_qcreport_source_path(self):
        self.assertEqual(cf.get_qcreport_source_path("KarlssonFH_2013"),
                         "/gfs/data/curated_metagenomes/KarlssonFH_2013")

        self.assertEqual(cf.get_qcreport_source_path("KarlssonFH_2013", True),
                         "/gfs/data/curated_metagenomes_preprocessed/KarlssonFH_2013/trimmed_files")

    def test_get_qcreport_destination_path(self):
        self.assertEqual(cf.get_qcreport_destination_path("KarlssonFH_2013"),
                         "/gfs/data/curated_metagenomes_qc/raw_qc/KarlssonFH_2013")
        self.assertEqual(cf.get_qcreport_destination_path("KarlssonFH_2013", True),
                         "/gfs/data/curated_metagenomes_qc/trimmed_qc/KarlssonFH_2013")

    def test_get_multiqc_destination_path(self):
        self.assertEqual(cf.get_qcreport_destination_path("KarlssonFH_2013"),
                         "/gfs/data/curated_metagenomes_qc/raw_qc/KarlssonFH_2013")
        self.assertEqual(cf.get_qcreport_destination_path("KarlssonFH_2013", True),
                         "/gfs/data/curated_metagenomes_qc/trimmed_qc/KarlssonFH_2013")


    def test_get_multiqc_output_destination_file_name(self):
        self.assertEqual(cf.get_multiqc_output_destination_file_name("KarlssonFH_2013"),
                         "/gfs/data/curated_metagenomes_qc/KarlssonFH_2013_raw_qc.index.html")
        self.assertEqual(cf.get_multiqc_output_destination_file_name("KarlssonFH_2013", True),
                         "/gfs/data/curated_metagenomes_qc/KarlssonFH_2013_trimmed_qc.index.html")

    def test_get_trimommatic_parameters(self):
        self.assertEqual(cf.get_trimommatic_parameters("KarlssonFH_2013"),
                         {
                            "gzip_streaming": True,
                            "illuminaclip" : ("TruSeq3-PE.fa", 2, 30, 10),
                            "leading": 20,
                            "trailing": 20,
                            "slidingwindow": (6, 20),
                            "minlen": 36,
                             "pair": True,
                             "process": False
                         })

    def test_get_trimmomatic_path(self):
        self.assertEqual(cf.get_trimmomatic_path(),
                         "/gfs/progs/Trimmomatic-0.36/trimmomatic-0.36.jar")

    def test_get_trimmomatic_destination_path(self):
         self.assertEqual(cf.get_trimmomatic_destination_path("KarlssonFH_2013"),
                         "/gfs/data/curated_metagenomes_preprocessed/KarlssonFH_2013/trimmed_files")

    def test_get_trimmomatic_thread_count(self):
        self.assertEqual(cf.get_trimmomatic_thread_count(),22)

    def test_get_trimmomatic_output_extensions(self):
        self.assertEqual(cf.get_trimmomatic_output_extensions(),
                         {
                             "pair_1_ext": ".trimmed",
                             "junk_pair_1_ext": ".junk.trimmed",
                             "pair_2_ext": ".trimmed",
                             "junk_pair_2_ext": ".junk.trimmed",
                             "single_ext": ".trimmed",
                             "single_junk_ext": ".junk.trimmed"
                         })

    def test_is_trimmable(self):
        self.assertFalse(cf.is_trimmable("Castro-NallarE_2015"))
        self.assertTrue(cf.is_trimmable("AsnicarF_2017"))

    def test_is_pair(self):
        self.assertFalse(cf.is_pair("Castro-NallarE_2015"))
        self.assertTrue(cf.is_pair("AsnicarF_2017"))

if __name__ == '__main__':
    unittest.main()