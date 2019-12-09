import unittest
import FastqFileLister as FastqFileLister

class FastqFileListerTest(unittest.TestCase):


    def test_check_if_file_single_end(self):
        self.assertFalse(FastqFileLister.check_if_file_single_end("ERR321066_2.fastq.gz"))
        self.assertFalse(FastqFileLister.check_if_file_single_end("ERR321089_1.fastq.gz"))
        self.assertTrue(FastqFileLister.check_if_file_single_end("SRR3313035.fastq.gz"))

        self.assertFalse(FastqFileLister.check_if_file_single_end("ERR321066_2.fastq.gz"))
        self.assertFalse(FastqFileLister.check_if_file_single_end("ERR321089_1.fastq.gz"))
        self.assertTrue(FastqFileLister.check_if_file_single_end("SRR3313035.fastq.gz"))

    def test_check_if_file_pair_end(self):
        self.assertTrue(FastqFileLister.check_if_file_pair_end("ERR321066_2.fastq.gz"))
        self.assertTrue(FastqFileLister.check_if_file_pair_end("ERR321089_1.fastq.gz"))
        self.assertFalse(FastqFileLister.check_if_file_pair_end("SRR3313035.fastq.gz"))

    def test_check_if_file_processable(self):
        self.assertTrue(FastqFileLister.check_if_file_processable("ERR321066_2.fastq.gz"))
        self.assertFalse(FastqFileLister.check_if_file_processable("sample.gz"))

    def test_is_second_pair(self):
        self.assertTrue(FastqFileLister.is_second_pair("ERR321066_2.fastq.gz"))
        self.assertFalse(FastqFileLister.is_second_pair("ERR321066_1.fastq.gz"))
        self.assertFalse(FastqFileLister.is_second_pair("ERR321066.fastq.gz"))

    def test_is_first_pair(self):
        self.assertFalse(FastqFileLister.is_first_pair("ERR321066_2.fastq.gz"))
        self.assertTrue(FastqFileLister.is_first_pair("ERR321066_1.fastq.gz"))
        self.assertFalse(FastqFileLister.is_first_pair("ERR321066.fastq.gz"))

    def test_get_key(self):
        self.assertEqual(FastqFileLister.get_key("ERR321066_2.fastq.gz"), "ERR321066")
        self.assertEqual(FastqFileLister.get_key("ERR321066_1.fastq.gz"), "ERR321066")
        self.assertEqual(FastqFileLister.get_key("ERR321066.fastq.gz"), "ERR321066")

    def test_order_files_to_dict(self):
        self.maxDiff = None
        self.assertEqual(FastqFileLister.order_files_to_dict(dataset_name="AsnicarF_2017",
                                                             file_list=["ALMA_1.fastq.gz", "ALMA_2.fastq.gz", "KORTE_1.fastq.gz", "KORTE_2.fastq.gz", "MALAC.fastq.gz"]),
                         {
                             "ALMA": {
                                 "key": "ALMA",
                                 "dataset": "AsnicarF_2017",
                                 "pair_1": "ALMA_1.fastq.gz",
                                 "pair_2": "ALMA_2.fastq.gz",
                                 "single": "",
                                 "is_pair_end": True,
                                 "path": "/gfs/data/curated_metagenomes/AsnicarF_2017"
                             },
                             "KORTE": {
                                 "key": "KORTE",
                                 "dataset": "AsnicarF_2017",
                                 "pair_1": "KORTE_1.fastq.gz",
                                 "pair_2": "KORTE_2.fastq.gz",
                                 "single": "",
                                 "is_pair_end": True,
                                 "path": "/gfs/data/curated_metagenomes/AsnicarF_2017"
                             },
                             "MALAC": {
                                 "key": "MALAC",
                                 "dataset": "AsnicarF_2017",
                                 "pair_1": "",
                                 "pair_2": "",
                                 "single": "MALAC.fastq.gz",
                                 "is_pair_end": False,
                                 "path": "/gfs/data/curated_metagenomes/AsnicarF_2017"
                             }
                         })

    def test_list_files(self):
        self.assertEqual(FastqFileLister.list_files("KarlssonFH_2013"),
                         ['ERR260176_2.fastq.gz', 'ERR260174_2.fastq.gz', 'ERR260175_1.fastq.gz', 'ERR260176_1.fastq.gz', 'ERR260177_1.fastq.gz', 'ERR260132_1.fastq.gz', 'ERR260132_2.fastq.gz', 'ERR260133_1.fastq.gz', 'ERR260133_2.fastq.gz', 'ERR260134_1.fastq.gz', 'ERR260134_2.fastq.gz', 'ERR260135_1.fastq.gz', 'ERR260135_2.fastq.gz', 'ERR260136_1.fastq.gz', 'ERR260136_2.fastq.gz', 'ERR260137_1.fastq.gz', 'ERR260137_2.fastq.gz', 'ERR260138_1.fastq.gz', 'ERR260138_2.fastq.gz', 'ERR260139_1.fastq.gz', 'ERR260139_2.fastq.gz', 'ERR260140_1.fastq.gz', 'ERR260140_2.fastq.gz', 'ERR260141_1.fastq.gz', 'ERR260141_2.fastq.gz', 'ERR260142_1.fastq.gz', 'ERR260142_2.fastq.gz', 'ERR260143_1.fastq.gz', 'ERR260143_2.fastq.gz', 'ERR260144_1.fastq.gz', 'ERR260144_2.fastq.gz','ERR260145_1.fastq.gz', 'ERR260145_2.fastq.gz', 'ERR260146_1.fastq.gz', 'ERR260146_2.fastq.gz', 'ERR260147_1.fastq.gz', 'ERR260147_2.fastq.gz', 'ERR260148_1.fastq.gz', 'ERR260148_2.fastq.gz', 'ERR260149_1.fastq.gz', 'ERR260149_2.fastq.gz', 'ERR260150_1.fastq.gz', 'ERR260150_2.fastq.gz', 'ERR260151_1.fastq.gz', 'ERR260151_2.fastq.gz', 'ERR260152_1.fastq.gz', 'ERR260152_2.fastq.gz', 'ERR260153_1.fastq.gz', 'ERR260153_2.fastq.gz', 'ERR260154_1.fastq.gz', 'ERR260154_2.fastq.gz', 'ERR260155_1.fastq.gz', 'ERR260155_2.fastq.gz', 'ERR260156_1.fastq.gz', 'ERR260156_2.fastq.gz', 'ERR260157_1.fastq.gz', 'ERR260157_2.fastq.gz', 'ERR260158_1.fastq.gz', 'ERR260158_2.fastq.gz', 'ERR260159_1.fastq.gz', 'ERR260159_2.fastq.gz', 'ERR260160_1.fastq.gz', 'ERR260160_2.fastq.gz', 'ERR260161_1.fastq.gz', 'ERR260161_2.fastq.gz', 'ERR260162_1.fastq.gz', 'ERR260162_2.fastq.gz', 'ERR260163_1.fastq.gz', 'ERR260163_2.fastq.gz', 'ERR260164_1.fastq.gz', 'ERR260164_2.fastq.gz', 'ERR260165_1.fastq.gz', 'ERR260165_2.fastq.gz', 'ERR260166_1.fastq.gz', 'ERR260166_2.fastq.gz', 'ERR260167_1.fastq.gz', 'ERR260167_2.fastq.gz', 'ERR260168_1.fastq.gz', 'ERR260168_2.fastq.gz', 'ERR260169_1.fastq.gz', 'ERR260169_2.fastq.gz', 'ERR260170_1.fastq.gz', 'ERR260170_2.fastq.gz', 'ERR260171_1.fastq.gz', 'ERR260171_2.fastq.gz', 'ERR260172_1.fastq.gz', 'ERR260172_2.fastq.gz', 'ERR260173_1.fastq.gz', 'ERR260173_2.fastq.gz', 'ERR260174_1.fastq.gz', 'ERR260175_2.fastq.gz', 'ERR260177_2.fastq.gz', 'ERR260178_1.fastq.gz', 'ERR260178_2.fastq.gz', 'ERR260179_1.fastq.gz', 'ERR260179_2.fastq.gz', 'ERR260180_1.fastq.gz', 'ERR260180_2.fastq.gz', 'ERR260181_1.fastq.gz', 'ERR260181_2.fastq.gz', 'ERR260182_1.fastq.gz', 'ERR260182_2.fastq.gz', 'ERR260183_1.fastq.gz', 'ERR260183_2.fastq.gz', 'ERR260184_1.fastq.gz', 'ERR260184_2.fastq.gz', 'ERR260185_1.fastq.gz', 'ERR260185_2.fastq.gz', 'ERR260186_1.fastq.gz', 'ERR260186_2.fastq.gz', 'ERR260187_1.fastq.gz', 'ERR260187_2.fastq.gz', 'ERR260188_1.fastq.gz', 'ERR260188_2.fastq.gz', 'ERR260189_1.fastq.gz', 'ERR260189_2.fastq.gz', 'ERR260190_1.fastq.gz', 'ERR260190_2.fastq.gz', 'ERR260191_1.fastq.gz', 'ERR260191_2.fastq.gz', 'ERR260192_1.fastq.gz', 'ERR260192_2.fastq.gz', 'ERR260193_1.fastq.gz', 'ERR260193_2.fastq.gz', 'ERR260194_1.fastq.gz', 'ERR260194_2.fastq.gz', 'ERR260195_1.fastq.gz', 'ERR260195_2.fastq.gz', 'ERR260196_1.fastq.gz', 'ERR260196_2.fastq.gz', 'ERR260197_1.fastq.gz', 'ERR260197_2.fastq.gz', 'ERR260198_1.fastq.gz', 'ERR260198_2.fastq.gz', 'ERR260199_1.fastq.gz', 'ERR260199_2.fastq.gz', 'ERR260200_1.fastq.gz', 'ERR260200_2.fastq.gz', 'ERR260201_1.fastq.gz', 'ERR260201_2.fastq.gz', 'ERR260202_1.fastq.gz', 'ERR260202_2.fastq.gz', 'ERR260203_1.fastq.gz', 'ERR260203_2.fastq.gz', 'ERR260204_1.fastq.gz', 'ERR260204_2.fastq.gz', 'ERR260205_1.fastq.gz', 'ERR260205_2.fastq.gz', 'ERR260206_1.fastq.gz', 'ERR260206_2.fastq.gz', 'ERR260207_1.fastq.gz', 'ERR260207_2.fastq.gz', 'ERR260208_1.fastq.gz', 'ERR260208_2.fastq.gz', 'ERR260209_1.fastq.gz', 'ERR260209_2.fastq.gz', 'ERR260210_1.fastq.gz', 'ERR260210_2.fastq.gz', 'ERR260211_1.fastq.gz', 'ERR260211_2.fastq.gz', 'ERR260212.fastq.gz', 'ERR260213.fastq.gz', 'ERR260214_1.fastq.gz', 'ERR260214_2.fastq.gz', 'ERR260215_1.fastq.gz', 'ERR260215_2.fastq.gz', 'ERR260216_1.fastq.gz', 'ERR260216_2.fastq.gz', 'ERR260217_1.fastq.gz', 'ERR260217_2.fastq.gz', 'ERR260218_1.fastq.gz', 'ERR260218_2.fastq.gz', 'ERR260219_1.fastq.gz', 'ERR260219_2.fastq.gz', 'ERR260220_1.fastq.gz', 'ERR260220_2.fastq.gz', 'ERR260221_1.fastq.gz', 'ERR260221_2.fastq.gz', 'ERR260222_1.fastq.gz', 'ERR260222_2.fastq.gz', 'ERR260223_1.fastq.gz', 'ERR260223_2.fastq.gz', 'ERR260224_1.fastq.gz', 'ERR260224_2.fastq.gz', 'ERR260225_1.fastq.gz', 'ERR260225_2.fastq.gz', 'ERR260226_1.fastq.gz', 'ERR260226_2.fastq.gz', 'ERR260227_1.fastq.gz', 'ERR260227_2.fastq.gz', 'ERR260228_1.fastq.gz', 'ERR260228_2.fastq.gz', 'ERR260229_1.fastq.gz', 'ERR260229_2.fastq.gz', 'ERR260230_1.fastq.gz', 'ERR260230_2.fastq.gz', 'ERR260231_1.fastq.gz', 'ERR260231_2.fastq.gz', 'ERR260232_1.fastq.gz', 'ERR260232_2.fastq.gz', 'ERR260233_1.fastq.gz', 'ERR260233_2.fastq.gz', 'ERR260234_1.fastq.gz', 'ERR260234_2.fastq.gz', 'ERR260235_1.fastq.gz','ERR260235_2.fastq.gz', 'ERR260236_1.fastq.gz', 'ERR260236_2.fastq.gz', 'ERR260237_1.fastq.gz', 'ERR260237_2.fastq.gz', 'ERR260238_1.fastq.gz', 'ERR260238_2.fastq.gz', 'ERR260239_1.fastq.gz', 'ERR260239_2.fastq.gz', 'ERR260240_1.fastq.gz', 'ERR260240_2.fastq.gz', 'ERR260241_1.fastq.gz', 'ERR260241_2.fastq.gz', 'ERR260242_1.fastq.gz', 'ERR260242_2.fastq.gz', 'ERR260243_1.fastq.gz', 'ERR260243_2.fastq.gz', 'ERR260244_1.fastq.gz', 'ERR260244_2.fastq.gz', 'ERR260245_1.fastq.gz', 'ERR260245_2.fastq.gz', 'ERR260246_1.fastq.gz', 'ERR260246_2.fastq.gz', 'ERR260247_1.fastq.gz', 'ERR260247_2.fastq.gz', 'ERR260248_1.fastq.gz', 'ERR260248_2.fastq.gz', 'ERR260249_1.fastq.gz', 'ERR260249_2.fastq.gz', 'ERR260250_1.fastq.gz', 'ERR260250_2.fastq.gz', 'ERR260251_1.fastq.gz', 'ERR260251_2.fastq.gz', 'ERR260252_1.fastq.gz', 'ERR260252_2.fastq.gz', 'ERR260253_1.fastq.gz', 'ERR260253_2.fastq.gz', 'ERR260254_1.fastq.gz', 'ERR260254_2.fastq.gz', 'ERR260255_1.fastq.gz', 'ERR260255_2.fastq.gz', 'ERR260256_1.fastq.gz', 'ERR260256_2.fastq.gz', 'ERR260257_1.fastq.gz', 'ERR260257_2.fastq.gz', 'ERR260258_1.fastq.gz', 'ERR260258_2.fastq.gz', 'ERR260259_1.fastq.gz', 'ERR260259_2.fastq.gz', 'ERR260260_1.fastq.gz', 'ERR260260_2.fastq.gz', 'ERR260261_1.fastq.gz', 'ERR260261_2.fastq.gz', 'ERR260262_1.fastq.gz', 'ERR260262_2.fastq.gz', 'ERR260263_1.fastq.gz', 'ERR260263_2.fastq.gz', 'ERR260264_1.fastq.gz', 'ERR260264_2.fastq.gz', 'ERR260265_1.fastq.gz', 'ERR260265_2.fastq.gz', 'ERR260266_1.fastq.gz', 'ERR260266_2.fastq.gz', 'ERR260267_1.fastq.gz', 'ERR260267_2.fastq.gz', 'ERR260268_1.fastq.gz', 'ERR260268_2.fastq.gz', 'ERR260269_1.fastq.gz', 'ERR260269_2.fastq.gz', 'ERR260270_1.fastq.gz', 'ERR260270_2.fastq.gz', 'ERR260271_1.fastq.gz', 'ERR260271_2.fastq.gz', 'ERR260272_1.fastq.gz', 'ERR260272_2.fastq.gz', 'ERR260273_1.fastq.gz', 'ERR260273_2.fastq.gz', 'ERR260274_1.fastq.gz', 'ERR260274_2.fastq.gz', 'ERR260275_1.fastq.gz', 'ERR260275_2.fastq.gz', 'ERR260276_1.fastq.gz', 'ERR260276_2.fastq.gz', 'ERR275251_1.fastq.gz', 'ERR275251_2.fastq.gz', 'ERR275252_1.fastq.gz', 'ERR275252_2.fastq.gz'])



if __name__ == '__main__':
    unittest.main()