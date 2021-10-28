from collections import OrderedDict
from os import path, remove
from shutil import copy
import pandas as pd
import unittest

from qebil.tools.fastq import get_fastq_head

from qebil.tools.util import (
    get_checksum,
    detect_qiita_study,
    get_ebi_ids,
    load_project_file,
    parse_document,
    parse_details,
    remove_index_read_file,
    scrape_ebi_ids,
    scrape_seq_method,
    setup_output_dir,
    unpack_fastq_ftp,
)


_THIS_DIR, _THIS_FILENAME = path.split(__file__)

_TEST_SUPPORT_DIR = path.join(_THIS_DIR, "..", "support_files")

_TEST_OUTPUT_DIR = path.join(_THIS_DIR, "..", "test_output/")

_CLEANUP_LIST = []

_TEST_STUDY_DICT = OrderedDict(
    [
        (
            "STUDY_SET",
            OrderedDict(
                [
                    (
                        "STUDY",
                        OrderedDict(
                            [
                                ("@accession", "SRP283872"),
                                ("@alias", "PRJNA660883"),
                                ("@center_name", "BioProject"),
                                ("@broker_name", "NCBI"),
                                (
                                    "IDENTIFIERS",
                                    OrderedDict(
                                        [
                                            (
                                                "PRIMARY_ID",
                                                "SRP283872",
                                            ),
                                            (
                                                "SECONDARY_ID",
                                                "PRJNA660883",
                                            ),
                                            (
                                                "EXTERNAL_ID",
                                                OrderedDict(
                                                    [
                                                        (
                                                            "@label",
                                                            "primary",
                                                        ),
                                                        (
                                                            "@namespace",
                                                            "BioProject",
                                                        ),
                                                        (
                                                            "#text",
                                                            "PRJNA660883",
                                                        ),
                                                    ]
                                                ),
                                            ),
                                            (
                                                "SUBMITTER_ID",
                                                OrderedDict(
                                                    [
                                                        (
                                                            "@namespace",
                                                            "BioProject",
                                                        ),
                                                        (
                                                            "#text",
                                                            "PRJNA660883",
                                                        ),
                                                    ]
                                                ),
                                            ),
                                        ]
                                    ),
                                ),
                                (
                                    "DESCRIPTOR",
                                    OrderedDict(
                                        [
                                            (
                                                "STUDY_TITLE",
                                                "human stool Metagenome from acute COVID19+ subjects",
                                            ),
                                            (
                                                "STUDY_TYPE",
                                                OrderedDict(
                                                    [
                                                        (
                                                            "@existing_study_type",
                                                            "Metagenomics",
                                                        )
                                                    ]
                                                ),
                                            ),
                                            (
                                                "STUDY_ABSTRACT",
                                                "16S and metagenomic sequences from individuals with acute COVID19",
                                            ),
                                            (
                                                "CENTER_PROJECT_NAME",
                                                "human stool",
                                            ),
                                        ]
                                    ),
                                ),
                                (
                                    "STUDY_LINKS",
                                    OrderedDict(
                                        [
                                            (
                                                "STUDY_LINK",
                                                [
                                                    OrderedDict(
                                                        [
                                                            (
                                                                "XREF_LINK",
                                                                OrderedDict(
                                                                    [
                                                                        (
                                                                            "DB",
                                                                            "ENA-SAMPLE",
                                                                        ),
                                                                        (
                                                                            "ID",
                                                                            "SRS7392524-SRS7392559",
                                                                        ),
                                                                    ]
                                                                ),
                                                            )
                                                        ]
                                                    ),
                                                    OrderedDict(
                                                        [
                                                            (
                                                                "XREF_LINK",
                                                                OrderedDict(
                                                                    [
                                                                        (
                                                                            "DB",
                                                                            "ENA-EXPERIMENT",
                                                                        ),
                                                                        (
                                                                            "ID",
                                                                            "SRX9152521-SRX9152556",
                                                                        ),
                                                                    ]
                                                                ),
                                                            )
                                                        ]
                                                    ),
                                                    OrderedDict(
                                                        [
                                                            (
                                                                "XREF_LINK",
                                                                OrderedDict(
                                                                    [
                                                                        (
                                                                            "DB",
                                                                            "ENA-RUN",
                                                                        ),
                                                                        (
                                                                            "ID",
                                                                            "SRR12672280-SRR12672315",
                                                                        ),
                                                                    ]
                                                                ),
                                                            )
                                                        ]
                                                    ),
                                                    OrderedDict(
                                                        [
                                                            (
                                                                "XREF_LINK",
                                                                OrderedDict(
                                                                    [
                                                                        (
                                                                            "DB",
                                                                            "ENA-SUBMISSION",
                                                                        ),
                                                                        (
                                                                            "ID",
                                                                            "SRA1127885",
                                                                        ),
                                                                    ]
                                                                ),
                                                            )
                                                        ]
                                                    ),
                                                    OrderedDict(
                                                        [
                                                            (
                                                                "XREF_LINK",
                                                                OrderedDict(
                                                                    [
                                                                        (
                                                                            "DB",
                                                                            "ENA-FASTQ-FILES",
                                                                        ),
                                                                        (
                                                                            "ID",
                                                                            "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=SRP283872&result=read_run&fields=run_accession,fastq_ftp,fastq_md5,fastq_bytes",
                                                                        ),
                                                                    ]
                                                                ),
                                                            )
                                                        ]
                                                    ),
                                                    OrderedDict(
                                                        [
                                                            (
                                                                "XREF_LINK",
                                                                OrderedDict(
                                                                    [
                                                                        (
                                                                            "DB",
                                                                            "ENA-SUBMITTED-FILES",
                                                                        ),
                                                                        (
                                                                            "ID",
                                                                            "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=SRP283872&result=read_run&fields=run_accession,submitted_ftp,submitted_md5,submitted_bytes,submitted_format",
                                                                        ),
                                                                    ]
                                                                ),
                                                            )
                                                        ]
                                                    ),
                                                ],
                                            )
                                        ]
                                    ),
                                ),
                                (
                                    "STUDY_ATTRIBUTES",
                                    OrderedDict(
                                        [
                                            (
                                                "STUDY_ATTRIBUTE",
                                                [
                                                    OrderedDict(
                                                        [
                                                            (
                                                                "TAG",
                                                                "ENA-SPOT-COUNT",
                                                            ),
                                                            (
                                                                "VALUE",
                                                                "166909868",
                                                            ),
                                                        ]
                                                    ),
                                                    OrderedDict(
                                                        [
                                                            (
                                                                "TAG",
                                                                "ENA-BASE-COUNT",
                                                            ),
                                                            (
                                                                "VALUE",
                                                                "50072960400",
                                                            ),
                                                        ]
                                                    ),
                                                    OrderedDict(
                                                        [
                                                            (
                                                                "TAG",
                                                                "ENA-FIRST-PUBLIC",
                                                            ),
                                                            (
                                                                "VALUE",
                                                                "2020-10-20",
                                                            ),
                                                        ]
                                                    ),
                                                    OrderedDict(
                                                        [
                                                            (
                                                                "TAG",
                                                                "ENA-LAST-UPDATE",
                                                            ),
                                                            (
                                                                "VALUE",
                                                                "2020-10-20",
                                                            ),
                                                        ]
                                                    ),
                                                ],
                                            )
                                        ]
                                    ),
                                ),
                            ]
                        ),
                    )
                ]
            ),
        )
    ]
)


class utilTest(unittest.TestCase):
    """
    Tests for functions in the util module.
    """

    def setUp(self):
        self.maxDiff = None  # helps with troubleshooting
        setup_output_dir(_TEST_OUTPUT_DIR)

    def tearDown(self):
        # clean up the directory at the end
        for c in _CLEANUP_LIST:
            if path.isfile(c):
                remove(c)

    def test_load_project_file(self):
        test_project_file = _TEST_SUPPORT_DIR + "/test_project_file.tsv"
        test_project_file_no_header = (
            _TEST_SUPPORT_DIR + "/test_project_file_no_header.tsv"
        )
        test_load_project_list = load_project_file(test_project_file)
        test_load_project_list_no_header = load_project_file(
            test_project_file_no_header
        )

        self.assertEqual(len(test_load_project_list), 2)
        self.assertEqual(len(test_load_project_list_no_header), 2)
        self.assertEqual(
            test_load_project_list, test_load_project_list_no_header
        )

    def test_get_checksum(self):
        test_local_fastq_path = _TEST_SUPPORT_DIR + "/SRR13874871.fastq.gz"
        expected_md5 = "06445ed5341e3779ac1d5230c787c538"
        self.assertTrue(path.isfile(test_local_fastq_path))
        md5checksum_res = get_checksum(test_local_fastq_path, expected_md5)
        self.assertEqual(expected_md5, md5checksum_res)

    def test_parse_document(self):
        tokens_html = parse_document(
            "https://www.nature.com/articles/s41598-020-60564-8"
        )
        tokens_pdf = parse_document(
            "https://www.nature.com/articles/s41598-021-83922-6.pdf"
        )

        self.assertEqual(len(tokens_html), 85235)  # this changes a lot...
        self.assertEqual(len(tokens_pdf), 15771)

    def test_scrape_ebi_ids(self):
        test_paths = [
            "https://bmccancer.biomedcentral.com/track/pdf/10.1186/s12885-020-6654-5.pdf",
            "https://bmccancer.biomedcentral.com/articles/10.1186/s12885-020-6654-5",
            _TEST_SUPPORT_DIR
            + "/test.pdf",  # contains PRJNA, but now invalid ID. Update?
        ]
        test_ids = [["PRJNA533024"], ["PRJNA533024"], []]
        test_ids_retrieved = []

        for f in test_paths:
            parsed_doc = parse_document(f)
            for p in parsed_doc:
                if "PRJNA" in p:
                    print(p)
            test_ids_retrieved.append(scrape_ebi_ids(parsed_doc))

        self.assertEqual(test_ids, test_ids_retrieved)

    def test_get_ebi_ids(self):
        expected_tuple = "SRP283872", "PRJNA660883"
        test_id_tuple = get_ebi_ids(_TEST_STUDY_DICT)
        self.assertEqual(expected_tuple, test_id_tuple)

    def test_setup_output_dir(self):
        # called in setUp so just make sure it worked
        self.assertTrue(path.isdir(_TEST_OUTPUT_DIR))

    def test_detect_qiita_study(self):
        expected_studies = ["test3", "test4"]
        test_dict_1 = {
            "sample_name": ["test1", "test2"],
            "qiita_study_id": ["test3", "test4"],
        }

        test_df_1 = pd.DataFrame(test_dict_1)
        test_df_2 = pd.DataFrame()
        self.assertEqual(
            set(expected_studies), set(detect_qiita_study(test_df_1))
        )
        self.assertFalse(detect_qiita_study(test_df_2))

    def test_parse_details(self):
        expected_dict = {
            "abstract": "16S and metagenomic sequences from individuals with acute COVID19",
            "description": "PRJNA660883",
            "title": "human stool Metagenome from acute COVID19+ subjects",
            "seq_method": ["16s"],
        }
        parsed_dict = parse_details(_TEST_STUDY_DICT)
        self.assertEqual(expected_dict, parsed_dict)

    def test_scrape_seq_method(self):
        test_string_1 = "16S and metagenomic sequences from individuals with acute COVID19"
        test_string_2 = "ITS a beautiful mix of methods 16s"
        test_string_3 = "ITS1 18S"
        self.assertEqual(
            scrape_seq_method(test_string_1), scrape_seq_method(test_string_2)
        )
        self.assertNotEqual(
            scrape_seq_method(test_string_2), scrape_seq_method(test_string_3)
        )

    def test_unpack_fastq_ftp(self):
        test_fastq_ftp_string = (
            "ftp.sra.ebi.ac.uk/vol1/fastq/SRR126/080/"
            + "SRR12672280/SRR12672280_1.fastq.gz;"
            + "ftp.sra.ebi.ac.uk/vol1/fastq/SRR126/080/"
            + "SRR12672280/SRR12672280_2.fastq.gz"
        )
        test_fastq_md5_string = (
            "5c1da3b86d2bbb0d09e1f05cef0107f2;"
            + "fe207ea59d80b5143e142050e37bbd11"
        )

        test_fastq_bytes_string = "411122459;446055180"
        test_unpack_result = unpack_fastq_ftp(
            test_fastq_ftp_string,
            test_fastq_md5_string,
            test_fastq_bytes_string,
            2,
        )
        test_read_dict = {
            "read_1": {
                "ftp": test_fastq_ftp_string.split(";")[0],
                "md5": test_fastq_md5_string.split(";")[0],
                "bytes": 411122459,
            },
            "read_2": {
                "ftp": test_fastq_ftp_string.split(";")[1],
                "md5": test_fastq_md5_string.split(";")[1],
                "bytes": 446055180,
            },
        }
        self.assertEqual(test_unpack_result[0], test_read_dict)

    def test_remove_index_read_file(self):
        # set up files
        test_dict = {}
        test_fastq_path = _TEST_SUPPORT_DIR + "/SRR13874871.fastq.gz"

        for r in range(1, 4):
            fq_src = (
                _TEST_SUPPORT_DIR
                + "/SAMN07663020.SRR6050387.R"
                + str(r)
                + ".ebi.fastq.gz"
            )
            fq_test = (
                _TEST_OUTPUT_DIR
                + "/SAMN07663020.SRR6050387.R"
                + str(r)
                + ".ebi.fastq.gz"
            )
            copy(fq_src, fq_test)
            test_dict["read" + str(r)] = {}
            test_dict["read" + str(r)]["fp"] = fq_test
            test_dict["read" + str(r)]["md5"] = get_checksum(fq_test)

        # print(test_fastq_path)
        """
        for r in range(1, 4):
            r_clone = test_fastq_path.replace(
                ".fastq.gz", ".R" + str(r) + ".ebi.fastq.gz"
            )
            get_fastq_head(test_fastq_path, r_clone, r * 10)
            test_dict["read"+str(r)]= r_clone
        """
        # run test
        remove_index_read_file(test_dict, "PAIRED")

        # confirm behavior
        self.assertTrue(path.isfile(test_dict["read1"]["fp"]))
        self.assertTrue(path.isfile(test_dict["read2"]["fp"]))
        self.assertFalse("read3" in test_dict.keys())


if __name__ == "__main__":
    # begin the unittest.main()
    unittest.main()
