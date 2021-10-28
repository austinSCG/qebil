import unittest
import pandas as pd
from os import path
from collections import OrderedDict
from pandas.testing import assert_frame_equal

# This is the class we want to test. So, we need to import it
from qebil.core import Study
from qebil.tools.metadata import load_metadata

_THIS_DIR, _THIS_FILENAME = path.split(__file__)

_TEST_SUPPORT_DIR = path.join(_THIS_DIR, "support_files")


def assert_frame_not_equal(*args, **kwargs):
    try:
        assert_frame_equal(*args, **kwargs)
    except AssertionError:
        # frames are not equal
        pass
    else:
        # frames are equal
        raise AssertionError


class StudyTest(unittest.TestCase):
    """
    The basic class that inherits unittest.TestCase
    """

    study = Study()  # instantiate the Study Class
    study_id_list = []  # variable that stores obtained study_id
    md_dict = {}  # variable that stores study metadata as dict

    md_dict_empty = {}
    md_dict_basic = {
        "study_id": ["TESTING123"],
        "study_accession": ["PRJXXXX"],
        "sample_name": ["sample1"],
    }
    md_dict_invalid_sample = {
        "study_id": ["TESTING123"],
        "study_accession": ["PRJXXXX"],
        "sample_name": ["sample-1"],
    }
    md_dict_list = [md_dict_basic, md_dict_invalid_sample]

    test_study_valid_id = "SRP283872"
    test_study_valid_tsv = _TEST_SUPPORT_DIR + "/test_study_valid.tsv"
    test_prep_valid_tsv = _TEST_SUPPORT_DIR + "/test_prep_valid.tsv"
    test_study_invalid_id = "NOT A STUDY"
    test_study_numeric = 123
    test_study_mixed = "SRP276821"

    test_modes = ["ebi"]  # 'ncbi', 'invalid']

    test_fields_invalid = ["Larry", "Curly", "Moe"]
    test_fields_string = "test example"
    test_fields_defaults = [
        "library_name",
        "secondary_sample_accession",
        "run_accession",
        "experiment_accession",
        "fastq_bytes" "fastq_ftp",
        "library_source",
        "instrument_platform",
        "submitted_format",
        "library_strategy",
        "library_layout",
        "tax_id",
        "scientific_name",
        "instrument_model",
        "library_selection",
        "center_name",
        "experiment_title",
        "study_title",
        "study_alias",
        "experiment_alias",
        "sample_alias",
        "sample_title",
        "fastq_md5",
    ]

    # test case function to check the Study.set_name function
    def test_study_empty(self):
        study_empty = Study()
        assert_frame_equal(study_empty.metadata, pd.DataFrame())

    """ Not used so consider removing
    def test_study_creation_from_dict(self):
        for m in self.md_dict_list:
            study_from_dict = Study.from_dict(m)
            assert_frame_equal(study_from_dict.metadata,
                                pd.DataFrame.from_dict(m))
    """

    def test_study_creation_from_dataframe(self):
        for m in self.md_dict_list:
            from_dict_df = pd.DataFrame.from_dict(m)
            study_from_dataframe = Study(from_dict_df)
            assert_frame_equal(
                study_from_dataframe.metadata.sort_index(axis=1),
                from_dict_df.sort_index(axis=1),
            )

    def test_study_creation_from_remote_valid(self):
        valid_study = Study.from_remote(
            self.test_study_valid_id, full_details=True, max_samples=2
        )
        # valid_study.metadata.to_csv('./new_test.tsv', sep='\t', index_label='sample_name')
        test_df = pd.read_csv(
            self.test_study_valid_tsv, header=0, sep="\t", dtype=str
        )
        test_df["sample_name"] = test_df["sample_accession"]
        test_df = test_df.set_index("sample_name")
        # print(valid_study.metadata.columns)
        # print(test_df.columns)

        assert_frame_equal(
            valid_study.metadata.sort_index(axis=1),
            test_df.sort_index(axis=1),
        )
        self.assertEqual(valid_study.study_id, self.test_study_valid_id)

    def test_study_creation_from_remote_invalid_fields(self):
        test_df = pd.read_csv(
            self.test_study_valid_tsv, header=0, sep="\t", dtype=str
        )
        valid_study = Study.from_remote(
            self.test_study_valid_id, fields=self.test_fields_invalid
        )
        assert_frame_not_equal(
            valid_study.metadata.sort_index(axis=1),
            test_df.sort_index(axis=1),
        )
        self.assertEqual(valid_study.ebi_id, self.test_study_valid_id)
        for f in self.test_fields_invalid:
            self.assertFalse(f in valid_study.metadata.columns)

    def test_populate_sample_names(self):
        test_study = Study.from_remote(self.test_study_valid_id)
        test_study.populate_sample_names()
        self.assertTrue("run_prefix" in test_study.metadata.columns)

        # TODO: add test to ensure name resolution when sample_name is same
        # for all samples need to find study ID with this issue to use

    def test_populate_sample_info(self):
        test_study = Study.from_remote(
            self.test_study_valid_id, full_details=True, max_samples=2
        )
        test_df = load_metadata(self.test_study_valid_tsv)
        # print(test_study.metadata.columns)
        # print(test_df.columns)
        self.assertEqual(test_study.metadata.shape, test_df.shape)

    def test_populate_expt_info_and_preps(self):
        test_study = Study.from_remote(
            self.test_study_valid_id, full_details=True, max_samples=2
        )
        test_study.populate_preps()
        test_df = pd.read_csv(
            self.test_prep_valid_tsv,
            header=0,
            sep="\t",
            dtype=str,
            index_col=0,
        )
        prep_cols = test_df.columns
        # print("prep_cols:" + str(prep_cols))
        prep_subset = test_study.metadata[prep_cols]
        prep_subset.index = prep_subset.index.rename("sample_name")
        print("prep_subset:" + prep_subset.head(1))
        print("prep_subset:" + test_df.head(1))
        assert_frame_equal(
            prep_subset.sort_index(axis=1), test_df.sort_index(axis=1)
        )

    def test_filter_samples(self):
        test_study = Study.from_remote(
            self.test_study_mixed, full_details=True
        )
        filter_dict = {"library_strategy": ["wgs"]}
        self.assertEqual(len(test_study.metadata), 8)
        test_study.filter_samples(filter_dict)
        self.assertEqual(len(test_study.metadata), 2)

    def test_populate_details(self):
        # Note, this is almost identical to test_fetch_ebi_info_study,
        # so maybe refactor to check for the log messages?
        accession = "SRP283872"
        test_study_dict = OrderedDict(
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
        test_study = Study()
        test_study.ebi_id = accession
        test_study.populate_details()
        self.assertEqual(test_study.details, test_study_dict)


if __name__ == "__main__":
    # begin the unittest.main()
    unittest.main()
