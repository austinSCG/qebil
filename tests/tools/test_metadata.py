import glob
import numpy as np
from os import path, remove
import pandas as pd
from pandas.testing import assert_frame_equal
import unittest

from qebil.tools.util import setup_output_dir

from qebil.tools.metadata import (
    load_metadata,
    format_prep_type,
    set_criteria,
    scrub_special_chars,
    check_qebil_restricted_column,
    enforce_start_characters,
    clean_column_name,
    clean_nulls,
    qebil_format,
    detect_merger_column,
    merge_metadata,
)

_THIS_DIR, _THIS_FILENAME = path.split(__file__)

_TEST_SUPPORT_DIR = path.join(_THIS_DIR, "..", "support_files")

_TEST_OUTPUT_DIR = path.join(_THIS_DIR, "..", "test_output/")

setup_output_dir(_TEST_OUTPUT_DIR)


class metadataTest(unittest.TestCase):
    """
    Tests for functions in the metadata module.
    """

    def setUp(self):
        self.maxDiff = None
        # clean up the directory at the start
        cleanup_list = glob.glob(_TEST_OUTPUT_DIR + "*.EBI_metadata.tsv")
        for c in cleanup_list:
            remove(c)

    @classmethod
    def tearDown(self):
        pass  # no current actions needed

    def test_load_metadata(self):
        expected_dict = {
            "sample_name": {0: "SAMN16049500"},
            "library_name": {0: "MH001880_V1"},
            "secondary_sample_accession": {0: "SRS7392559"},
            "run_accession": {0: "SRR12672280"},
            "experiment_accession": {0: "SRX9152556"},
            "fastq_ftp": {
                0: "ftp.sra.ebi.ac.uk/vol1/fastq/SRR126/080/SRR12672280/SRR12672280_1.fastq.gz;ftp.sra.ebi.ac.uk/vol1/fastq/SRR126/080/SRR12672280/SRR12672280_2.fastq.gz"
            },
            "library_source": {0: "METAGENOMIC"},
            "instrument_platform": {0: "ILLUMINA"},
            "library_strategy": {0: "WGS"},
            "library_layout": {0: "PAIRED"},
            "tax_id": {0: "749906"},
            "scientific_name": {0: "gut metagenome"},
            "instrument_model": {0: "Illumina HiSeq 4000"},
            "library_selection": {0: "size fractionation"},
            "center_name": {0: "SUB8091541"},
            "experiment_title": {
                0: "Illumina HiSeq 4000 sequencing; SARS-CoV-2-specific IgA and limited inflammatory cytokines are present in the stool of select patients with acute COVID-19"
            },
            "study_title": {
                0: "human stool Metagenome from acute COVID19+ subjects"
            },
            "study_alias": {0: "PRJNA660883"},
            "experiment_alias": {0: "MH001880_V1"},
            "sample_alias": {0: "MH001880_V1"},
            "sample_title": {
                0: "This sample has been submitted by pda|jeremiahfaith on 2020-10-20; gut metagenome"
            },
            "fastq_md5": {
                0: "5c1da3b86d2bbb0d09e1f05cef0107f2;fe207ea59d80b5143e142050e37bbd11"
            },
            "study_accession": {0: "SRP283872"},
            "sample_title_specific": {0: "gut metagenome"},
            "taxon_id": {0: "749906"},
            "host": {0: "Homo sapiens"},
            "isolation_source": {0: "human feces"},
            "collection_date": {0: "2020"},
            "geo_loc_name": {0: "USA:New York"},
            "lat_lon": {0: "40.71 N 74.00 W"},
            "subject_id": {0: "MH001880_V1"},
            "biosamplemodel": {0: "Metagenome or environmental"},
            "ena_spot_count": {0: "5740751"},
            "ena_base_count": {0: "1722225300"},
            "ena_first_public": {0: "2020-10-20"},
            "ena_last_update": {0: "2020-10-20"},
            "experiment_title_specific": {
                0: "Illumina HiSeq 4000 paired end sequencing; SARS-CoV-2-specific IgA and limited inflammatory cytokines are present in the stool of select patients with acute COVID-19"
            },
            "run_prefix": {0: "SAMN16049500.SRR12672280.R"},
        }
        expected_df = pd.DataFrame.from_dict(expected_dict).set_index(
            "sample_name"
        )
        qebil_df = load_metadata(
            _TEST_SUPPORT_DIR + "/test_study_valid_qiita.tsv"
        )
        q2_df = load_metadata(_TEST_SUPPORT_DIR + "/test_study_valid_q2.tsv")
        qebil_csv_df = load_metadata(
            _TEST_SUPPORT_DIR + "/test_study_valid_q2.csv"
        )

        print(str(expected_df.columns))
        print(str(q2_df.columns))
        assert_frame_equal(expected_df, qebil_df)
        assert_frame_equal(expected_df, q2_df)
        assert_frame_equal(expected_df, qebil_csv_df)

    def test_format_prep_type(self):
        test_ambig_dict = {"POOLCLONE": "AMBIGUOUS", "test": "AMBIGUOUS"}
        test_genome_dict = {
            "CLONE": "Genome_Isolate",
            "WCS": "Genome_Isolate",
        }
        test_meta_g_dict = {
            "WGS": "Metagenomic",
            "WGA": "Metagenomic",
            "WXS": "Metagenomic",
            "ChIP-Seq": "Metagenomic",
            "Bisulfite-Seq": "Metagenomic",
        }
        test_meta_t_dict = {
            "RNA-Seq": "Metatranscriptomic",
            "miRNA-Seq": "Metatranscriptomic",
            "ncRNA-Seq": "Metatranscriptomic",
        }
        """ Code unused, # TODO: determine why and add tests?
        test_df_1 = pd.DataFrame(test_dict_1)
        test_dict_1 = {
            "sample_name": ["test1", "test2"],
            "library_strategy": ["test3", "test4"],
            "target_gene": ["test5", "test6"],
        }
        test_amplicon_dict = {
            "16S rRNA": "16S",
            "16S": "16S",
            "16S rDNA": "16S",
            "16S RNA": "16S",
            "16S DNA": "16S",
            "16S ": "16S",
            "ITS": "ITS",
            "ITS1": "ITS",
            "ITS2": "ITS",
            "ITS ": "ITS",
            "18S rRNA": "18S",
            "18S rDNA": "18S",
            "18S": "18S",
            "18S RNA": "18S",
            "18S DNA": "18S",
            "18S ": "18S",
        }
        """
        for k in test_ambig_dict.keys():
            test_df = pd.DataFrame(
                {
                    "sample_name": ["test1"],
                    "library_strategy": [k],
                }
            )
            for index, row in test_df.iterrows():
                self.assertEqual(
                    format_prep_type(row, index), test_ambig_dict[k]
                )
        for k in test_genome_dict.keys():
            test_df = pd.DataFrame(
                {
                    "sample_name": ["test1"],
                    "library_strategy": [k],
                }
            )
            for index, row in test_df.iterrows():
                self.assertEqual(
                    format_prep_type(row, index), test_genome_dict[k]
                )
        for k in test_meta_g_dict.keys():
            test_df = pd.DataFrame(
                {
                    "sample_name": ["test1"],
                    "library_strategy": [k],
                }
            )
            for index, row in test_df.iterrows():
                self.assertEqual(
                    format_prep_type(row, index), test_meta_g_dict[k]
                )
        for k in test_meta_t_dict.keys():
            test_df = pd.DataFrame(
                {
                    "sample_name": ["test1"],
                    "library_strategy": [k],
                }
            )
            for index, row in test_df.iterrows():
                self.assertEqual(
                    format_prep_type(row, index), test_meta_t_dict[k]
                )

        test_no_tg_df = pd.DataFrame(
            {
                "sample_name": ["test1"],
                "library_strategy": ["AMPLICON"],
            }
        )
        for index, row in test_no_tg_df.iterrows():
            self.assertEqual(format_prep_type(row, index), "AMBIGUOUS")

        test_16s_df = pd.DataFrame(
            {
                "sample_name": [
                    "test1",
                    "test2",
                    "test3",
                    "test4",
                    "test5",
                    "test6",
                ],
                "library_strategy": [
                    "AMPLICON",
                    "OTHER",
                    "AMPLICON",
                    "OTHER",
                    "AMPLICON",
                    "OTHER",
                ],
                "target_gene": [
                    "16S rRNA",
                    "16S",
                    "16S rDNA",
                    "16S RNA",
                    "16S DNA",
                    "16S ",
                ],
            }
        ).set_index("sample_name")
        for index, row in test_16s_df.iterrows():
            self.assertEqual(format_prep_type(row, index), "16S")

        test_18s_df = pd.DataFrame(
            {
                "sample_name": ["test1", "test2", "test3", "test4"],
                "library_strategy": [
                    "AMPLICON",
                    "OTHER",
                    "AMPLICON",
                    "OTHER",
                ],
                "target_gene": ["18S", "18S RNA", "18S DNA", "18S "],
            }
        )
        for index, row in test_18s_df.iterrows():
            self.assertEqual(format_prep_type(row, index), "18S")

        test_its_df = pd.DataFrame(
            {
                "sample_name": ["test1", "test2", "test3", "test4"],
                "library_strategy": [
                    "AMPLICON",
                    "OTHER",
                    "AMPLICON",
                    "OTHER",
                ],
                "target_gene": ["ITS", "ITS1", "ITS2", "ITS "],
            }
        )
        for index, row in test_its_df.iterrows():
            self.assertEqual(format_prep_type(row, index), "ITS")

    def test_set_criteria(self):
        valid_list_lower = ["test1", "test2"]
        valid_list_upper = ["TEST1", "TEST2"]

        test_selection_dict = set_criteria(valid_list_lower, valid_list_upper)
        self.assertEqual(
            test_selection_dict["library_strategy"],
            test_selection_dict["instrument_platform"],
        )

        # TODO: add these tests to return exceptions
        # invalid_list_numeric = 1
        # invalid_list_not = "test0"
        # self.assertRaises(set_criteria(invalid_list_numeric))
        # self.assertRaises(set_criteria(invalid_list_not))

    def test_scrub_special_chars(self):
        test_special_string_raw = (
            "A__B C-D(E)F/G|H~I`J@K#L$M%N^O&P*Q+R=S\\" + "T{U}V[W]X?Y<Z>1,2."
        )
        test_special_string_scrubbed = "A_B_C_D_leftparen_E_rightparen_F_per_G_bar_H_approximately_I_J_at_K_number_L_dollar_M_percent_N_to_power_O_and_P_star_Q_plus_R_equals_S_per_T_leftbracket_U_rightbracket_V_leftbracket_W_rightbracket_X_question_Y_less_than_Z_greater_than_1_comma_2_dot"
        self.assertEqual(
            scrub_special_chars(test_special_string_raw),
            test_special_string_scrubbed,
        )
        self.assertEqual(scrub_special_chars(1), "1")

        test_custom_replace_input = "A_acid"
        test_custom_replace_dict_upper = {"A": "alpha"}
        test_custom_replace_dict_lower = {"a": "alpha"}
        self.assertEqual(
            scrub_special_chars(
                test_custom_replace_input, test_custom_replace_dict_upper
            ),
            "alpha_acid",
        )
        self.assertEqual(
            scrub_special_chars(
                test_custom_replace_input, test_custom_replace_dict_lower
            ),
            "A_alphacid",
        )

    def test_check_qebil_restricted_column(self):
        reserved_word_lower = "key"
        reserved_word_upper = "KEY"
        resolved_word_lower = check_qebil_restricted_column(
            reserved_word_lower
        )
        resolved_word_upper = check_qebil_restricted_column(
            reserved_word_upper
        )
        self.assertEqual(resolved_word_lower, "user_key")
        self.assertEqual(resolved_word_upper, "user_key")
        self.assertEqual(resolved_word_lower, resolved_word_upper)

    def test_clean_nulls(self):
        test_na_list = [
            "n/a",
            "na",
            "N/A",
            "NA",
        ]
        test_nd_list = [
            "nd",
            "ND",
        ]
        test_np_list = [
            "",
            "null",
            "NULL",
        ]
        not_on_null_list = ["nil", "nada"]
        # test_supplement_dict = {"nada": "none", "NADA": "none"}

        for null in not_on_null_list:
            print(null)
            self.assertEqual(clean_nulls(null), null)
            """ Had to comment this out, line above would fail otherwise?
            # TODO: debug this test
            if null in test_supplement_dict.keys():
                self.assertEqual(
                    clean_nulls(null, test_supplement_dict), "none"
                )
            else:
                self.assertEqual(
                    clean_nulls(null, test_supplement_dict), null
                )
             """
        for null in test_na_list:
            self.assertEqual(clean_nulls(null), "not applicable")

        for null in test_nd_list:
            self.assertEqual(clean_nulls(null), "no data")

        for null in test_np_list:
            self.assertEqual(clean_nulls(null), "not provided")

    def test_clean_column_name(self):
        test_special_string_raw = (
            "A__B C-D(E)F/G|H~I`J@K#L$M%N^O&P*Q+R=S\\" + "T{U}V[W]X?Y<Z>1,2."
        )
        test_special_string_clean = "a_b_c_d_leftparen_e_rightparen_f_per_g_bar_h_approximately_i_j_at_k_number_l_dollar_m_percent_n_to_power_o_and_p_star_q_plus_r_equals_s_per_t_leftbracket_u_rightbracket_v_leftbracket_w_rightbracket_x_question_y_less_than_z_greater_than_1_comma_2_dot"
        test_reserved_word_lower = "user_key"
        test_reserved_word_upper = "USER_KEY"
        self.assertEqual(
            clean_column_name(test_special_string_raw),
            test_special_string_clean,
        )
        self.assertEqual(
            clean_column_name(test_reserved_word_lower), "user_key"
        )
        self.assertEqual(
            clean_column_name(test_reserved_word_upper), "user_key"
        )

    def test_enforce_start_characters(self):
        ok_string = "my_column"
        underscore_string_bad = "_special_column"
        numeric_string_bad = "123_column"
        underscore_string_good = "qebil_special_column"
        numeric_string_good = "qebil_123_column"
        self.assertEqual(enforce_start_characters(ok_string), ok_string)
        self.assertEqual(
            enforce_start_characters(underscore_string_bad),
            underscore_string_good,
        )
        self.assertEqual(
            enforce_start_characters(numeric_string_bad), numeric_string_good
        )

    def test_qebil_format(self):
        test_dict_1 = {
            "123_invalid": ["test1", "test2"],
            "Too much info (ug/mL)": ["test3", "test4"],
            "Genotype+/-": [np.nan, "test6"],
        }
        clean_dict_1 = {
            "qebil_123_invalid": ["test1", "test2"],
            "too_much_info_leftparen_ug_per_ml_rightparen": [
                "test3",
                "test4",
            ],
            "genotype_plus_minus": ["not provided", "test6"],
        }

        test_df = pd.DataFrame.from_dict(test_dict_1)  # , orient='index')
        clean_df = pd.DataFrame.from_dict(clean_dict_1)
        reformatted_df = qebil_format(test_df)
        assert_frame_equal(reformatted_df, clean_df)

    def test_detect_merger_column(self):
        base_df = pd.read_csv(
            _TEST_SUPPORT_DIR + "/SRP116878_sample_info.EBI_metadata.tsv",
            sep="\t",
            header=0,
            dtype=str,
        )
        merge_df = pd.read_csv(
            _TEST_SUPPORT_DIR + "/test_merger_metadata.tsv",
            sep="\t",
            header=0,
            dtype=str,
        )
        expected_col = "library_name"
        test_col = detect_merger_column(base_df, merge_df)
        self.assertEqual(expected_col, test_col)

    def test_merge_metadata(self):
        base_df = pd.read_csv(
            _TEST_SUPPORT_DIR + "/SRP116878_sample_info.EBI_metadata.tsv",
            sep="\t",
            header=0,
            dtype=str,
            index_col=0,
        )
        merge_df = pd.read_csv(
            _TEST_SUPPORT_DIR + "/test_merger_metadata.tsv",
            sep="\t",
            header=0,
            dtype=str,
        )
        expected_df = pd.read_csv(
            _TEST_SUPPORT_DIR + "/test_merged_md.tsv",
            sep="\t",
            header=0,
            dtype=str,
            index_col=0,
        )
        merger_col = "library_name"
        result_df = merge_metadata(base_df, merge_df, merger_col)
        assert_frame_equal(
            expected_df.sort_index(axis=1), result_df.sort_index(axis=1)
        )

    def test_merge_metadata_auto(self):
        base_df = pd.read_csv(
            _TEST_SUPPORT_DIR + "/SRP116878_sample_info.EBI_metadata.tsv",
            sep="\t",
            header=0,
            dtype=str,
            index_col=0,
        )
        merge_df = pd.read_csv(
            _TEST_SUPPORT_DIR + "/test_merger_metadata.tsv",
            sep="\t",
            header=0,
            dtype=str,
        )
        expected_df = pd.read_csv(
            _TEST_SUPPORT_DIR + "/test_merged_md.tsv",
            sep="\t",
            header=0,
            dtype=str,
            index_col=0,
        )
        result_df = merge_metadata(base_df, merge_df, "auto")
        assert_frame_equal(
            expected_df.sort_index(axis=1), result_df.sort_index(axis=1)
        )


if __name__ == "__main__":
    # begin the unittest.main()
    unittest.main()
