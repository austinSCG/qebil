from qebil.process import deplete_on_the_fly
from qebil.tools.fastq import get_read_count
from qebil.core import Study
import unittest
from os import path

from qebil.tools.util import setup_output_dir

_THIS_DIR, _THIS_FILENAME = path.split(__file__)

_TEST_SUPPORT_DIR = path.join(_THIS_DIR, "support_files")

_TEST_OUTPUT_DIR = path.join(_THIS_DIR, "test_output/")

setup_output_dir(_TEST_OUTPUT_DIR)
test_study_id = "SRP283872"


test_local_fastq_path = _TEST_SUPPORT_DIR + "/SRR13874871.fastq.gz"


class ProcessTest(unittest.TestCase):

    # https://click.palletsprojects.com/en/7.x/testing/#basic-testing
    def test_deplete_on_the_fly(self):
        """Tests on-the-fly host depletion and recovery"""
        test_study_id = "PRJNA629885"
        test_study = Study.from_remote(
            test_study_id, full_details=True, max_samples=2
        )
        dep_study = deplete_on_the_fly(
            test_study, output_dir=_TEST_OUTPUT_DIR, keep_files=True
        )
        md = dep_study.metadata
        for index in md.index:
            run_prefix = _TEST_OUTPUT_DIR + md.at[index, "run_prefix"]
            expected_raw_reads = str(md.at[index, "qebil_raw_reads"])
            expected_filtered_reads = str(
                md.at[index, "qebil_quality_filtered_reads"]
            )
            expected_mb_reads = str(md.at[index, "qebil_non_host_reads"])
            expected_filter_ratio = str(
                md.at[index, "qebil_frac_reads_passing_filter"]
            )
            expected_mb_ratio = str(md.at[index, "qebil_frac_non_host_reads"])
            actual_raw_reads = str(
                get_read_count(run_prefix + "1.ebi.fastq.gz")
            )
            actual_filtered_reads = str(
                get_read_count(run_prefix + "1.fastp.fastq.gz")
            )
            actual_mb_reads = str(
                get_read_count(run_prefix + "1.filtered.fastq.gz")
            )
            self.assertEqual(expected_raw_reads, actual_raw_reads)
            self.assertEqual(expected_filtered_reads, actual_filtered_reads)
            self.assertEqual(expected_mb_reads, actual_mb_reads)
            self.assertEqual(
                expected_filter_ratio,
                str(int(actual_filtered_reads) / int(actual_raw_reads)),
            )
            self.assertEqual(
                expected_mb_ratio,
                str(int(actual_mb_reads) / int(actual_filtered_reads)),
            )

        # now remove the host filtered reads, and run again, make sure the fastp files go away

        # need to create a study and supply the test_output dir, with and without keeping files
        # and make sure that the returned result dictionary has been populated with the expected
        # read counts and fractions of quality and microbial reads
        # for now, will read in expected final dataframe to compare, and run from remote with
        # small study. Need to also load from local partial files with only the first file finished
        # and second file still to do as well as with only a partially complete first file
        # need to also make sure that this behaves for mixed data types so that only metaG an metaT
        # files are processed, not amplicon. To do this, we'll load from a local file that is a mix
        # of the short metaG study and the short amplicon study used elsewhere in the tests
        # raise NotImplementedError

    def test_run_fastp(self):
        """Tests the running of fastp

        # need to supply a run_prefix, a list of raw_reads, instrument model and an output location
        # test keeping the raw files and removing them
        #raise NotImplementedError
        copy('test file','temp test file')
        paired_reads_prefix = ''
        read_1_only_prefix = ''
        too_many_reads_prefix = ''
        paired_read_count = 1
        read_1_only_count = 1
        expected_filtered_paired_read_count = 1
        expected_filtered_read_1_only_count = 1
        corrupt_file_prefix = ''

        self.assertEqual(None, run_fastp(too_many_reads_prefix,0,keep_files=True))
        self.assertEqual(None, run_fastp(paired_reads_prefix,0,keep_files=True))
        self.assertEqual('not determined',run_fastp(corrupt_file_prefix,0,keep_files=True))

        self.assertEqual(expected_filtered_paired_read_count,run_fastp(paired_reads_prefix,paired_read_count,keep_raw=True))
        self.assertEqual(expected_filtered_read_1_only_count,run_fastp(read_1_only_prefix,read_1_only_count,keep_raw=True))
        self.assertTrue(path.isdir(_TEST_OUTPUT_DIR + "fastp_reports/"))
        self.assertTrue(path.isfile('temp test file'))
        self.assertTrue(path.isfile('temp test file'))

        self.assertEqual(expected_filtered_paired_read_count,run_fastp(paired_reads_prefix,paired_read_count))
        self.assertEqual(expected_filtered_read_1_only_count,run_fastp(read_1_only_prefix,read_1_only_count))
        self.assertFalse(path.isfile('temp test file'))
        self.assertFalse(path.isfile('temp test file'))
        """

    def test_run_host_depletion(self):
        """Tests host depletion of fastq files

        #raise NotImplementedError
        copy('test file','temp test file')
        paired_reads_prefix = ''
        read_1_only_prefix = ''
        too_many_reads_prefix = ''
        paired_read_count = 1
        read_1_only_count = 1
        expected_filtered_paired_read_count = 1
        expected_filtered_read_1_only_count = 1
        corrupt_file_prefix = ''

        self.assertEqual(None, run_host_depletion(too_many_reads_prefix,0))
        self.assertEqual(None, run_host_depletion(paired_reads_prefix,0))
        self.assertEqual('not determined',run_host_depletion(corrupt_file_prefix,0))
        self.assertTrue(path.isdir(_TEST_OUTPUT_DIR + "fastp_reports/"))

        self.assertEqual(expected_filtered_paired_read_count,run_host_depletion(paired_reads_prefix,paired_read_count, keep=True))
        self.assertEqual(expected_filtered_read_1_only_count,run_host_depletion(read_1_only_prefix,read_1_only_count, keep=True))
        self.assertTrue(path.isfile('temp test file'))
        self.assertTrue(path.isfile('temp test file'))
        self.assertFalse(path.isfile(paired_reads_prefix + "_minimap2.start"))

        self.assertEqual(expected_filtered_paired_read_count,run_fastp(paired_reads_prefix,paired_read_count))
        self.assertEqual(expected_filtered_read_1_only_count,run_fastp(read_1_only_prefix,read_1_only_count))
        self.assertFalse(path.isfile('temp test file'))
        self.assertFalse(path.isfile('temp test file'))
        """


if __name__ == "__main__":
    # begin the unittest.main()
    unittest.main()

    """
    def deplete_on_the_fly(study, cpus=4, output_dir="./", keep_files=True):
    Processes a study sample-by-sample to filter and deplete human reads

    For human studies, it is convenient to on-the-fly host deplete for
    heavily host-contaminated samples. This function allows for this while
    performing lots of checks to enable interrupted runs to be resumed. For
    other sample types, please use either Qiita's built in host depletion
    pipeline or the Knight Lab/CMI scripts for bulk qc and host depletion

    Parameters
    ----------
    study: string
        EBI/ENA study ID to retrieve
    cpus: int
        number of processors available for running fastp and minimap2
    output_dir: string
         path to write the files to, and where already processed files should be
    keep_files: bool
        whether to retain raw (.ebi.) and intermediate (.fastp.) fastq files

    Returns
    ----------
    md: pd.DataFrame
        updated metadata with human and microbial read information

    ---------

    
    # We'll need to make sure the data type makes sense to host deplete
    # and if we're doing this, assume we'll want the metadata in Qiita
    # format
    if "qebil_prep_file" not in study.metadata.columns:
        study.populate_preps()

    md = study.metadata

    for index in md.index:
        run_prefix = md.at[index, "run_prefix"]
        strategy = md.at[index, "library_strategy"]
        model = md.at[index, "instrument_model"]
        try:
            raw_reads = str(md.at[index, "qebil_raw_reads"])
        except:
            raw_reads = "not determined"
        try:
            filtered_reads = str(md.at[index, "qebil_quality_filtered_reads"])
        except:
            filtered_reads = "not determined"

        try:
            mb_reads = str(md.at[index, "qebil_non_host_reads"])
        except:
            mb_reads = "not determined"

        prep_file_type = md.at[index, "qebil_prep_file"].split("_")[0]
        ebi_dict = unpack_fastq_ftp(
            md.at[index, "fastq_ftp"], md.at[index, "fastq_md5"]
        )

        if not raw_reads.isnumeric():
            # this value should only be set if downloaded
            md.at[index, "qebil_raw_reads"] = fetch_fastq_files(
                run_prefix, ebi_dict, output_dir
            )
            raw_reads = md.at[index, "qebil_raw_reads"]
        if not filtered_reads.isnumeric():
            # this value should only be set if filtered
            md.at[index, "qebil_quality_filtered_reads"] = run_fastp(
                run_prefix, raw_reads, model, output_dir, cpus, keep_files
            )
            filtered_reads = str(md.at[index, "qebil_quality_filtered_reads"])
            md.at[index, "qebil_frac_reads_passing_filter"] = (
                filtered_reads / raw_reads
            )
        # check to make sure the data can be host_depleted
        if prep_file_type not in ["Metagenomic", "Metatranscriptomic"]:
            logger.warning(
                "Skipping host depletion for "
                + run_prefix
                + " library_strategy: "
                + prep_file_type
                + " invalid for host depletion."
            )
        else:
            # this value should only be set if host filtering finished
            if not mb_reads.isnumeric():
                md.at[index, "qebil_non_host_reads"] = run_host_depletion(
                    run_prefix,
                    filtered_reads,
                    output_dir,
                    cpus,
                    keep_files,
                )
                md_reads = md.at[index, "qebil_non_host_reads"]
                md.at[index, "qebil_frac_non_host_reads"] = (
                    mb_reads / filtered_reads
                )

    # update study metadata
    return md


def run_fastp(
    run_prefix,
    raw_reads,
    model="",
    output_dir="./",
    cpus=4,
    keep_raw=False,
    min_length=45,
    raw_suffix=".ebi",
):
    Runs quality filtering of raw fastq(.gz) files

    This method takes in a list of fastq(.gz) files to quality filter using
    fastp. The minimum length of 45 bp is set to allow variance for short
    read RNAseq data. The instrument model is used to determine polyG
    filtering.

    Parameters
    ----------
    run_prefix: str
        prefix of fastq file names
    raw_reads: int
        number of reads expected in the input files
    model: str
        model of the instrument used to generate the fastq(s)
    cpus : int
        number of threads to use for fastp and fastq
    output_dir : string
        the path to save the output files
    keep : boolean
        whether to keep the raw files after filtering
    min_length : int
        minimum sequencing length for quality filtering
        default set to 45 to permit RNA-Seq data

    Returns
    -------
    filtered_reads: int
        the number of reads in the filtered file(s)

    
    polyG_model_list = [
        "",  # when in doubt, better to drop polyG..
        "Illumina MiniSeq",
        "Illumina NovaSeq 6000",
        "NextSeq 500",
        "NextSeq 550",
    ]

    fastp_dir = output_dir + "fastp_reports/"

    if not path.exists(fastp_dir):
        makedirs(fastp_dir)

    forward_reads = glob.glob(
        output_dir + run_prefix + ".R1" + raw_suffix + ".fastq.gz"
    )
    reverse_reads = glob.glob(
        output_dir + run_prefix + ".R2" + raw_suffix + ".fastq.gz"
    )

    if len(forward_reads) == 0:
        logger.warning(
            "No file found matching "
            + output_dir
            + "/"
            + run_prefix
            + ".R1"
            + raw_suffix
            + ".fastq.gz"
        )
    elif len(forward_reads) >= 2:
        logger.warning(
            "Too many files matching"
            + output_dir
            + "/"
            + run_prefix
            + ".R1"
            + raw_suffix
            + ".fastq.gz"
        )
    else:
        paired = False
        raw_read_1 = forward_reads[0]
        filtered_read_1 = raw_read_1.replace(raw_suffix, ".fastp")
        raw_read_2 = ""  # dummy file to help with checks and cleanup

        qc_json = fastp_dir + filtered_read_1.split("/")[-1].replace(
            ".R1.fastp.fastq.gz", ".fastp.json"
        )
        qc_html = fastp_dir + filtered_read_1.split("/")[-1].replace(
            ".R1.fastp.fastq.gz", ".fastp.html"
        )

        fastp_args = [
            "fastp",
            "-l",
            str(min_length),
            "-w",
            str(cpus),
            "-i",
            raw_read_1,
            "-j",
            qc_json,
            "-h",
            qc_html,
        ]
        if len(reverse_reads) >= 2:
            logger.warning(
                "Too many files matching"
                + output_dir
                + "/"
                + run_prefix
                + ".R2"
                + raw_suffix
                + ".fastq.gz"
            )
        elif len(reverse_reads) == 1:
            raw_read_2 = reverse_reads[0]
            filtered_read_2 = raw_read_2.replace(raw_suffix, ".fastp")
            fastp_args += ["-I", raw_read_2]
            paired = True

        if model in polyG_model_list:
            fastp_args += [
                "-g",
                "--poly_g_min_len",
                "10",  # polyG filtering, 10 is default
            ]

        # stage args for with and without filtering
        fastp_filter_args = fastp_args + ["-o", filtered_read_1]

        if paired:
            fastp_filter_args += ["-O", filtered_read_2]

        # now that files are staged check that the input matches expectations
        input_read_count = get_read_count(raw_read_1, raw_read_2)

        if str(input_read_count) != raw_reads:
            logger.warning(
                "Number of raw reads in file(s) "
                + raw_read_1
                + " "
                + raw_read_2
                + " "
                + str(input_read_count)
                + "does not match expected: "
                + str(raw_reads)
                + "Skipping fastp."
            )
        else:
            logger.info(fastp_filter_args)
            fastp_filter_ps = Popen(fastp_filter_args)
            fastp_filter_ps.wait()
            if not path.isfile(qc_json):
                logger.warning(
                    "No quality filtered file found for: "
                    + run_prefix
                    + " Quality filtering evaluation."
                )
                return "not determined"
            else:
                logger.info("Fastp finished for " + run_prefix)
                with open(qc_json) as json_file:
                    fastp_results = json.load(json_file)
                summary = fastp_results["summary"]
                filtered_reads = summary["after_filtering"]["total_reads"]

                if not keep_raw:
                    for f in [raw_read_1, raw_read_2]:
                        if path.isfile(f):
                            remove(f)

                return filtered_reads


def run_host_depletion(
    run_prefix,
    filtered_reads,
    output_dir="./",
    cpus=4,
    keep=False,
    filter_db="/databases/minimap2/human-phix-db.mmi",
):
    Quality filter and host deplete list of files

    Takes fastp filtered files and runs minimap2 on them
    to remove human reads. Other databases can be supplied
    instead so long as there is a valid path to the .mmi file

    N.B. the default parameters will only work on barnacle and
    assume the location of the human+PhiX database as shown above.

    Parameters
    ----------
    run_prefix: str
         prefix of fastq file names
    filtered_reads: int
         number of reads expected in the input files
    output_dir : string
        the path to save the output files
    cpus : int
        number of threads to use for fastp and fastq
    keep : boolean
        whether to keep the raw files after quality filtering,
        and fastp files after host filtering
    filter_db : string
        path to the mmi database to be used for filtering

    Returns
    -------
    non_host_reads_dict : dict
        dict of of non human reads per sample
    

    # list to be used to confirm things worked as expected
    confirm_list = []

    forward_reads = glob.glob(
        output_dir + "/" + run_prefix + ".R1.fastp.fastq.gz"
    )
    reverse_reads = glob.glob(
        output_dir + "/" + run_prefix + ".R2.fastp.fastq.gz"
    )

    if len(forward_reads) == 0:
        logger.warning(
            "No file found matching "
            + output_dir
            + "/"
            + run_prefix
            + ".R1.fastp.fastq.gz"
        )
    elif len(forward_reads) >= 2:
        logger.warning(
            "Too many files matching"
            + output_dir
            + "/"
            + run_prefix
            + ".R1.fastp.fastq.gz"
        )
    else:
        paired = False
        fastp_read_1 = forward_reads[0]
        hd_read_1 = fastp_read_1.replace(".fastp", ".filtered")
        confirm_list.append(hd_read_1)
        fastp_read_2 = ""  # dummy file to help with checks and cleanup

        # set minimap2 arguments
        minimap2_args = [
            "minimap2",
            "-ax",
            "sr",
            "-t",
            str(cpus),
            filter_db,
            "-a",
            fastp_read_1,
        ]

        if len(reverse_reads) >= 2:
            logger.warning(
                "Too many files matching"
                + output_dir
                + "/"
                + run_prefix
                + ".R2.fastp.fastq.gz"
            )
        elif len(reverse_reads) == 1:
            fastp_read_2 = reverse_reads[0]
            hd_read_2 = fastp_read_2.replace(".fastp", ".filtered")
            confirm_list.append(hd_read_2)
            paired = True

        stf_args = ["samtools", "fastq", "-@", str(cpus), "-F", "256", "-"]

        if paired:
            minimap2_args += [fastp_read_2]
            stf_args = stf_args + [
                "-f",
                "12",
                "-1",
                hd_read_1,
                "-2",
                hd_read_2,
            ]
        else:  # different samtools parameter for unpaired
            stf_args = stf_args + ["-f", "4", "-0", hd_read_1]

        # now that files are staged check that the input matches expectations
        input_read_count = get_read_count(fastp_read_1, fastp_read_2)
        if str(input_read_count) != str(filtered_reads):
            logger.warning(
                "Number of reads in file(s) "
                + fastp_read_1
                + " "
                + fastp_read_2
                + " "
                + str(input_read_count)
                + "does not match expected: "
                + str(filtered_reads)
                + "Skipping minimap2."
            )
        else:
            logger.info("Starting minimap2 depletion.")
            f = open(run_prefix + "_minimap2.start", "a")
            f.close()

            minimap2_ps = Popen(minimap2_args, stdout=PIPE, stderr=PIPE)
            stf_ps = Popen(
                stf_args, stdin=minimap2_ps.stdout, stdout=PIPE, stderr=PIPE
            )
            stf_ps.wait()
            remove(run_prefix + "_minimap2.start")

        # now make sure the output files are valid with fastqs
        confirm = True
        for c in confirm_list:
            if confirm:
                confirm = check_valid_fastq(c, cpus, output_dir)

        if confirm:
            if not keep:
                for r in [fastp_read_1, fastp_read_2]:
                    if path.isfile(r):
                        remove(r)
            logger.info("minimap2 success.")
            if paired:
                mb_reads = get_read_count(hd_read_1, hd_read_2)
            else:
                mb_reads = get_read_count(hd_read_1)
        else:
            logger.warning(
                "minimap2 failed."
                + " Corrupt or invalid fastq"
                + " generated. Removing files"
            )
            for c in confirm_list:
                if path.isfile(c):
                    remove(c)
            mb_reads = "minimap2 error"

    return mb_reads
    """
