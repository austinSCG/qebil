import hashlib
from os import makedirs, path, remove
import pandas as pd
import PyPDF2
import re
import requests
from requests.exceptions import RequestException
import urllib
from urllib.request import urlretrieve

from qebil.log import logger

THIS_DIR, THIS_FILENAME = path.split(__file__)


def load_project_file(filepath):
    """Reads in a project file to populate a list of study ids

    This method looks for a 'study_id' column in a tsv-formatted
    text file, as generated by the search function of QEBIL, and
    checks the first column if not found. This latter catch allows
    the user to supply a list of project ids, one per line to fetch.

    Parameters
    ----------
    filepath: string:
        path to the project file with a 'study_id' column

    Returns
    ---------
    add_list: list
        list of EBI project IDs found in the file

    """
    proj_df = pd.read_csv(filepath, sep="\t", header=0)
    if "study_id" not in proj_df.columns:
        proj_df = pd.read_csv(filepath, sep="\t", header=None)
        proj_df.columns = ["study_id"]
    add_list = list(proj_df["study_id"].unique())
    add_list = [x.strip() for x in add_list]

    return add_list


def retrieve_ftp_file(ftp_path, filepath, remote_checksum, overwrite=False):
    """Method to retrieve an ftp file and check accuracy
    of the download. If not overwriting, check the local copy for validity
    before downloading.

    Parameters
    ----------
    ftp_path: string:
        the ftp url to download from
    filepath: string
        the local path to save the file to
    remote_checksum: string
        hexadecimal md5 checksum to validate the downlod
    overwrite : bool
        whether to overwrite the local copy of the file

    Returns
    ---------
    checksum: str or bool
        returns either the string for the checksum or False
        if there is an issue with the download or integrity

    """
    # see if the file exists and make sure it's valid if so
    local = False
    if not overwrite:
        if path.isfile(filepath):
            checksum = check_download_integrity(filepath, remote_checksum)
            if checksum:
                local = True
                logger.info(
                    "Valid file found. Skipping download of file: "
                    + "ftp://"
                    + ftp_path
                    + " to "
                    + filepath
                )
                return checksum
            else:
                logger.warning(
                    "Local file found but invalid checksum."
                    + " Downloading again: "
                    + "ftp://"
                    + ftp_path
                    + " to "
                    + filepath
                )

    if not local:
        # add catch in case there is an issue with the connection
        try:
            urlretrieve("ftp://" + ftp_path, filepath)
            checksum = check_download_integrity(filepath, remote_checksum)
            return checksum
        except urllib.error.URLError: #RequestException:
            logger.warning(
                "Issue with urlretrieve for "
                + "download of file:"
                + "ftp://"
                + ftp_path
                + " to "
                + filepath
            )
            # cleanup file if present
            if path.isfile(ftp_path):
                remove(ftp_path)
            return False


def check_download_integrity(filepath, md5_value):
    """Compares the checksum of local and remote files for validation

    This function compares the expected md5 check sum based on the
    remote file information with the md5 checksum from the specified
    locatl file. Returns False if they do not match or the checksum
    value if they do.

    Parameters
    ----------
    filepath: string
        the path to the file or URL to be scraped

    Returns
    ---------
    False OR local_checksum:string
        checksum of local file if valid, otherwise False
    """

    if path.isfile(filepath):
        fq = open(filepath, "rb")
        fq_contents = fq.read()
        fq.close()
        local_checksum = hashlib.md5(fq_contents).hexdigest()

        if md5_value != local_checksum:
            return False
        else:
            return local_checksum
    else:
        logger.warning(filepath + " not found.")
        return False


def parse_document(filepath):
    """This method allows local files or URLs to be parsed

    Uses basic NLP to parse a html page or file (PDF or text)
    for its contents

    Parameters
    ----------
    filepath: string
        the path to the file or URL to be scraped

    Returns:
    ----------
    tokens: list
        list of parsed tokens
    """
    full_text = ""
    tokens = []

    if filepath[0:3] == "10.1":
        # this is a doi, so just prepend and go with it
        filepath = "https://doi.org/" + filepath

    if filepath[0:3] == "htt":
        request = requests.get(filepath)
        if request.status_code == 200:
            if filepath[-3:] == "pdf":
                file_name, headers = urllib.request.urlretrieve(filepath)
                document = PyPDF2.PdfFileReader(open(file_name, "rb"))
                for i in range(document.numPages):
                    page_to_print = document.getPage(i)
                    full_text += page_to_print.extractText()
                tokens = [
                    t for t in re.split(r"\; |\, |\. | |\n|\t", full_text)
                ]
            else:
                full_text = request.text
                tokens = [
                    t
                    for t in re.split(
                        r"\;|\,|\.| |\n|\t|<|>|\/|\"|'", full_text
                    )
                ]
    else:
        if filepath[-3:] == "pdf":
            document = PyPDF2.PdfFileReader(open(filepath, "rb"))
            for i in range(document.numPages):
                page_to_print = document.getPage(i)
                full_text += page_to_print.extractTexit()
            tokens = [t for t in re.split(r"\; |\, |\. | |\n|\t", full_text)]
        else:
            text_file = open(filepath, "r")
            full_text = text_file.read()
            text_file.close()
            tokens = [
                t
                for t in re.split(r"\;|\,|\.| |\n|\t|<|>|\/|\"|'", full_text)
            ]

    return tokens


def scrape_ebi_ids(tokens, proj_id_stems=["PRJEB", "PRJNA", "ERP", "SRP"]):
    """This method allows local files or URLs to be stripped for
    EBI/ENA project IDs using tokenization

    Parameters
    ----------
    tokens: list
        list of parsed tokens from a file or webpage
    proj_id_stems: list
        list of expected project prefixes to look for

    Returns
    ---------
    found_ebi_ids : list
        list of EBI IDs found in the document
    """

    found_ebi_ids = []

    potential_tokens = {
        stem: t for stem in proj_id_stems for t in tokens if stem in t
    }
    potential_tokens = list(potential_tokens.values())

    for t in potential_tokens:
        for stem in proj_id_stems:
            # locate the start of the prefix within the token if present
            start = t.find(stem)  # will return -1 if not found

            if start != -1 and len(t) >= start + len(stem) + 1:
                # if a prefix is found, see if the next character is a number
                next_char = t[start + len(stem)]
                if next_char.isnumeric():
                    test_study = t[
                        start:
                    ]  # for now, assume that the token contains the ID
                    ebi_url = (
                        "https://www.ebi.ac.uk/ena/browser/view/"
                        + test_study
                        + "&display=xml"
                    )
                    request2 = requests.get(ebi_url)
                    invalid = True
                    if request2.status_code == 200:
                        if len(test_study) >= 6:  # min size filter
                            found_study = (
                                test_study.strip()
                                .strip(",")
                                .strip(".")
                                .strip("-")
                            )
                            found_ebi_ids.append(found_study)
                            invalid = False
                    if invalid:
                        logger.warning(
                            "Found stem: "
                            + test_study
                            + " in "
                            + test_study
                            + " but EBI project URL "
                            + ebi_url
                            + " does not exist."
                        )

    return found_ebi_ids


def get_ebi_ids(study_xml_dict):
    """Simple method to retrieve the paired study and project accessions

    If the the EBI accessions are not found, will return a tuple of
    False, False

    Parameters
    -----------
    study_xml_dict: dict
        details retrieved from fetch_ebi_info

    Returns
    ----------
    study_accession, proj_accession: (string, string) or (False, False)
        stud and project accessions if found
    """
    # TODO: does not catch case where STUDY_SET and PROJECT_SET in keys
    # TODO: consider refactor with xml object instead of dict

    study_accession = False
    proj_accession = False

    if "STUDY_SET" in study_xml_dict:
        id_dict = study_xml_dict["STUDY_SET"]["STUDY"]["IDENTIFIERS"]
        study_accession = id_dict["PRIMARY_ID"]
        if "SECONDARY_ID" in id_dict:
            proj_accession = id_dict["SECONDARY_ID"]
            # need to catch for the case when there are two secondary IDs...
            # all project accessions should start with PRJ
            if type(proj_accession) == list:
                proj_accession = "".join(
                    [p for p in proj_accession if "PRJ" in p]
                )
        else:
            logger.warning("No project ID for study: " + study_accession)
            proj_accession = study_accession
    elif "PROJECT_SET" in study_xml_dict:
        id_dict = study_xml_dict["PROJECT_SET"]["PROJECT"]["IDENTIFIERS"]
        proj_accession = id_dict["PRIMARY_ID"]
        if "SECONDARY_ID" in id_dict:
            study_accession = id_dict["SECONDARY_ID"]
        else:
            logger.warning("No study ID for project: " + proj_accession)
    else:
        logger.warning("No study information found.")

    return study_accession, proj_accession


def setup_output_dir(output_dir):
    """Helper function to ensure output path exists and is valid"""

    # output_dir directory
    if output_dir[-1] != "/":
        output_dir += "/"

    if not path.exists(output_dir):
        makedirs(output_dir)

    return output_dir
