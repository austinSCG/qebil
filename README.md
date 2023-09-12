# Table of Contents  

- [Table of Contents](#table-of-contents)
- [Overview](#overview)
- [Installation](#installation)
- [Basic Usage](#basic-usage)
  * [Example Default EBI search](#example-default-ebi-search)
  * [Example Custom EBI search](#example-custom-ebi-search)
- [Downloading studies](#downloading-studies)
  * [Basic study retrieval](#basic-study-retrieval)
  * [Testing settings and peeking at studies](#testing-settings-and-peeking-at-studies)
  * [Downloading fastq files](#downloading-fastq-files)
  * [Default processing of downloaded data](#default-processing-of-downloaded-data)
    + [Qiita-compatibility processing](#qiita-compatibility-processing)
    + [Compatibility check](#compatibility-check)
    + [Fastq file processing](#fastq-file-processing)
  * [Supplementing EBI metadata](#supplementing-ebi-metadata)
  * [Dealing with host contamination](#dealing-with-host-contamination)
- [Stanalone Operations](#stanalone-operations)
  * [Normalizing metadata](#normalizing-metadata)
  * [Merging metadata files](#merging-metadata-files)
  * [Processing fastq files](#processing-fastq-files)
- [FAQ](#faq)
    + [Now what do I do with these files?](#now-what-do-i-do-with-these-files-)
    + [This seems overwhelming, I just want to be able to work with a public study in Qiita](#this-seems-overwhelming--i-just-want-to-be-able-to-work-with-a-public-study-in-qiita)
    + [But I am qiita-help](#but-i-am-qiita-help)
    + [I don't want to use Qiita, can I still use this?](#i-don-t-want-to-use-qiita--can-i-still-use-this-)
    + [What about NCBI/SRA?](#what-about-ncbi-sra-)
    + [I tried normalizing my metadata, but my sample_type isn't part of the default list](#i-tried-normalizing-my-metadata--but-my-sample-type-isn-t-part-of-the-default-list)
    + [I don't like what it does to my file.](#i-don-t-like-what-it-does-to-my-file)
    + [I want more columns to be automatically validated.](#i-want-more-columns-to-be-automatically-validated)
    + [I want to use the same validation across all sample types.](#i-want-to-use-the-same-validation-across-all-sample-types)
    + [I found a bug.](#i-found-a-bug)
    + [I found a lot of bugs,](#i-found-a-lot-of-bugs-)
    + [I wish QEBIL could do something more, better, or different.](#i-wish-qebil-could-do-something-more--better--or-different)
    + [I need help.](#i-need-help)

<small><i><a href='http://ecotrust-canada.github.io/markdown-toc/'>Table of contents generated with markdown-toc</a></i></small>

# Overview

The Qiita-EBI-Liaison, QEBIL (canonically pronounced, 'kibble'), is a tool for assisting with the movement of sequencing file and their associated metadata into Qiita. In addition, QEBIL allows for the rapid searching and retrieval of metadata and fastq files from the EBI/ENA repository with automatic reformatting into Qiita-compliant files. Two additional features under development are automatic quality-filtering and host-depletion following standard CMI/Knight Lab practices and automatic normalization using Qiimp-formatted templates.

Projects downloaded this way can further be automatically imported into new studies in Qiita with assistance from Qiita Help (qiita.help@gmail.com)

There are three primary uses for QEBIL:
1) Find studies in EBI based on search criteria including terms, scientific names, library preparation strategies (e.
g. AMPLICON, WGS, RNA-Seq),and sequencing platform
2) Downloading metadata (and optionally fastq files) from one or more studies in public repositories
3) Metadata augmentation and normalization using Qiimp-format templates for import into Qiita

Results from (1) can be fed into (2) which will download the identified studies.

QEBIL is best suited for the download of intentially generated metagenomic/metatranscriptomic data, but metagenomic/metatranscriptomic data can also be scraped from studies aimed solely at host sequencing in the manner of Poore et al. (Nature 2020). Support for downloading isolates and clones has also been added, but not extensively tested. Several use cases are illustrated below in addition to an explanation of the many options/parameters that can be used to configure the output to your needs. 

Additional features of QEBIL without reference to EBI to peform common metadata tasks including normalization, validation, merging, etc. as well as to process fastq files for quality filtering or host-depletion are under active development. See the section on [Standalone Operations](#Standalone-Operations) below.

# Installation

Installation assumes you have installed Anaconda or Miniconda. If not, please see the instructions provided on the [QIIME 2 website](https://docs.qiime2.org/2021.4/install/native/#miniconda)

After that, clone the GitHub repository, and pip install:
```
git clone https://github.com/ucsd-cmi/qebil

conda env create -n qebil --file qebil/qebil_env.yml

conda activate qebil

pip install -e qebil
```

If you choose to install in another environment, be sure to run the following command to ensure you have all the non-pip installable dependencies:
```
conda install fastp minimap2 samtools fqtools
```

# Basic Usage

For many users the default settings will be suitable.
```
Usage: qebil [OPTIONS] COMMAND [ARGS]...

Options:
  --version  Show the version and exit.
  --help     Show this message and exit.

Commands:
  fetch     Fetch metadata from EBI/ENA for studies of interest matching...
  metadata  Normalize metadata
  search    Search EBI/ENA for studies of interest matching criteria
```
  
  

## Example Default EBI search

```
Usage: qebil search ebi [OPTIONS]

  Search EBI for studies that match the criteria specified

  This method performs a search in EBI for the provided terms, applies the
  selection criteria as a filter, and then aggregates and summarizing the
  results found.

Options:
  --query TEXT                    EBI search string to identify studies to
                                  download.

  --output-dir TEXT               Directory for output_dir files.  Default is
                                  working directory.

  --quiet / --verbose             Sets logging level for information. If
                                  --quiet (default) logs will only be written
                                  to log file.

  --prefix TEXT                   Prefix to prepend to results and log files.
                                  Defaults to timestamp.

  --source [GENOMIC|GENOMIC SINGLE CELL|TRANSCRIPTOMIC|TRANSCRIPTOMIC SINGLE CELL|METAGENOMIC|METATRANSCRIPTOMIC|SYNTHETIC|VIRAL RNA|OTHER]
                                  Library sources for restricting search.
  --strategy [POOLCLONE|CLONE|CLONEEND|WGS|WGA|WCS|WXS|AMPLICON|ChIP-Seq|RNA-Seq|MRE-Seq|MeDIP-Seq|MBD-Seq|MNase-Seq|DNase-Hypersensitivity|Bisulfite-Seq|EST|FL-cDNA|miRNA-Seq|ncRNA-Seq|FINISHING|TS|Tn-Seq|VALIDATION|FAIRE-seq|SELEX|RIP-Seq|ChIA-PET|RAD-Seq|Other]
                                  Library strategy to restrict search.
  --platform [LS454|Illumina|Ion_Torrent|PacBio_SMRT|OXFORD_NANOPORE]
                                  Instrument platform to restrict search.
  --selection [RANDOM|PCR|RANDOM PCR|RT-PCR|HMPR|MF|repeat fractionation|size fractionation|MSLL|cDNA|cDNA_randomPriming|cDNA_oligo_dT|PolyA|Oligo-dT|Inverse rRNA|Inverse rRNA selection|ChIP|ChIP-Seq|MNase|DNase|Hybrid Selection|Reduced Representation|Restriction Digest|5-methylcytidine antibody|MBD2 protein methyl-CpG binding domain|CAGE|RACE|MDA|padlock probes capture method|unspecified|other]
                                  Library selection method to restrict search.
  --scientific-name TEXT          Scientific names to restrict search.
  --no-filter                     Ignore all defaults and do not filter
                                  samples with any default selection criteria.

  --summarize [center_name|host_scientific_name|instrument_model|instrument_platform|library_layout|library_selection|library_source|library_strategy|scientific_name|study_title|submitted_format|tax_id]
                                  List of metadata to group summary by, in
                                  order.

  --help                          Show this message and exit.
```
Let's do a basic search for COVID19 studies:
``` 
qebil search ebi --query "COVID19" --output-dir example/ --prefix COVID19_defaults
```
This will produce a file with tab-separated results detailing the Qiita-compatible studies found:

[COVID19_defaults_search_results.tsv](examples/COVID19_defaults_search_results.tsv)

Adding further terms can limit the search results:
```
qebil search ebi --query "COVID19" --output-dir example/ --prefix COVID19_gut_metaG --scientific-name "gut metagenome"
```

This file is reduced to those studies with only "gut metagenome" samples.

[COVID19_gut_metaG_search_results.tsv](examples/COVID19_gut_metaG_search_results.tsv)


Multiple criteria for the same and different parameters can be passed:
```
qebil search ebi --query "COVID19" --output-dir example/ --prefix COVID19_human_metaG_WGS --scientific-name "gut metagenome" --scientific-name "human metagenome" --strategy "WGS"
```
   
[COVID19_human_metaG_WGS_search_results.tsv](examples/COVID19_human_metaG_WGS_search_results.tsv)
   
Note: for any term not passed to filter the results, the default terms are used. These include:
--source:
    -genomic
    -genomic single cell
    -transcriptomic
    -transcriptomic single cell
    -metagenomic
    -metatranscriptomic
    -viral rna
    -other
--strategy
    -amplicon
    -other
    -wgs
    -rna-seq
    -wcs
    -poolclone
    -clone
--platform
    -illumina
--selection
    -random
    -pcr
    -random pcr
    -rt-pcr
    -cdna
    -cdna_randompriming
    -inverse rrna
    -inverse rrna selection
    -unspecified
    -size fractionation
    -repeat fractionation
    -race
    -other

None of the parameters are case-sensitive, so Illumina, illumina, and ILLUMINA are all treated the same.

Finally, the metadata categories to summarize by can also be set, provided they are in the metadata. By default the following categories are used for summarizing:
-instrument_platform
-library_selection
-library_source
-library_strategy
-scientific_name
-study_title

And the additional options are:
-center_name
-host_scientific_name
-instrument_model
-library_layout
-submitted_format
-tax_id


## Example Custom EBI search

We can remove the Qiita-default criteria to see more studies with the --no-filter flag:
 ```
qebil search ebi --query "COVID19" --output-dir example/ --prefix COVID19_no_filter --no-filter
 ```
 
[COVID19_no_filter_search_results.tsv](examples/COVID19_no_filter_search_results.tsv)

Multiple terms can be added to further impact the search results, and the summarize list is reduced to only the filtering terms unless otherwise specified:
``` 
qebil search ebi --query "COVID19" --output-dir example/ --prefix COVID19_gut_long_read --no-filter \
--scientific-name "gut metagenome" --scientific-name "human gut metagenome" \
--strategy WGS --strategy RNA-seq \
--platform "OXFORD_NANOPORE"
```

# Downloading studies

```
Usage: qebil fetch project [OPTIONS]

  Retrieves metadata and fastq files if requested, for samples from studies
  specified

  This method retrieves metadata (and data) for a study in EBI/ENA as
  requested either by --ebi-id or listed in the --project-file location,
  output from "qebil search ebi" or manually created. If a metadata file is
  supplied, only the samples in that file will be retreived rather than the
  whole study. Numerous options are provided to enable customization of the
  output, filtering for samples that match certain criteria, and
  augmentation with either additional metadata and/or with EMP protocol
  information.
  
  The script is designed to be robust to restarts and performs multiple
  safety checks to ensure the files downloaded match those in the EBI/ENA
  repository.

  Several additional parameters to handle common situations are available:
  a) --max-samples and --random-subsample allow for users to probe a study
      without accessing all of the meatadata and data
  b) --add-metadata-file allows users to supply additional metadata to be
      merged with the metadata on EBI using --merge-column as the shared
      key. The term "auto" can be supplied to --merge-column to attempt to
      detect automatically
  c) --emp-protocol allows users to add standard EMP protocol metadata fields
      and valuesfor 16S V4 sequencing to preparation info files
  d) --correct-index allows users to automatically resolve studies with three
      reads by renaming and removing the index file (presumed to be the first
      in the fastq_ftp list
  e) Untested/under development:
      On-the-fly (per sample) quality filtering with fastp and human
      read removal and quantification is in place but needs unittesting
      before going live.

Options:
  --ebi-id TEXT                   EBI/ENA project or study accession(s) to
                                  retrieve

  --download-fastq                Whether to download the associated fastq
                                  files.

  --human-removal                 On-the-fly human read removal. See
                                  documentation for settings and information
                                  for other organisms.

  --qiita / --raw                 Whether to format the output metadata as
                                  sample and prep info.

  --overwrite                     Do not check for existing metadata,
                                  overwrite results.

  --correct-index                 Fix the common issue of three reads by
                                  removing the index file.

  --project-file FILE             File with list of EBI/ENA projects or study
                                  accession(s) to retrieve, one per line.

  --publication FILE              Publications (pdf, url, or plain text)
                                  containing list of  NCBI/SRA or EBI/ENA
                                  projects or study accession(s) to retrieve.

  --metadata-file FILE            Metadata file(s) to use for processing
                                  specific samples.

  --prep-max INTEGER              Max number of samples per prep info file.
                                  Default is 250

  --prefix TEXT                   Prefix to prepend to results and log files.
  --output-dir TEXT               Directory for output_dir files.  Default is
                                  working directory.

  --quiet / --verbose             Sets logging level for information. If
                                  --quiet (default) logs will only be written
                                  to log file.

  --cpus INTEGER                  Number of processors to use during host
                                  depletion. Default is 4.

  --keep-files                    Whether or not to retain raw and
                                  intermediate fastq files.

  --source [GENOMIC|GENOMIC SINGLE CELL|TRANSCRIPTOMIC|TRANSCRIPTOMIC SINGLE CELL|METAGENOMIC|METATRANSCRIPTOMIC|SYNTHETIC|VIRAL RNA|OTHER]
                                  Library sources for restricting search.
  --strategy [POOLCLONE|CLONE|CLONEEND|WGS|WGA|WCS|WXS|AMPLICON|ChIP-Seq|RNA-Seq|MRE-Seq|MeDIP-Seq|MBD-Seq|MNase-Seq|DNase-Hypersensitivity|Bisulfite-Seq|EST|FL-cDNA|miRNA-Seq|ncRNA-Seq|FINISHING|TS|Tn-Seq|VALIDATION|FAIRE-seq|SELEX|RIP-Seq|ChIA-PET|RAD-Seq|Other]
                                  Library strategy to restrict search.
  --platform [LS454|Illumina|Ion_Torrent|PacBio_SMRT|OXFORD_NANOPORE]
                                  Instrument platform to restrict search.
  --selection [RANDOM|PCR|RANDOM PCR|RT-PCR|HMPR|MF|repeat fractionation|size fractionation|MSLL|cDNA|cDNA_randomPriming|cDNA_oligo_dT|PolyA|Oligo-dT|Inverse rRNA|Inverse rRNA selection|ChIP|ChIP-Seq|MNase|DNase|Hybrid Selection|Reduced Representation|Restriction Digest|5-methylcytidine antibody|MBD2 protein methyl-CpG binding domain|CAGE|RACE|MDA|padlock probes capture method|unspecified|other]
                                  Library selection method to restrict search.
  --scientific-name TEXT          Scientific names to restrict search.
  --no-filter                     Ignore all defaults and do not filter
                                  samples with any default selection criteria.

  --add-metadata-file FILE        Supply additional metadata file for merging.
                                  Uses sample_name unless otherwise specified
                                  with --merge-column

  --merge-column TEXT             Column for merging supplemental metadata
                                  with downloaded metadata.

  --emp-protocol                  Update prep information with EMP protocol
                                  standards  for 16S rRNA sequencing.

  --max-samples TEXT              Max number of samples to grab from the
                                  study.

  --random-subsample              When sampling, randomly select subset for
                                  processing.N.B. must supply a number with
                                  --max_samples).

  --help                          Show this message and exit.
```
## Basic study retrieval
The basic usage is to provide one or more EBI/ENA project or study IDs either as a list, in a file, or in the form of a paper/publication to then download and format to be ready for Qiita. Supplying the study or project ID will produce the same result. N.B. if no --prefix is provided the study/project ID supplied will be used, so in theory the same study could be accidentally downloaded twice.

Each of these commands would produce the same result:
```
qebil fetch project --ebi-id PRJNA660883 --output-dir example/ --prefix COVID19_gut_metaG
qebil fetch project --ebi-id SRP283872 --output-dir example/ --prefix COVID19_gut_metaG
qebil fetch project --project-file example/COVID19_gut_metaG_EBI_search_results.tsv --output-dir example/ --prefix COVID19_gut_metaG
qebil fetch project --publication https://www.medrxiv.org/content/10.1101/2020.09.03.20183947v2.full-text --output-dir example/ --prefix COVID19_gut_metaG 
```
Any plain text file of IDs (one per line) or tsv with study/project ids in the first column, or in a column labelled 'study_id' can also be read.

## Testing settings and peeking at studies

Sometimes we may just want to test whether we've entered the parameters correctly, or to just peek at the samples in the study to see how good the metadata is, or the percent microbial reads before downloading. In this case we can supply the paramater --max-samples `[N]` which will grab only the first N samples from the study. You may also supply the --random-subsample flag to instead choose a random subset of the samples. For example:
```
qebil fetch project --ebi-id SRP283872 --output-dir example/ --prefix COVID19_gut_metaG --max-samples 2 --random-subsample
```

## Downloading fastq files

If we want to have the fastq files downloaded, we simply add the flag --download-fastq:
```
qebil fetch project --ebi-id SRP283872 --output-dir example/ --prefix COVID19_gut_metaG --download-fastq
```

## Default processing of downloaded data

Following download, the metadata will go through several levels of processing:

### Qiita-compatibility processing
1) columns will be added to translate between EBI and Qiita terms for data types, e.g. WGS -> Metagenomic, RNA-Seq -> Metatranscriptomic. 2) For AMPLICON or Other data, QEBIL will look for target_gene as a metadata field, and then normalize the target_gene to Qiita-accepted options: 16S, 18S, or ITS.
3) If target_gene is not provided for the study, or the strategy is not recognized the samples will be labelled as AMBIGUOUS. If you are certain from an associated publication, that the EMP protocol was followed, you can add the flag --empo-protocol and the appropriate prep_information terms will be added to the metadata in part D below.

### Compatibility check
1) Samples will be checked for compatibility with the supplied scientific name(s), strategies, platform(s), and source(s)
2) Raw sample and preparation information files will be written for the valid and invalid samples identified.

### Fastq file processing
1) If the --download-fastq flag is passed, fastq files will be downloaded. If download of read 1 or read 2 fails for paired-end studies, the other file will be skipped or removed.
2) All downloaded files will have their remote and local md5 checksums compared and then be validated using fqtools


## Supplementing EBI metadata

Sometimes you may have received metadata from the study authors, or downloaded this metadata from the Supplemental Information in the publication. You may provide this metadata for merging with the metadata from the public repository using the --add-metadata-file flag. You may add multiple files and all of them will be merge using sample_name by default or you can specify your own column with --merge-column. You may also specify --merge-column "auto" to allow QEBIL to try to automatically merge the metadata for you.

## Dealing with host contamination
**UNDER DEVELOPMENT**
In some cases, we'll want to download files where we also want to filter out host (typically human) reads. This is recommended when downloading data from sample types other than feces and is vital when accessing data intended to sequence the human rather than microbial genomes.

To minimize the storage impact of so much data which will subsequently be removed, QEBIL is designed to processs reads in singles (or pairs for paired-end data) to perform quality filtering with fastp followed by human read removal with minimap2. For this to work, you must have access to a minimap2 .mmi file with the indexed host genome which can be specified by --host-db. by default, the tool assumes you have access to the barnacle supercomputer at UCSD and points to the shared footpath, however a copy of this file is located here for consistency:

[Human+PhiX.mmi](https://drive.google.com/file/d/1rufLD4jZlmvitsLC7-SpT6YvavwH4a5L/view?usp=sharing)

To perform host depletion, simply add --host-deplete. For an example, lets run it on a small project:
```
qebil fetch project --ebi-id PRJNA593867 --output-dir example/ --human-removal
```
By default, following host deletion, the raw and quality filtered intermediate fastq files are removed. If you do want to keep the intermediate files, add --keep-files

**NOTE: the metadata for human scraped metageonomic/metatranscriptomic data is typically Qiita-non compliant as the scientific_name indicates Homo sapiens rather than human metagenome with Homo sapiens in host_scientific_name. QEBIL's normalization feature is designed
to help fix this**

 
# Stanalone Operations

Not every routine task for feeding the Qiita requires moving data from public repositories. This section deals with these applications.

## Normalizing metadata
**UNDER DEVELOPMENT**
Metadata is a pain, so leave it to Qiimp and QEBIL. Make sure you have a metadata file with a sample_name column and at least a column for sample_type, scientific_name, and host_scientific name. If you don't have a host, used 'not applicable' for host_scientific name. You can supply tab- or comma-separated data.

Then give QEBIL your file:
qebil metadata normalize --metadata-file your_file.tsv

QEBIL will automatically search its in-built set of Qiimp-generated validators and complete a Qiita compliant metadata file. if you have sequencing data associated with your samples, put this in a column labelled run_prefix and enter the Qiita data type (Metagenomic, Metatranscriptomic, 16S, 18S, or ITS) in a column titled prep_file. You may add the suffixes `'_0'`, `'_1'`, etc. to divide your samples into multiple prep files of the same type. 

## Merging metadata files
As above in the section on supplementing EBI metadata, you may provide one or more metadata files for merging. At least one metadata file should be specified with --metadata-file which will be the 'base' for merging the remainder. The rest can be added with the preceding --add-metadata-file flag. You may add multiple files and all of them will be merge using sample_name by default or you can specify your own column with --merge-column. You may also specify --merge-column "auto" to allow QEBIL to try to automatically merge the metadata for you.

## Processing fastq files
If your fastq files aren't coming from EBI, you can specify where the files are located with --fastq-dir . The --quality-filter, --host-deplete, etc. can then be used to process the files in a standalone format.


# FAQ

### Now what do I do with these files?
 
 Create a study in Qiita and upload [following the instructions online](https://qiita.ucsd.edu/static/doc/html/gettingstartedguide/index.html#creating-a-study). 

### This seems overwhelming, I just want to be able to work with a public study in Qiita

Good news! If you email complete the form here:
[EBI/NCBI Study Import Request Form](https://docs.google.com/forms/d/1SIq_JNWai7cZ2wwjD8xZpTab7qBifLKu3TZm2B363CE/edit?ts=5fbe8c0b&gxids=7628#responses)
the study can likely be automatically imported for you. If you don't get a reply after one week, email qiita-help@gmail.com
 
### But I am qiita-help
Great, run the script with the path to the specified files and automatically create the study.
 
### I don't want to use Qiita, can I still use this?
Sure, you're missing out, but adding --raw will omit the Qiita-specific files, but basic metadata compliance checks will still be performed to make files compliant with Qiime2. 

### What about NCBI/SRA?
Qiita interfaces best with EBI/ENA, and most data and metadata in NCBI/SRA is duplicated in EBI/ENA. In the event that it is not, you can try to use pysradb to see if the files can be obtained, but we do not support or endorse that tool. Metadata and/or files from that tool may possibly still be supplied to QEBIL for further processing.

### I tried normalizing my metadata, but my sample_type isn't part of the default list
### I don't like what it does to my file. 
### I want more columns to be automatically validated.

Great make it a [Qiimp template](https://qiita.ucsd.edu/iframe/?iframe=qiimp) for your sample type and/or additional columns. If you don't like the presets, choose the highest-level custom template 'base other' and add your own columns.

Now specify your new validation files with --validator. You may pass as many validators as you'd like by repeating this flag to create a library of options or you may simply pass --validators-dir to have QEBIL read all .xlsx or .yml files in that folder.

### I want to use the same validation across all sample types.

Add the --global-validator file. An example is provided in support_files/validator/default_mapping.yml

### I found a bug.
### I found a lot of bugs,
### I wish QEBIL could do something more, better, or different.
### I need help.

QEBIL is still in active development. Please email qiita-help@gmail.com with any questions or visit the Qiita Forum.
