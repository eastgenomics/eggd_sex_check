<!-- dx-header -->
# eggd_sex_check (DNAnexus Platform App)

Verifies the reported sex of given sample

This is the source code for an app that runs on the DNAnexus Platform.
For more information about how to run or modify it, see
https://documentation.dnanexus.com/.
<!-- /dx-header -->

## What does this app do?
This app determines the sex of a sample based on mapped read counts to chromosomes Y in a given BAM file. It employs thresholds (ie number of reads mapped to chromosome Y) to predict the sex and compares it against the reported sex embedded in the sample name.
<br></br>

## What are typical use cases for this app?
The app can be used as one of the quality control steps of sequencing data to confirm the reported sex against the sex inferred from the genomic data.
<br></br>

## What data are required for this app to run?
Required input files:

1. A BAM file (.bam),
2. A corresponding index file (.bai),
3. Male threshold (float) - A threshold for determining male based on chromosome Y reads.
4. Female threshold (float) - A threshold for determining female based on chromosome Y reads.

## What does this app output?
This app outputs two main files:

- {prefix}_idxstat.tsv: A file containing the output of samtools idxstat.
- {prefix}_mqc.json: A MultiQC compatible file summarising the sex check results. This file includes the sample name, mapped reads for chromosomes 1 and Y, normalised score, reported sex, and predicted sex.
<br></br>

## How to run this app from the command
To run this app, you would use a command like:

```
dx run eggd_sex_check \
-iinput_bam=file-xxxx \
-iindex_file=file-yyyy \
-imale_threshold={} \
-ifemale_threshold={}

```
Note: Replace file-xxxx and file-yyyy with the actual file IDs for the BAM file and its index file. Adjust male_threshold and female_threshold according to your requirements.
<br></br>