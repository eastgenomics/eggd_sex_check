#!/usr/bin/env python

import os
import subprocess
import logging
import shutil
import json
import dxpy

logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')


def run_samtools_idxstat(bamfile, bamfile_prefix):
    """
    Run samtools idxstat on a BAM file and write the output to a TSV file.

    Args:
        bamfile (str): local name of the BAM file.
        bamfile_prefix (str): Prefix for the output file name.

    Returns:
        str: The name of the output file containing the idxstat results.
    """
    # Define the output file name
    output_file = bamfile_prefix + '_idxstat.tsv'

    # Run samtools idxstat and capture any errors
    with open(output_file, 'w') as outfile:
        subprocess.run(['samtools', 'idxstat', bamfile],
                        stdout=outfile, check=True)
        logging.info("idxstat finished running")
    return output_file

def get_mapped_reads(filename):
    """
    Reads a file containing idxstat output and extracts the mapped reads for
    chromosomes 1 and Y.

    Args:
        filename (str): The path to the idxstat output file.

    Returns:
        tuple: A tuple containing the number of mapped reads for chromosome 1,
        chromosome Y, and the normalized chromosome Y mapped reads 
        (normalized to chromosome 1 reads).
    """
    chr1_mapped_reads = chrY_mapped_reads = 0

    with open(filename) as file:
        found_chr1 = found_chrY = False
        for line in file:
            if line.startswith("1"):
                chr1_mapped_reads = int(line.split("\t")[-2])
                found_chr1 = True
            elif line.startswith("Y"):
                chrY_mapped_reads = int(line.split("\t")[-2])
                found_chrY = True

        if not (found_chr1 and found_chrY):
            logging.warning("File does not contain mapped reads for chromosome \
                1 and/or Y. Using 0 instead")

    n_chrY = chrY_mapped_reads / chr1_mapped_reads if chr1_mapped_reads != 0 else 0
    return chr1_mapped_reads, chrY_mapped_reads, n_chrY

def get_reported_sex(sample_name):
    """
    Extracts the reported sex from a sample name based on its naming convention.
    Returns 'N' if the sex cannot be determined or is not 'M', 'F', or 'U'.

    Args:
        sample_name (str): The name of the sample, expected to contain the 
        sex information.

    Returns:
        str: The reported sex extracted from the sample name or 'N' if 
        undetermined or invalid.
    """
    parts = sample_name.split('-')
    if len(parts) < 3:
        logging.warning("Sample name '%s' is too short to determine sex. \
                        Returning 'N'.", sample_name)
        return "N"

    sex = parts[-2].upper()
    if sex not in ["M", "F", "U"]:
        logging.warning("Extracted sex '%s' from sample name '%s' is not valid. \
                        Returning 'N'.", sex, sample_name)
        return "N"

    return sex

def get_predicted_sex(nChrY, male_threshold, female_threshold):
    """
    Determines the predicted sex based on nChrY and defined thresholds.

    Args:
        nChrY (float): The normalised reads count for chromosome Y.
        male_threshold (float): The threshold count above which the sample is 
        considered male.
        female_threshold (float): The threshold count below which the sample is
        considered female.

    Returns:
        str: The predicted sex ('M' for male, 'F' for female, 'U' for undetermined).

    Raises:
        ValueError: If the thresholds are not logically set
        (ie male_threshold should be higher than female_threshold).
    """
    # Validate thresholds
    if male_threshold <= female_threshold:
        raise ValueError("Male threshold must be greater than female threshold.")

    # Determine sex based on thresholds
    if nChrY >= male_threshold:
        return "M"
    elif nChrY <= female_threshold:
        return "F"
    else:
        return "U"


@dxpy.entry_point('main')
def main(input_bam, index_file, male_threshold, female_threshold):
    """
    Main function for the DNAnexus app to perform sex determination based on BAM file data.

    Args:
        input_bam (str): The ID of the input BAM file in DNAnexus.
        index_file (str): The ID of the index file associated with the BAM file.
        male_threshold (int): Threshold for determining male based on chromosome Y reads.
        female_threshold (int): Threshold for determining female based on chromosome Y reads.
        
    Returns:
        dict: Dictionary of output file links in DNAnexus.
    """

    inputs = dxpy.download_all_inputs()

    shutil.move(inputs['input_bam_path'][0], os.getcwd())
    shutil.move(inputs['index_file_path'][0], os.getcwd())

    bam_file_name = inputs['input_bam_name'][0]
    bam_file_prefix = inputs['input_bam_prefix'][0].rstrip('_markdup')

    idxstat_output = run_samtools_idxstat(bam_file_name, bam_file_prefix)
    chr1, chrY, nChrY = get_mapped_reads(idxstat_output)
    predicted_sex = get_predicted_sex(nChrY, male_threshold, female_threshold)
    reported_sex = get_reported_sex(bam_file_name)
    matched = str(reported_sex==predicted_sex) if reported_sex != "N" else "NA"

    # format output to mqc json
    data = {
        bam_file_prefix: {
            "matched": matched,
            "reported_sex": reported_sex,
            "predicted_sex": predicted_sex,
            "nChrY": nChrY,
            "mapped_chrY": chrY,
            "mapped_chr1": chr1
        }
    }

    multiqc_config = {
        "id": "sex_check",
        "section_name": "Sex Check",
        "description": "Table comparing reported and predicted sex",
        "plot_type": "table",
        "pconfig": {
            "id": "sex_check_table",
            "title": "Sex Check Table",
            "format": "{:.0f}",
            "min": 0
        },
        "headers": {
            "matched": {
                "title": "Matched",
                "description": "Whether reported sex is same as predicted sex",
                "cond_formatting_rules": {
                "pass": [{"s_eq": "True"}],
                "warn": [{"s_eq": "NA"}],
                "fail": [{"s_eq": "False"}]
                }
            },
            "reported_sex": {
                "title": "Reported Sex",
                "description": "Expected sex reported in sample name",
                "cond_formatting_rules": {
                "warn": [{"s_eq": "N"}, {"s_eq": "U"}]
                }
            },
            "predicted_sex": {
                "title": "Predicted Sex",
                "description": "Sex inferred from mapped_chrY",
                "cond_formatting_rules": {
                "warn": [{"s_eq": "U"}]
                }
            },
            "nChrY": {
                "title": "Normalized ChrY Reads",
                "description": "Ratio of mapped_chrY:mapped_chr1",
                "format": "{:.4f}"
            },
            "mapped_chrY": {
                "title": "Mapped Reads ChrY",
                "description": "Number of reads mapped to chromosome Y"
            },
            "mapped_chr1": {
                "title": "Mapped Reads Chr1",
                "description": "Number of reads mapped to chromosome 1"
            }
        },
        "data": data
    }

    out_file_name = bam_file_prefix + '_mqc.json'

    with open(out_file_name, "w") as file:
        json.dump(multiqc_config, file, indent=2)

    idxstat_output = dxpy.upload_local_file(idxstat_output)
    sex_check_result = dxpy.upload_local_file(out_file_name)

    output = {}
    output["idxstat_output"] = dxpy.dxlink(idxstat_output)
    output["sex_check_result"] = dxpy.dxlink(sex_check_result)

    return output


dxpy.run()
