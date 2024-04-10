"""source code for eggd_sex_check
"""
#!/usr/bin/env python

import os
import subprocess
import math
import shutil
import json
import dxpy


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
    with open(output_file, 'w', encoding="utf-8") as outfile:
        subprocess.run(
            ['samtools', 'idxstat', bamfile],
            stdout=outfile, check=True
            )
    print("samtools finished running successfully")

    return output_file

def get_mapped_reads(filename):
    """
    Reads a file containing idxstat output and extracts the mapped reads for
    chromosomes 1 and Y. Then calculates a normalised score (-log(chrY/chr1))

    Args:
        filename (str): The path to the idxstat output file.

    Returns:
        tuple: (number of mapped reads for chromosome 1,
                number of mapped reads for chromosome Y,
                normalised score)
    """
    chr_1 = chr_y = 0
    epsilon = 1e-9 #small value to avoid log(0)

    with open(filename, encoding="utf-8") as file:
        for line in file:
            if line.startswith("1"):
                chr_1 = int(line.split("\t")[-2])
            elif line.startswith("Y"):
                chr_y = int(line.split("\t")[-2])

    if not chr_1:
        print("No mapped reads for chromosome 1. Using 0 instead.")
    if not chr_y:
        print("No mapped reads for chromosome Y. Using 0 instead.")

    n_chr_y = chr_y / chr_1 if chr_1 != 0 else 0
    score = -math.log(n_chr_y + epsilon)

    return chr_1, chr_y, score

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
        print(f"{sample_name} is too short to determine sex. Returning N")
        return "N"

    sex = parts[-2].upper()
    if sex not in ["M", "F", "U"]:
        print(f"Extracted {sex} from {sample_name} is invalid. Returning N")
        return "N"

    return sex

def get_predicted_sex(score, male_threshold, female_threshold):
    """
    Determines the predicted sex based on score and defined thresholds.
    N/B: Higher score = fewer reads

    Args:
        score (float): The normalised reads count for chromosome Y.
        male_threshold (float): The threshold below which the sample is 
        considered male.
        female_threshold (float): The threshold above which the sample is
        considered female.

    Returns:
        str: The predicted sex ('M' for male, 'F' for female, 'U' for unknown).

    Raises:
        ValueError: If the thresholds are not logically set
        (ie male_threshold should be lower than female_threshold).
    """
    # Validate thresholds
    if male_threshold >= female_threshold:
        raise ValueError("Male threshold must be less than female threshold.")

    # Determine sex based on thresholds
    if score <= male_threshold:
        return "M"
    elif score >= female_threshold:
        return "F"
    else:
        return "U"


def check_sex_match(reported_sex, predicted_sex):
    """
    Checks if the predicted sex matches the reported sex, handling cases where
    reported sex is unknown.

    Args:
        reported_sex (str): The reported sex.
        predicted_sex (str): The predicted sex.

    Returns:
        str: "True" if predicted sex matches reported sex, "False" otherwise.
             "NA" if reported sex is "N" or "U".
    """

    if reported_sex in ["N", "U"]:
        return "NA"

    return str(reported_sex == predicted_sex)



@dxpy.entry_point('main')
def main(input_bam, index_file, male_threshold, female_threshold):
    """
    Main function for the DNAnexus app to perform sex determination
    based on BAM file data.

    Args:
        input_bam (str): The ID of the input BAM file in DNAnexus.
        index_file (str): The ID of the index file associated with the BAM file.
        male_threshold (float): Value below which the sample is considered male.
        female_threshold (float): Value above which the sample is considered female.       
    Returns:
        dict: Dictionary of output file links in DNAnexus.
    """

    inputs = dxpy.download_all_inputs()

    shutil.move(inputs['input_bam_path'][0], os.getcwd())
    shutil.move(inputs['index_file_path'][0], os.getcwd())

    bam_file_name = inputs['input_bam_name'][0]
    bam_file_prefix = inputs['input_bam_prefix'][0].rstrip('_markdup')

    idxstat_output = run_samtools_idxstat(bam_file_name, bam_file_prefix)
    chr_1, chr_y, score = get_mapped_reads(idxstat_output)
    predicted_sex = get_predicted_sex(score, male_threshold, female_threshold)
    reported_sex = get_reported_sex(bam_file_name)
    matched = check_sex_match(reported_sex, predicted_sex)

    # format output to mqc json
    data = {
        bam_file_prefix: {
            "matched": matched,
            "reported_sex": reported_sex,
            "predicted_sex": predicted_sex,
            "score": score,
            "mapped_chrY": chr_y,
            "mapped_chr1": chr_1
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
            "format": "{:.0f}"
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
                "description": "Sex inferred from normalised score",
                "cond_formatting_rules": {
                "warn": [{"s_eq": "U"}]
                }
            },
            "score": {
                "title": "Normalised ChrY Reads",
                "description": "Negative log of mapped_chrY/mapped_chr1",
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

    with open(out_file_name, "w", encoding="utf-8") as file:
        json.dump(multiqc_config, file, indent=2)

    idxstat_output = dxpy.upload_local_file(idxstat_output)
    sex_check_result = dxpy.upload_local_file(out_file_name)

    output = {}
    output["idxstat_output"] = dxpy.dxlink(idxstat_output)
    output["sex_check_result"] = dxpy.dxlink(sex_check_result)

    return output


dxpy.run()
