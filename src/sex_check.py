#!/usr/bin/env python
# sex_check 0.0.1


import os
import dxpy
import subprocess
#import pandas as pd
import logging
import shutil

logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')

def run_samtools_idxstat(bamfile, bamfile_prefix):
    """
    Run samtools idxstat on a BAM file and write the output to a TSV file.

    Args:
        bamfile (str): local name of the BAM file.
        bamfile_prefix (str): Prefix for the output file name.

    Returns:
        str: The name of the output file containing the idxstat results.
    
    Raises:
        subprocess.CalledProcessError: If an error occurs in the samtools idxstat process.
        Exception: For any other unexpected errors.
    """
    try:
        # Define the output file name
        output_file = bamfile_prefix + '_idxstat.tsv'

        # Run samtools idxstat and capture any errors
        with open(output_file, 'w') as outfile:
            subprocess.run(['samtools', 'idxstat', bamfile], stdout=outfile, check=True)
            logging.info("idxstat finished running")
        return output_file

    except subprocess.CalledProcessError as e:
        logging.error("An error occurred while running samtools idxstat: %s", e)
        raise 
    except Exception as e:
        logging.error("An unexpected error occurred: %s", e)
        raise



def get_mapped_reads(filename):
    """
    Reads a file containing idxstat output and extracts the mapped reads for chromosomes 1 and Y.

    Args:
        filename (str): The path to the idxstat output file.

    Returns:
        tuple: A tuple containing the number of mapped reads for chromosome 1, chromosome Y, 
               and the normalized chromosome Y mapped reads (normalized to chromosome 1 reads).

    Raises:
        FileNotFoundError: If the specified file does not exist.
        ValueError: If the file format is incorrect or missing required data.
    """
    try:
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
                raise ValueError("File does not contain required chromosome data.")

        normalized_chrY = chrY_mapped_reads / chr1_mapped_reads if chr1_mapped_reads != 0 else 0
        return chr1_mapped_reads, chrY_mapped_reads, normalized_chrY

    except FileNotFoundError as e:
        logging.error("File not found: %s", e)
        raise
    except ValueError as e:
        logging.error("Data error in file: %s", e)
        raise
    except Exception as e:
        logging.error("An unexpected error occurred: %s", e)
        raise


def get_reported_sex(sample_name):
    """
    Extracts the reported sex from a sample name based on its naming convention.
    Returns 'N' if the sex cannot be determined or is not 'M', 'F', or 'U'.

    Args:
        sample_name (str): The name of the sample, expected to contain the sex information.

    Returns:
        str: The reported sex extracted from the sample name or 'N' if undetermined or invalid.
    """
    try:
        parts = sample_name.split('-')
        if len(parts) < 3:
            logging.warning("Sample name '%s' is too short to determine sex. Returning 'N'.", sample_name)
            return "N"
        
        sex = parts[-2].upper()
        if sex not in ["M", "F", "U"]:
            logging.warning("Extracted sex '%s' from sample name '%s' is not valid. Returning 'N'.", sex, sample_name)
            return "N"

        return sex

    except Exception as e:
        logging.error("An unexpected error occurred while extracting sex from sample name '%s': %s", sample_name, e)
        return "N"



def get_predicted_sex(chrY, male_threshold, female_threshold):
    """
    Determines the predicted sex based on the count of chromosome Y reads and defined thresholds.

    Args:
        chrY (float): The count of mapped reads for chromosome Y.
        male_threshold (float): The threshold count above which the sample is considered male.
        female_threshold (float): The threshold count below which the sample is considered female.

    Returns:
        str: The predicted sex ('M' for male, 'F' for female, 'U' for undetermined).

    Raises:
        ValueError: If the thresholds are not logically set (male_threshold should be higher than female_threshold).
    """
    # Validate thresholds
    if male_threshold <= female_threshold:
        raise ValueError("Male threshold must be greater than female threshold.")

    # Determine sex based on thresholds
    if chrY >= male_threshold:
        return "M"
    elif chrY <= female_threshold:
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

    # Fill in your application code here.
    idxstat_output = run_samtools_idxstat(bam_file_name, bam_file_prefix)
    chr1, chrY, nChrY = get_mapped_reads(idxstat_output)
    predicted_sex = get_predicted_sex(chrY, male_threshold, female_threshold)
    reported_sex = get_reported_sex(bam_file_name)
    
    out_file_name = bam_file_prefix + '_mcq.txt'

    with open(out_file_name, "w") as file:
        file.write("# plot_type: 'table'\n")
        file.write("# section_name: 'Sex Check Results'\n")
        file.write("Sample\tMapped Reads Chr1\tMapped Reads ChrY\tNormalized ChrY Reads\tReported Sex\tPredicted Sex\n")
        file.write(f"{bam_file_prefix}\t{chr1}\t{chrY}\t{nChrY}\t{reported_sex}\t{predicted_sex}\n")


    idxstat_output = dxpy.upload_local_file(idxstat_output)
    sex_check_result = dxpy.upload_local_file(out_file_name)

    output = {}
    output["idxstat_output"] = dxpy.dxlink(idxstat_output)
    output["sex_check_result"] = dxpy.dxlink(sex_check_result)

    return output

dxpy.run()
