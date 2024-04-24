"""sex_check 1.1.0 test suite
"""
#!/usr/bin/env python

import os
import subprocess
import sys
import unittest
from unittest import mock

sys.path.append(os.path.abspath(
    os.path.join(os.path.dirname(__file__), "..")
))

from src.sex_check import (
    check_sex_match,
    get_predicted_sex,
    get_reported_sex,
    get_mapped_reads,
    run_samtools_idxstat
)


class TestRunSamtoolsIdxstat(unittest.TestCase):
    """
    Test case for the run_samtools_idxstat function.
    """

    def setUp(self):
        """
        Set up test data.
        """
        self.bamfile = "test_bamfile.bam"
        self.bamfile_prefix = "test_bamfile"
        self.expected_output_file = self.bamfile_prefix + "_idxstat.tsv"

    def tearDown(self):
        """
        Clean up test data.
        """
        # Remove the expected output file if it exists
        if os.path.exists(self.expected_output_file):
            os.remove(self.expected_output_file)

    @mock.patch("src.sex_check.subprocess.run")
    def test_run_samtools_idxstat_success(self, mock_subprocess_run):
        """
        Test case for successful execution of samtools idxstat.
        """
        # Mock subprocess.run to avoid actual execution
        mock_subprocess_run.return_value.returncode = 0

        # Call the actual function under test
        output_file = run_samtools_idxstat(self.bamfile, self.bamfile_prefix)

        # Assert that subprocess.run is called with the correct arguments
        mock_subprocess_run.assert_called_once_with(
            ['samtools', 'idxstat', self.bamfile],
            stdout=mock.ANY, check=True
        )

        # Assert that the output file name is constructed correctly
        self.assertEqual(output_file, self.expected_output_file)

    @mock.patch("src.sex_check.subprocess.run")
    def test_failed_run(self, mock_subprocess_run):
        """
        Test case for handling failure of samtools idxstat.
        """
        # Mock subprocess.run to simulate a failure
        mock_subprocess_run.side_effect = subprocess.CalledProcessError(
            returncode=1, cmd=['samtools', 'idxstat', self.bamfile],
            stderr="Error: failed to run samtools idxstat"
        )

        # Assert that function raises CalledProcessError
        with self.assertRaises(subprocess.CalledProcessError):
            run_samtools_idxstat(self.bamfile, self.bamfile_prefix)


class TestGetMappedReads(unittest.TestCase):
    """
    Test cases for the get_mapped_reads function.

    For expected scores:
    - When both chr1 and chrY are present, the expected score is calculated as
        -log(chrY/chr1 + epislon).
    - When either chr1 or chrY is zero, the expected score is set to
        20.72326583694641, which is -log(1e-9), considering epsilon.
    """

    def test_correct_data(self):
        """Test case for correct data."""
        filename = os.path.join("test_data", "correct_data.tsv")
        chr_1, chr_y, score = get_mapped_reads(filename)
        self.assertEqual(chr_1, 100)
        self.assertEqual(chr_y, 60)
        self.assertAlmostEqual(score, 0.5108, places=4)

    def test_without_chr1(self):
        """Test case for data without chromosome 1."""
        filename = os.path.join("test_data", "without_chr1.tsv")
        chr_1, chr_y, score = get_mapped_reads(filename)
        self.assertEqual(chr_1, 0)
        self.assertEqual(chr_y, 60)
        self.assertAlmostEqual(score, 20.7233, places=4)

    def test_without_chry(self):
        """Test case for data without chromosome Y."""
        filename = os.path.join("test_data", "without_chry.tsv")
        chr_1, chr_y, score = get_mapped_reads(filename)
        self.assertEqual(chr_1, 100)
        self.assertEqual(chr_y, 0)
        self.assertAlmostEqual(score, 20.7233, places=4)

    def test_with_chr11(self):
        """Test case for data without chr1 but with chr11."""
        filename = os.path.join("test_data", "with_chr11.tsv")
        chr_1, chr_y, score = get_mapped_reads(filename)
        self.assertEqual(chr_1, 0)
        self.assertEqual(chr_y, 60)
        self.assertAlmostEqual(score, 20.7233, places=4)


class TestGetReportedSex(unittest.TestCase):
    """
    Unit tests for the get_reported_sex function.
    """

    def test_valid_reported_sex(self):
        """Test for valid reported sex."""
        # When the reported sex is valid
        sample_name = "X12345-GM1234567-23xxxx4-1234-F-12345678"
        self.assertEqual(get_reported_sex(sample_name), "F")

    def test_undetermined_sex(self):
        """Test case for no reported sex in sample name."""
        # When the sex cannot be determined
        sample_name = "X12345-GM1234567"
        self.assertEqual(get_reported_sex(sample_name), "N")

    def test_invalid_reported_sex(self):
        """Test for invalid reported sex."""
        # When the reported sex is invalid
        sample_name = "X12345-GM1234567-23xxxx4-1234-X-12345678"
        self.assertEqual(get_reported_sex(sample_name), "N")

    def test_short_sample_name(self):
        """Test for short sample name."""
        # When the sample name is too short
        sample_name = "X12345"
        self.assertEqual(get_reported_sex(sample_name), "N")


class TestGetPredictedSex(unittest.TestCase):
    """
    Unit tests for the get_predicted_sex function.
    """

    def test_male_prediction(self):
        """Test for when the score is below the male threshold."""
        # When score is lower than male threshold
        self.assertEqual(get_predicted_sex(0.5, 1.0, 2.0), "M")

    def test_female_prediction(self):
        """Test for when the score is above the female threshold."""
        # When score is higher than female threshold
        self.assertEqual(get_predicted_sex(3.0, 1.0, 2.0), "F")

    def test_unknown_prediction(self):
        """Test for when the score falls between male and female thresholds."""
        # When score is between male and female thresholds
        self.assertEqual(get_predicted_sex(1.5, 1.0, 2.0), "U")

    def test_equal_thresholds(self):
        """Test for when male and female thresholds are equal."""
        # When male and female thresholds are equal
        with self.assertRaises(ValueError):
            get_predicted_sex(1.5, 2.0, 2.0)

    def test_invalid_thresholds(self):
        """Test for when male threshold is greater than female threshold."""
        # When male threshold is greater than female threshold
        with self.assertRaises(ValueError):
            get_predicted_sex(1.5, 2.0, 1.0)


class TestCheckSexMatch(unittest.TestCase):
    """
    Unit test for the function check_sex_match. 
    """
    def test_unknown_reported_sex(self):
        """Test case for when the reported sex is unknown."""
        self.assertEqual(check_sex_match("N", "M"), "NA")
        self.assertEqual(check_sex_match("U", "F"), "NA")

    def test_matching_sex(self):
        """Test case for when the reported sex matches the predicted sex."""
        self.assertEqual(check_sex_match("M", "M"), "True")
        self.assertEqual(check_sex_match("F", "F"), "True")

    def test_non_matching_sex(self):
        """Test case for sex mismatch."""
        self.assertEqual(check_sex_match("M", "F"), "False")
        self.assertEqual(check_sex_match("F", "M"), "False")   


if __name__ == '__main__':
    unittest.main()
