"""sex_check 1.1.0 integration test suite

This test validates the functionality of the app as a whole on the DNAnexus 
platform.
"""
#!/usr/bin/env python

from __future__ import print_function

import os
import time
import unittest

import dxpy
import dxpy.app_builder

from dxpy.exceptions import DXAPIError


src_dir = os.path.join(os.path.dirname(__file__), "..")
test_resources_dir = os.path.join(src_dir, "test", "resources")


def make_inputs():
    """basic inputs with real DNAnexus sample to test app.
    Using a control sample randomly selected from
    a CEN run; Thresholds are for testing functionality only.
    """
    return {
        "male_threshold": 4.45,
        "female_threshold": 5.40,
        "input_bam": {
            "$dnanexus_link": {
                "project": "project-Ggyb2G84zJ4363x2JqfGgb6J",
                "id": "file-GgybG1j4VvFvj4zqg20px2gY"
            }
        },
        "index_file": {
            "$dnanexus_link": {
                "project": "project-Ggyb2G84zJ4363x2JqfGgb6J",
                "id": "file-GgybG1j4VvFp7xFf2F6gzQjP"
            }
        }
    }


class TestDXApp(unittest.TestCase):
    """test functionality of app on DNAnexus
    """
    @classmethod
    def setUpClass(cls):
        """Uploads the app to DNAnexus Platform.
        """
        cls.base_input = make_inputs()
        bundled_resources = dxpy.app_builder.upload_resources(src_dir)
        try:
            app_name = os.path.basename(os.path.abspath(src_dir)) + "_test"
        except OSError:
            app_name = "test_app"
        applet_basename = app_name + "_" + str(int(time.time()))
        cls.applet_id, _ignored_applet_spec = dxpy.app_builder.upload_applet(
            src_dir, bundled_resources, override_name=applet_basename
            )

    @classmethod
    def tearDownClass(cls):
        """Clean up by removing the app we created.
        """
        try:
            dxpy.api.container_remove_objects(dxpy.WORKSPACE_ID,
                                              {"objects": [cls.applet_id]})
        except DXAPIError as e:
            print(f"Error removing {cls.applet_id} during cleanup; ignoring.")
            print(e)

    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_base_input(self):
        """
        Tests if app runs successfuly on DNAnexus and returns expected output.
        """
        job = dxpy.DXApplet(self.applet_id).run(self.base_input)
        print(f"Waiting for {job.get_id()} to complete")
        job.wait_on_done()

        output = job.describe()["output"]

        # Assert output contains expected fields defined in dxapp.json
        self.assertIn("idxstat_output", output)
        self.assertIn("sex_check_result", output)

        # Assert values are valid DNAnexus file links
        self.assertTrue(dxpy.is_dxlink(output["idxstat_output"]))
        self.assertTrue(dxpy.is_dxlink(output["sex_check_result"]))


if __name__ == '__main__':
    unittest.main()
