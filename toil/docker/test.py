#!/usr/bin/env python2.7
# John Vivian
import subprocess
import tempfile
import unittest


class TestKallisto(unittest.TestCase):

    def test_docker_call(self):
        out, err = check_docker_output(tool='quay.io/ucsc_cgl/kallisto_sc')
        # this is part of an error message complaining about no 'config' file
        self.assertTrue('"/opt/single_cell/source/10xDetect_and_Prep.py"' in out)

def check_docker_output(tool):
    command = 'docker run ' + tool
    process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    output = process.communicate()
    return output

if __name__ == '__main__':
    unittest.main()