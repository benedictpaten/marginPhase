import logging
import os
from toil_lib import require, UserError
from toil_marginphase.marginphase_pipeline import main

import shlex
import shutil
import subprocess
import tempfile
import textwrap
from contextlib import closing
from unittest import TestCase
from urlparse import urlparse
from uuid import uuid4

import posixpath

log = logging.getLogger(__name__)


class MarginPhaseTest(TestCase):
    IN_REF_FA = "hg19.chr3.9mb.fa"
    IN_REF_VCF = "NA12878.PG.chr3.2mb.vcf"
    IN_BAM = "NA12878.pb1.chr3.1mb.bam"
    IN_PARAMS = "params.json"

    @classmethod
    def setUpClass(cls):
        super(MarginPhaseTest, cls).setUpClass()
        logging.basicConfig(level=logging.INFO)

    def setUp(self):
        self.toil_test_directory = os.path.dirname(os.path.abspath(__file__))
        self.margin_phase_root = os.path.abspath(os.path.join(self.toil_test_directory, "../../../.."))
        self.margin_phase_docker = os.path.join(self.margin_phase_root, "toil/docker")
        self.margin_phase_test = os.path.join(self.margin_phase_root, "tests")
        self.workdir = tempfile.mkdtemp()
        self.outputdir = os.path.join(self.workdir, "output")
        jobStore = os.path.join(self.workdir, 'jobstore-%s' % uuid4())

    def tearDown(self):
        print("tearing down {}".format(self.workdir))
        shutil.rmtree(self.workdir)

    def test_stub(self):
        self._generate_manifest("test")
        self._generate_config()
        self._make_docker_margin_phase()
        print("stub tested")

    # def test_samples_option(self):
    #     self._run(self.base_command, '--samples', self.sample.geturl())
    #     self._assertOutput()
    #
    # def test_manifest(self):
    #     num_samples = int(os.environ.get('TOIL_SCRIPTS_TEST_NUM_SAMPLES', '1'))
    #     self._run(self.base_command, '--manifest', self._generate_manifest(num_samples))
    #     self._assertOutput(num_samples)
    #
    # def _run(self, *args):
    #     args = list(concat(*args))
    #     log.info('Running %r', args)
    #     subprocess.check_call(args)

    def _make_docker_margin_phase(self):
        subprocess.check_call(['make'], cwd=self.margin_phase_docker)

    def _generate_config(self, partition_size=100000, partition_margin=5000):
        path = os.path.join(self.workdir, 'config-toil-rnaseq.yaml')
        with open(path, 'w') as f:
            f.write(textwrap.dedent("""
                    reference-vcf: file://{reference_vcf}
                    output-dir: {output_dir}
                    partition-size: {partition_size}
                    partition-margin: {partition_margin}
                    min-merge-ratio: .8
                    margin-phase-tag: latest
                    default-reference: {reference_fa}
                    default-contig: chr3
                    default-params: {params}
                    save-intermediate-files: False
                    """[1:]).format(
                reference_vcf=os.path.join(self.margin_phase_test, MarginPhaseTest.IN_REF_VCF),
                output_dir=self.outputdir,
                partition_size=partition_size,
                partition_margin=partition_margin,
                reference_fa=os.path.join(self.margin_phase_test, MarginPhaseTest.IN_REF_FA),
                params=os.path.join(self.toil_test_directory, MarginPhaseTest.IN_PARAMS)
            ))
        with open(path, 'r') as f:
            for line in f:
                print(line)
        return path

    def _generate_manifest(self, uuid):
        path = os.path.join(self.workdir, 'manifest-toil-rnaseq.tsv')
        with open(path, 'w') as f:
            f.write("{}\t{}".format(uuid, os.path.join(self.margin_phase_test, MarginPhaseTest.IN_BAM)))
        with open(path, 'r') as f:
            for line in f:
                print(line)
        return path

