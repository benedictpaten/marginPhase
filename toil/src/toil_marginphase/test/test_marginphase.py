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

    OUT_TOIL_VCF_FORMAT = "{}.merged.full.vcf"
    OUT_EXEC_VCF_FORMAT = "{}.vcf"

    DOCKER_MARGIN_PHASE = 'quay.io/ucsc_cgl/margin_phase:latest'

    @classmethod
    def setUpClass(cls):
        super(MarginPhaseTest, cls).setUpClass()
        logging.basicConfig(level=logging.INFO)

    def setUp(self):
        # paths for marginPhase files
        self.toil_test_directory = os.path.dirname(os.path.abspath(__file__))
        self.margin_phase_root = os.path.abspath(os.path.join(self.toil_test_directory, "../../../.."))
        self.margin_phase_docker = os.path.join(self.margin_phase_root, "toil/docker")
        self.margin_phase_test = os.path.join(self.margin_phase_root, "tests")
        # paths for test run
        self.uuid = str(uuid4())
        self.workdir_root = tempfile.mkdtemp()
        self.workdir = os.path.join(self.workdir_root, self.uuid)
        self.toil_outputdir = os.path.join(self.workdir, "toil_output")
        self.exec_outputdir = os.path.join(self.workdir, "exec_output")
        os.mkdir(self.workdir)
        os.mkdir(self.toil_outputdir)
        os.mkdir(self.exec_outputdir)
        # output paths
        self.toil_full_merged_vcf = None
        self.exec_output_vcf = None



    def tearDown(self):
        print("tearing down {}".format(self.workdir_root))
        shutil.rmtree(self.workdir_root)

    # def test_stub(self):
    #     self._generate_manifest("test")
    #     self._generate_config()
    #     self._make_docker_margin_phase()
    #     print("stub tested")

    def test_chunking(self):
        self._make_docker_margin_phase()
        log.info("Docker image created")
        docker_vcf = self._run_docker_marginPhase()
        log.info("Got docker output VCF '{}'".format(docker_vcf))
        toil_vcf = self._run_toil_marginPhase()
        log.info("Got toil output VCF '{}'".format(toil_vcf))
        shutil.copy(docker_vcf, os.path.join(self.toil_test_directory, "DOCKER.unittest.vcf"))
        shutil.copy(toil_vcf, os.path.join(self.toil_test_directory, "TOIL.unittest.vcf"))

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

    def _run_toil_marginPhase(self):
        #prep
        jobStore = os.path.join(self.workdir, 'toil-jobstore')
        work_dir = os.path.join(self.workdir, 'toil-workdir')

        # run toil
        toil_command = ['toil-marginphase', 'run',
                              '--config', self._generate_config(),
                              '--manifest', self._generate_manifest(self.uuid),
                              '--workDir', work_dir,
                              jobStore]
        log.info('Running %r', toil_command)
        subprocess.check_call(toil_command)

        # validate output
        extract_command = ['tar', 'xvf', "{}.tar.gz".format(self.uuid)]
        subprocess.check_call(extract_command, cwd=self.toil_outputdir)
        output_vcf_name = MarginPhaseTest.OUT_TOIL_VCF_FORMAT.format(self.uuid)
        full_merged_vcf = os.path.join(self.toil_outputdir, output_vcf_name)
        if not os.path.isfile(full_merged_vcf):
            contents = subprocess.check_output(['ls', '-la'], cwd=self.toil_outputdir)
            raise UserError("toil output vcf '{}' not found in directory '{}' with contents:\n{}"
                            .format(output_vcf_name, self.toil_outputdir, contents))

        # save and return
        self.toil_full_merged_vcf = full_merged_vcf
        return full_merged_vcf

    def _run_docker_marginPhase(self):
        # prep
        shutil.copy(os.path.join(self.margin_phase_test, MarginPhaseTest.IN_REF_FA),
                    os.path.join(self.exec_outputdir, MarginPhaseTest.IN_REF_FA))
        shutil.copy(os.path.join(self.margin_phase_test, MarginPhaseTest.IN_REF_VCF),
                    os.path.join(self.exec_outputdir, MarginPhaseTest.IN_REF_VCF))
        shutil.copy(os.path.join(self.margin_phase_test, MarginPhaseTest.IN_BAM),
                    os.path.join(self.exec_outputdir, MarginPhaseTest.IN_BAM))
        shutil.copy(os.path.join(self.toil_test_directory, MarginPhaseTest.IN_PARAMS),
                    os.path.join(self.exec_outputdir, MarginPhaseTest.IN_PARAMS))

        # run docker
        docker_command = ['docker', 'run', '--rm', '-v', "{}:/data".format(self.exec_outputdir), # '-it',
                          MarginPhaseTest.DOCKER_MARGIN_PHASE,
                          "/data/{}".format(MarginPhaseTest.IN_BAM), "/data/{}".format(MarginPhaseTest.IN_REF_FA),
                          "-p", "/data/{}".format(MarginPhaseTest.IN_PARAMS), "-o", "/data/{}".format(self.uuid),
                          # "-r", "/data/{}".format(MarginPhaseTest.IN_REF_VCF),
                          "-a", "info", "-v", "0"]
        log.info('Running %r', docker_command)
        subprocess.check_call(docker_command)

        # validate output
        output_vcf_name = MarginPhaseTest.OUT_EXEC_VCF_FORMAT.format(self.uuid)
        output_vcf = os.path.join(self.exec_outputdir, output_vcf_name)
        if not os.path.isfile(output_vcf):
            contents = subprocess.check_output(['ls', '-la'], cwd=self.exec_outputdir)
            raise UserError("exec output vcf '{}' not found in directory '{}' with contents:\n{}"
                            .format(output_vcf_name, self.exec_outputdir, contents))

        # save and return
        self.exec_output_vcf = output_vcf
        return output_vcf

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
                output_dir=self.toil_outputdir,
                partition_size=partition_size,
                partition_margin=partition_margin,
                reference_fa=os.path.join(self.margin_phase_test, MarginPhaseTest.IN_REF_FA),
                params=os.path.join(self.toil_test_directory, MarginPhaseTest.IN_PARAMS)
            ))
        # with open(path, 'r') as f:
        #     for line in f:
        #         print(line)
        return path

    def _generate_manifest(self, uuid):
        path = os.path.join(self.workdir, 'manifest-toil-rnaseq.tsv')
        with open(path, 'w') as f:
            f.write("{}\t{}".format(uuid, os.path.join(self.margin_phase_test, MarginPhaseTest.IN_BAM)))
        # with open(path, 'r') as f:
        #     for line in f:
        #         print(line)
        return path

