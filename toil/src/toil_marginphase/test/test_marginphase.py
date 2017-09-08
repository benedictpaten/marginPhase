import logging
import os
import sys
from toil_lib import require, UserError
import shutil
import subprocess
import tempfile
import textwrap
from unittest import TestCase
from uuid import uuid4

log = logging.getLogger(__name__)
#todo make log.* actually append to the console output
log.addHandler(logging.StreamHandler(sys.stdout)) #this doesn't work

# paths for marginPhase files
TOIL_TEST_DIR = os.path.dirname(os.path.abspath(__file__))
MARGIN_PHASE_ROOT = os.path.abspath(os.path.join(TOIL_TEST_DIR, "../../../.."))
MARGIN_PHASE_DOCKER = os.path.join(MARGIN_PHASE_ROOT, "toil/docker")
MARGIN_PHASE_TEST = os.path.join(MARGIN_PHASE_ROOT, "tests")

class MarginPhaseTest(TestCase):
    IN_REF_FA = "hg19.chr3.9mb.fa"
    IN_REF_VCF = "NA12878.PG.chr3.2mb.vcf"
    IN_BAM = "NA12878.pb1.chr3.1mb.bam"
    IN_PARAMS = "params.json"

    OUT_TOIL_VCF_FORMAT = "{}.merged.full.vcf"
    OUT_EXEC_VCF_FORMAT = "{}.vcf"

    DOCKER_MARGIN_PHASE = 'quay.io/ucsc_cgl/margin_phase:latest'

    DEBUG = True

    @classmethod
    def setUpClass(cls):
        super(MarginPhaseTest, cls).setUpClass()
        logging.basicConfig(level=logging.INFO)
        subprocess.check_call(['make'], cwd=MARGIN_PHASE_DOCKER)
        log.info("Docker image created")

    def setUp(self):
        self.uuid = str(uuid4())
        # paths for test run
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
        
    def test_five_percent_margin(self):
        self._run("5", 100000, 5000)

    def test_fifty_percent_margin(self):
        self._run("50", 100000, 50000)
        
    def test_hundred_percent_margin(self):
        self._run("100", 100000, 50000)

    def _run(self, identifier, partition_size=100000, partition_margin=5000):
        identifier = "{}-{}".format(identifier, self.uuid)
        docker_vcf = self._run_docker_marginPhase(identifier)
        log.info(identifier + " Got docker output VCF '{}'".format(docker_vcf))
        toil_vcf = self._run_toil_marginPhase(identifier, partition_size, partition_margin)
        log.info(identifier + " Got toil output VCF '{}'".format(toil_vcf))
        docker_name = "DOCKER.test.{}.vcf".format(identifier)
        toil_name = "TOIL.test.{}.vcf".format(identifier)
        shutil.copy(docker_vcf, os.path.join(self.workdir, docker_name))
        shutil.copy(toil_vcf, os.path.join(self.workdir, toil_name))
        if MarginPhaseTest.DEBUG:
            shutil.copy(docker_vcf, os.path.join(TOIL_TEST_DIR, docker_name))
            shutil.copy(toil_vcf, os.path.join(TOIL_TEST_DIR, toil_name))
        self._compare_output(self.workdir, docker_name, toil_name)

    def _run_toil_marginPhase(self, identifier, partition_size, partition_margin):
        #prep
        jobStore = os.path.join(self.workdir, 'toil-jobstore')
        work_dir = os.path.join(self.workdir, 'toil-workdir')
        os.mkdir(work_dir)

        # run toil
        toil_command = ['toil-marginphase', 'run',
                              '--config', self._generate_config(partition_size, partition_margin),
                              '--manifest', self._generate_manifest(identifier),
                              '--workDir', work_dir,
                              jobStore]
        log.info('Running %r', toil_command)
        subprocess.check_call(toil_command)

        # validate output
        extract_command = ['tar', 'xvf', "{}.tar.gz".format(identifier)]
        subprocess.check_call(extract_command, cwd=self.toil_outputdir)
        output_vcf_name = MarginPhaseTest.OUT_TOIL_VCF_FORMAT.format(identifier)
        full_merged_vcf = os.path.join(self.toil_outputdir, output_vcf_name)
        if not os.path.isfile(full_merged_vcf):
            contents = subprocess.check_output(['ls', '-la'], cwd=self.toil_outputdir)
            raise UserError("toil output vcf '{}' not found in directory '{}' with contents:\n{}"
                            .format(output_vcf_name, self.toil_outputdir, contents))

        # save and return
        self.toil_full_merged_vcf = full_merged_vcf
        return full_merged_vcf

    def _run_docker_marginPhase(self, identifier):
        # prep
        shutil.copy(os.path.join(MARGIN_PHASE_TEST, MarginPhaseTest.IN_REF_FA),
                    os.path.join(self.exec_outputdir, MarginPhaseTest.IN_REF_FA))
        shutil.copy(os.path.join(MARGIN_PHASE_TEST, MarginPhaseTest.IN_REF_VCF),
                    os.path.join(self.exec_outputdir, MarginPhaseTest.IN_REF_VCF))
        shutil.copy(os.path.join(MARGIN_PHASE_TEST, MarginPhaseTest.IN_BAM),
                    os.path.join(self.exec_outputdir, MarginPhaseTest.IN_BAM))
        shutil.copy(os.path.join(TOIL_TEST_DIR, MarginPhaseTest.IN_PARAMS),
                    os.path.join(self.exec_outputdir, MarginPhaseTest.IN_PARAMS))

        # run docker
        docker_command = ['docker', 'run', '--rm', '-v', "{}:/data".format(self.exec_outputdir), # '-it',
                          MarginPhaseTest.DOCKER_MARGIN_PHASE,
                          "/data/{}".format(MarginPhaseTest.IN_BAM), "/data/{}".format(MarginPhaseTest.IN_REF_FA),
                          "-p", "/data/{}".format(MarginPhaseTest.IN_PARAMS), "-o", "/data/{}".format(identifier),
                          # "-r", "/data/{}".format(MarginPhaseTest.IN_REF_VCF),
                          "-a", "info", "-v", "0"]
        log.info('Running %r', docker_command)
        subprocess.check_call(docker_command)

        # validate output
        output_vcf_name = MarginPhaseTest.OUT_EXEC_VCF_FORMAT.format(identifier)
        output_vcf = os.path.join(self.exec_outputdir, output_vcf_name)
        if not os.path.isfile(output_vcf):
            contents = subprocess.check_output(['ls', '-la'], cwd=self.exec_outputdir)
            raise UserError("exec output vcf '{}' not found in directory '{}' with contents:\n{}"
                            .format(output_vcf_name, self.exec_outputdir, contents))

        # save and return
        self.exec_output_vcf = output_vcf
        return output_vcf

    def _compare_output(self, work_dir, reference_vcf_name, evaluated_vcf_name):
        # run the vcfCompare executable
        docker_command = ['docker', 'run', '--rm', '-v', "{}:/data".format(work_dir),
                          '--entrypoint', '/opt/marginPhase/vcfCompare', MarginPhaseTest.DOCKER_MARGIN_PHASE,
                          '-r', "/data/{}".format(reference_vcf_name), '-e', "/data/{}".format(evaluated_vcf_name)]
        log.info('Running %r', docker_command)
        vcf_compare_output = subprocess.check_output(docker_command)

        # log output
        log.info("VCF Comparison Output:\n{}".format(vcf_compare_output))
        
        #todo validate
        return vcf_compare_output

    def _generate_config(self, partition_size, partition_margin):
        path = os.path.join(self.workdir, 'config-toil-rnaseq.yaml')
        with open(path, 'w') as f:
            f.write(textwrap.dedent("""
                    reference-vcf: file://{reference_vcf}
                    output-dir: {output_dir}
                    partition-size: {partition_size}
                    partition-margin: {partition_margin}
                    min-merge-ratio: .8
                    margin-phase-tag: latest
                    default-reference: file://{reference_fa}
                    default-contig: chr3
                    default-params: file://{params}
                    save-intermediate-files: False

                    """[1:]).format(
                reference_vcf=os.path.join(MARGIN_PHASE_TEST, MarginPhaseTest.IN_REF_VCF),
                output_dir=self.toil_outputdir,
                partition_size=partition_size,
                partition_margin=partition_margin,
                reference_fa=os.path.join(MARGIN_PHASE_TEST, MarginPhaseTest.IN_REF_FA),
                params=os.path.join(TOIL_TEST_DIR, MarginPhaseTest.IN_PARAMS)
            ))
        # with open(path, 'r') as f:
        #     for line in f:
        #         print(line)
        return path

    def _generate_manifest(self, uuid):
        path = os.path.join(self.workdir, 'manifest-toil-rnaseq.tsv')
        with open(path, 'w') as f:
            f.write("{}\tfile://{}".format(uuid, os.path.join(MARGIN_PHASE_TEST, MarginPhaseTest.IN_BAM)))
        # with open(path, 'r') as f:
        #     for line in f:
        #         print(line)
        return path

