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
TOIL_TEST_STORAGE_DIR = os.path.abspath(os.path.join(TOIL_TEST_DIR, "test_results"))
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
    DOCKER_RTG_TOOLS = 'quay.io/ucsc_cgl/rtg_tools:latest'

    ACCEPTABLE_SENSITIVITY_DIFFERENCE = .01
    ACCEPTABLE_PRECISION_DIFFERENCE = .01

    DEBUG = True

    @classmethod
    def setUpClass(cls):
        super(MarginPhaseTest, cls).setUpClass()
        logging.basicConfig(level=logging.INFO)
        subprocess.check_call(['make'], cwd=MARGIN_PHASE_DOCKER)
        log.info("Docker image created")
        if MarginPhaseTest.DEBUG and not os.path.isdir(TOIL_TEST_STORAGE_DIR):
            os.mkdir(TOIL_TEST_STORAGE_DIR)

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
        self._run("005", 100000, 5000)

    # def test_fifty_percent_margin(self):
    #     self._run("050", 100000, 50000)
    #
    # def test_hundred_percent_margin(self):
    #     self._run("100", 100000, 50000)

    def _run(self, identifier, partition_size=100000, partition_margin=5000):
        identifier = "{}-{}".format(self.uuid, identifier)
        docker_vcf = self._run_docker_marginPhase(identifier)
        log.info(identifier + " Got docker output VCF '{}'".format(docker_vcf))
        toil_vcf = self._run_toil_marginPhase(identifier, partition_size, partition_margin)
        log.info(identifier + " Got toil output VCF '{}'".format(toil_vcf))
        docker_name = "{}.DOCKER.vcf".format(identifier)
        toil_name = "{}.TOIL.vcf".format(identifier)
        shutil.copy(docker_vcf, os.path.join(self.workdir, docker_name))
        shutil.copy(toil_vcf, os.path.join(self.workdir, toil_name))
        if MarginPhaseTest.DEBUG:
            shutil.copy(docker_vcf, os.path.join(TOIL_TEST_STORAGE_DIR, docker_name))
            shutil.copy(toil_vcf, os.path.join(TOIL_TEST_STORAGE_DIR, toil_name))
        self._compare_output(self.workdir, identifier, docker_name, toil_name)

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

    def _compare_output(self, work_dir, identifier, docker_vcf_name, toil_vcf_name):
        # prep - get required files
        shutil.copy(os.path.join(MARGIN_PHASE_TEST, MarginPhaseTest.IN_REF_VCF),
                    os.path.join(work_dir, MarginPhaseTest.IN_REF_VCF))
        shutil.copy(os.path.join(MARGIN_PHASE_TEST, MarginPhaseTest.IN_REF_FA),
                    os.path.join(work_dir, MarginPhaseTest.IN_REF_FA))
        reference_sdf_name = "SDF"

        #bgzip
        vcf_bgzip_command = ['docker', 'run', '--rm', '-v', "{}:/data".format(work_dir),
                            MarginPhaseTest.DOCKER_RTG_TOOLS, "bgzip",
                            "/data/{}".format(docker_vcf_name),
                             "/data/{}".format(toil_vcf_name),
                            "/data/{}".format(MarginPhaseTest.IN_REF_VCF)]
        log.info('Running %r', vcf_bgzip_command)
        subprocess.check_call(vcf_bgzip_command)

        #index
        vcf_index_command = ['docker', 'run', '--rm', '-v', "{}:/data".format(work_dir),
                            MarginPhaseTest.DOCKER_RTG_TOOLS, "index",
                            "/data/{}.gz".format(docker_vcf_name),
                             "/data/{}.gz".format(toil_vcf_name),
                             "/data/{}.gz".format(MarginPhaseTest.IN_REF_VCF)]
        log.info('Running %r', vcf_index_command)
        subprocess.check_call(vcf_index_command)

        #sdf
        ref_sdf_command = ['docker', 'run', '--rm', '-v', "{}:/data".format(work_dir),
                            MarginPhaseTest.DOCKER_RTG_TOOLS, "format",
                            "-o", "/data/{}".format(reference_sdf_name), "/data/{}".format(MarginPhaseTest.IN_REF_FA)]
        log.info('Running %r', ref_sdf_command)
        subprocess.check_call(ref_sdf_command)

        # vcf eval prep
        toil_to_docker_eval_identifier = "{}-vcfeval_t2d".format(identifier)
        toil_to_ref_eval_identifier = "{}-vcfeval_t2r".format(identifier)
        docker_to_ref_eval_identifier = "{}-vcfeval_d2r".format(identifier)
        vcf_eval_base = ['docker', 'run', '--rm', '-v', "{}:/data".format(work_dir),
                         MarginPhaseTest.DOCKER_RTG_TOOLS, "vcfeval", "-t", os.path.join("/data", reference_sdf_name)]

        # EVAL: toil to docker
        vcf_eval_command = list(vcf_eval_base)
        vcf_eval_command.extend(
            ["-o", os.path.join("/data" if MarginPhaseTest.DEBUG else "/tmp", toil_to_docker_eval_identifier),
             "-b", "/data/{}.gz".format(docker_vcf_name), "-c", "/data/{}.gz".format(toil_vcf_name)])
        log.info('Running %r', vcf_eval_command)
        t2d_vcf_eval_output = subprocess.check_output(vcf_eval_command)
        if MarginPhaseTest.DEBUG:
            shutil.copytree(os.path.join(work_dir, toil_to_docker_eval_identifier),
                            os.path.join(TOIL_TEST_STORAGE_DIR, toil_to_docker_eval_identifier))

        # EVAL: toil to reference
        vcf_eval_command = list(vcf_eval_base)
        vcf_eval_command.extend(
            ["-o", os.path.join("/data" if MarginPhaseTest.DEBUG else "/tmp", toil_to_ref_eval_identifier),
             "-b", "/data/{}.gz".format(MarginPhaseTest.IN_REF_VCF), "-c", "/data/{}.gz".format(toil_vcf_name)])
        log.info('Running %r', vcf_eval_command)
        t2r_vcf_eval_output = subprocess.check_output(vcf_eval_command)
        if MarginPhaseTest.DEBUG:
            shutil.copytree(os.path.join(work_dir, toil_to_ref_eval_identifier),
                            os.path.join(TOIL_TEST_STORAGE_DIR, toil_to_ref_eval_identifier))

        # EVAL: docker to reference
        vcf_eval_command = list(vcf_eval_base)
        vcf_eval_command.extend(
            ["-o", os.path.join("/data" if MarginPhaseTest.DEBUG else "/tmp", docker_to_ref_eval_identifier),
             "-b", "/data/{}.gz".format(MarginPhaseTest.IN_REF_VCF), "-c", "/data/{}.gz".format(docker_vcf_name)])
        log.info('Running %r', vcf_eval_command)
        d2r_vcf_eval_output = subprocess.check_output(vcf_eval_command)
        if MarginPhaseTest.DEBUG:
            shutil.copytree(os.path.join(work_dir, docker_to_ref_eval_identifier),
                            os.path.join(TOIL_TEST_STORAGE_DIR, docker_to_ref_eval_identifier))

        # now we analyze docker and toil as compared to the reference
        t2r_vcf_eval = t2d_vcf_eval_output.split("\n")
        d2r_vcf_eval = d2r_vcf_eval_output.split("\n")
        if len(t2r_vcf_eval) < 3 or len(d2r_vcf_eval) < 3:
            raise UserError("Incorrect format for vcf eval output: len {}/{} (expected 3)".format(
                len(t2r_vcf_eval), len(d2r_vcf_eval)))
        header = t2r_vcf_eval[0].split()
        precision_idx = None
        sensitivity_idx = None
        idx = 0
        while idx < len(header):
            if header[idx] == "Precision":
                precision_idx = idx
            if header[idx] == "Sensitivity":
                sensitivity_idx = idx
            idx += 1
        t2r_precision = float(t2r_vcf_eval[2].split()[precision_idx])
        t2r_sensitivity = float(t2r_vcf_eval[2].split()[sensitivity_idx])
        d2r_precision = float(d2r_vcf_eval[2].split()[precision_idx])
        d2r_sensitivity = float(d2r_vcf_eval[2].split()[sensitivity_idx])

        precision_diff = abs(t2r_precision - d2r_precision)
        sensitivity_diff = abs(t2r_sensitivity - d2r_sensitivity)
        if precision_diff > MarginPhaseTest.ACCEPTABLE_PRECISION_DIFFERENCE \
                or sensitivity_diff > MarginPhaseTest.ACCEPTABLE_SENSITIVITY_DIFFERENCE:
            self.fail(("Toil and Docker marginPhase runs have unacceptable difference when compared to the reference:\n"
                      "\tPRECISION  \tToil:%5f\tDocker:%5f\n"
                      "\tSENSITIVITY\tToil:%5f\tDocker:%5f") %
                      (t2r_precision, d2r_precision, t2r_sensitivity, d2r_sensitivity))

        return "\nTOIL to DOCKER:\n{}\nTOIL to REFERENCE:\n{}\nDOCKER to REFERENCE:\n{}".format(
            t2d_vcf_eval_output, t2r_vcf_eval_output, d2r_vcf_eval_output)


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
                    unittest: True

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

