
from __future__ import print_function
import sys
import subprocess

# read.0
# read.1
# read.2
# chr
# subs
# ptype

def assert_exists(loc):
    if not loc.startswith("s3://"): return True
    if len(subprocess.check_output(['s3cmd', 'ls', loc])) == 0:
        print("DNE: {}".format(loc), file=sys.stderr)
        return False
    return True

UUID = "NA12878.hg38.{read0}.{chr}.{sub}.{ptype}"
BAM = ""
CHR = "{chr}"
FA = "s3://margin-phase/fasta/hg38.{chr}.fa"
PARAMS = "s3://margin-phase/params/averaged/params.{read1}.averaged.{sub}.{ptype}.json"
VCF = "s3://margin-phase/vcf/NA12878.hg38.PG.{chr}.vcf"

reads = [
    ['np2', 'np', 's3://margin-phase/bam/nanopore2/NA12878.hg38.np2.{chr}.bam'],
    ['np1-2', 'np', 's3://margin-phase/bam/nanopore1-2/NA12878.hg38.np1-2.mm.{chr}.bam'],
    ['pb', 'pb', 's3://margin-phase/bam/realigned/NA12878.hg38.pb.mm.{chr}.bam']
]
chrs = ['chr21']
subs = ['998', '9995']
ptypes = ['rp', 'rp-train1', 'plain']


have_everything = True
for read in reads:
    for chr in chrs:
        for sub in subs:
            for ptype in ptypes:
                uuid = UUID.format(read0=read[0], chr=chr, sub=sub, ptype=ptype)
                bam = read[2].format(read0=read[0], chr=chr)
                chr = CHR.format(chr=chr)
                fa = FA.format(chr=chr)
                params = PARAMS.format(read1=read[1], sub=sub, ptype=ptype)
                vcf = VCF.format(chr=chr)

                sample = [uuid, bam, chr, fa, params, vcf]
                print("\t".join(sample))
                for s in sample:
                    have_everything = have_everything and assert_exists(s)
                if not have_everything: sys.exit(1)
