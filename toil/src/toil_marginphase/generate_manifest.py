
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

UUID = "NA12878.hg38.{read0}.{chr}.{sub}.{read_err}.{ptype}"
BAM = ""
CHR = "{chr}"
FA = "s3://margin-phase/fasta/hg38.{chr}.fa"
# PARAMS = "s3://margin-phase/params/averaged/params.{read1}.averaged.{sub}.{ptype}.json"
PARAMS = "s3://margin-phase/params/q30-test/params.{sub}.{read_err}.gap.{ptype}.json"
VCF = "s3://margin-phase/vcf/NA12878.hg38.PG.{chr}.vcf"

reads = [
    ['np2', 'np', 's3://margin-phase/bam/nanopore2.q30/NA12878.hg38.np2.q30.{chr}.bam'],
    ['pb', 'pb', 's3://margin-phase/bam/pacbio.q30/NA12878.hg38.pb.mm.q30.{chr}.bam']
]

# chrs = ["chr{}".format(x) for x in range(1,23)]
# chrs.append("chrX")
# chrs.append("chrY")
chrs = ['chr18']

# subs = ['998', '9995']
# subs = ['9995']
subs = ['hs993']

read_err = ['re98', 're998']

# ptypes = ['rp', 'plain']
# ptypes = ['rp', 'rp-train1', 'plain']
# ptypes = ['rp']
# ptypes = ['plain']
ptypes = ['t1','t2','t4']

# have_everything = True
# for read in reads:
#     for chr in chrs:
#         for sub in subs:
#             for ptype in ptypes:
#                 uuid = UUID.format(read0=read[0], chr=chr, sub=sub, ptype=ptype)
#                 bam = read[2].format(read0=read[0], chr=chr)
#                 chr = CHR.format(chr=chr)
#                 fa = FA.format(chr=chr)
#                 params = PARAMS.format(read1=read[1], sub=sub, ptype=ptype)
#                 vcf = VCF.format(chr=chr)
#
#                 sample = [uuid, bam, chr, fa, params, vcf]
#                 print("\t".join(sample))
#                 for s in sample:
#                     have_everything = have_everything and assert_exists(s)
#                 if not have_everything: sys.exit(1)



read_errs = ['re98', 're998']

have_everything = True
unique_uuids = set()
for read in reads:
    for chr in chrs:
        for sub in subs:
            for read_err in read_errs:
                for ptype in ptypes:
                    uuid = UUID.format(read0=read[0], chr=chr, sub=sub, ptype=ptype, read_err=read_err)
                    bam = read[2].format(read0=read[0], chr=chr)
                    chr = CHR.format(chr=chr)
                    fa = FA.format(chr=chr)
                    params = PARAMS.format(read1=read[1], sub=sub, read_err=read_err, ptype=ptype)
                    vcf = VCF.format(chr=chr)

                    sample = [uuid, bam, chr, fa, params, vcf]
                    print("\t".join(sample))
                    if uuid in unique_uuids:
                        print("Duplicated UUID: {}".format(uuid), file=sys.stderr)
                    unique_uuids.add(uuid)
                    for s in sample:
                        have_everything = have_everything and assert_exists(s)
                    if not have_everything: sys.exit(1)