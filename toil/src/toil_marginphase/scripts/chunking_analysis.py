#!/usr/bin/env python
from __future__ import print_function
import argparse
import glob
import gzip
import math
import numpy as np
import os

def parse_args():
    parser = argparse.ArgumentParser("Analyzes chunks from curated marginPhase output and logfiles")
    parser.add_argument('--log_dir', '-l', dest='log_dir', required=True, type=str,
                       help='Location where logs are')
    parser.add_argument('--full_vcf_dir', '-f', dest='full_vcf_dir', required=True, type=str,
                       help='Location where full vcfs are')
    parser.add_argument('--merged_vcf_dir', '-m', dest='merged_vcf_dir', required=True, type=str,
                       help='Location where merged vcfs are')

    return parser.parse_args()



def main():
    args = parse_args()
    assert False not in map(os.path.isdir, [args.log_dir, args.full_vcf_dir, args.merged_vcf_dir])

    uuids = map(lambda x: x.rstrip('merged.full.vcf'), os.listdir(args.full_vcf_dir))
    uuids.sort()

    for uuid in uuids:
        print(uuid)
        assert ("np" in uuid) != ("pb" in uuid)
        nanopore = "np" in uuid
        chunk_type = uuid.split(".")[-1]

        full_vcf = os.path.join(args.full_vcf_dir, "{}.merged.full.vcf".format(uuid))
        merged_vcf_files = glob.glob(os.path.join(args.merged_vcf_dir, "{}/{}.merged.*.vcf"
                                                  .format("np" if nanopore else "pb", uuid)))
        log_file = os.path.join(args.log_dir, "merge_chunks/{}.toil-marginPhase.19q.{}.merge_chunks.log"
                             .format("np" if nanopore else "pb", chunk_type))

        assert (False not in map(os.path.isfile, [full_vcf, log_file])) and len(merged_vcf_files) != 0






if __name__ == "__main__":
    main()