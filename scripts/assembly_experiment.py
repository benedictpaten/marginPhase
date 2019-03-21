#!/usr/bin/env python2
from __future__ import print_function
import argparse
import subprocess
import sys
import os
import numpy as np
import pysam
import datetime

percent = lambda small, big: int(100.0 * small / big) if big != 0 else 0.0


def parse_args(args = None):
    parser = argparse.ArgumentParser("Runs marginPolish on assembly, compares polished assembly to true reference")
    parser.add_argument('--true_reference', '-t', dest='true_reference', default=None, required=True, type=str,
                        help='True reference, polished assembly will be compared to')
    parser.add_argument('--assembly_fasta', '-a', dest='assembly_fasta', default=None, required=False, type=str,
                        help='FASTA file, assembled reference')
    parser.add_argument('--assembly_bam', '-b', dest='assembly_bam', default=None, required=False, type=str,
                        help='BAM file, reads aligned to assembled reference')
    parser.add_argument('--polisher_params', '-p', dest='polisher_params', default=None, required=False, type=str,
                        help='Parameters for marginPolish')

    parser.add_argument('--output_name', '-o', dest='output_name', default="out", required=False, type=str,
                        help='Name to use for output from marginPolish. ')
    parser.add_argument('--read_size', '-r', dest='read_size', default=10000, required=False, type=int,
                        help='Size to chunk reads into for post-polishing alignment')
    parser.add_argument('--force_realign', '-R', dest='force_realign', default=False, required=False, action='store_true',
                        help='Will force overwrite the post-polishing reads and alignment')

    parser.add_argument('--marginPolish_invocation', '-P', dest='marginPolish_invocation', default="./marginPolish",
                        required=False, type=str,
                        help='Invocation for marginPolish.  Default: ./marginPolish (ie: assumes in cwd)')
    parser.add_argument('--minimap2_invocation', '-M', dest='minimap2_invocation', default="minimap2", required=False, type=str,
                        help='Invocation for minimap2.  Default: minimap2 (ie: assumes in path)')

    # parser.add_argument('--verbose', '-v', dest='verbose', action='store_true', default=False,
    #                     help='Print histograms for length and depth')

    return parser.parse_args() if args is None else parser.parse_args(args)


def log(msg):
    print(msg, file=sys.stderr)

# I'm a bastard
global_outfiles = None

def log_result(msg, outfile=None):
    print(msg)
    global global_outfiles
    if outfile is not None:
        with open(outfile, 'a') as outstream:
            print(msg, file=outstream)
    elif global_outfiles is not None:
        for go in global_outfiles:
            with open(go, 'a') as outstream:
                print(msg, file=outstream)


def run_command(command, out_file=None, log_prefix=None, include_pmp=True):

    # include poor man's profiler?
    if include_pmp:
        cmd = ["/usr/bin/time", "-f", "\nDEBUG_MAX_MEM:%M\nDEBUG_RUNTIME:%E\n"]
    else:
        cmd = []

    # finish command
    cmd.extend(command)
    log("Running command:\n\t{}".format(" ".join(list(map(lambda x: x.replace("\n", "\\n"), cmd)))))

    # handle streams or logs
    out_stream = sys.stderr
    err_stream = sys.stderr
    try:
        if out_file is not None:
            out_stream = open(out_file, 'w')
        if log_prefix is not None:
            if out_file is None:
                out_stream = open("{}.stdout.log".format(log_prefix), 'w')
            err_stream = open("{}.stderr.log".format(log_prefix), 'w')
        subprocess.check_call(cmd, stdout=out_stream, stderr=err_stream)
    finally:
        if out_file is not None or log_prefix is not None:
            out_stream.close()
        if log_prefix is not None:
            err_stream.close()


def run_polisher(assembly_fasta, assembly_bam, polisher_params, output_name, margin_polish_invocation):
    # prep
    polish_cmd = [margin_polish_invocation, assembly_bam, assembly_fasta, polisher_params, '-o', output_name]

    # run
    run_command(polish_cmd, log_prefix="{}.marginPolish".format(output_name))


def run_polished_output_to_reads(polished_reference, polished_reads, read_size, quality_char='?'):
    # loggit
    log("Creating reads of size {} from {} to {}".format(read_size, polished_reference, polished_reads))

    # generate reads from polished output
    linenr = -1
    with open(polished_reference, 'r') as input, open(polished_reads, 'w') as output:
        # how to save a read
        def write_read(contig, idx, start_pos, read_lines):
            read = "".join(read_lines)
            read_len = len(read)
            output.write("@{}.i{}.s{}.e{}\n".format(contig, idx, start_pos, start_pos + read_len))
            output.write("{}\n".format(read))
            output.write("+\n")
            output.write("{}\n".format(quality_char * read_len))
            return read_len
        contig_count = 0
        read_count = 0

        current_contig = None
        total_read_len = 0
        current_read_idx = 0
        current_read_len = 0
        current_read_lines = list()
        for line in input:
            # prep
            linenr += 1
            line = line.strip()
            if len(line) == 0:
                log("Found empty line ({}) in {}".format(linenr, polished_reference))
                continue

            # new assembled contig
            if line.startswith(">"):
                contig_count += 1
                if current_contig is not None:
                    write_read(current_contig, current_read_idx, total_read_len, current_read_lines)
                    read_count += 1
                current_contig = line[1:]
                total_read_len = 0
                current_read_idx = 0
                current_read_len = 0
                current_read_lines = list()

            # new read line
            else:
                current_read_lines.append(line)
                current_read_len += len(line)
                if current_read_len >= read_size:
                    total_read_len += write_read(current_contig, current_read_idx, total_read_len, current_read_lines)
                    read_count += 1
                    current_read_idx += 1
                    current_read_lines = list()
                    current_read_len = 0

        # write last read
        if current_contig is not None:
            write_read(current_contig, current_read_idx, total_read_len, current_read_lines)
        else:
            log("ERROR: found no reads in {}".format(polished_reference))

        # loggit
        log("Read {} lines, wrote {} reads out of {} contigs".format(linenr, read_count, contig_count))


def run_polished_read_alignment(true_reference, polished_reads, polished_alignment, minimap2_invocation):
    # get minimap command
    mm2_cmd = [minimap2_invocation, '-a', '-Y', '--eqx', '-x', 'asm5', true_reference, polished_reads]

    # runnit
    run_command(mm2_cmd, out_file=polished_alignment, log_prefix="{}.minimap2".format(
        polished_alignment.replace(".reads.bam", "")))


def run_alignment_analysis(true_reference, polished_alignment, args, outfile=None, type_identifier=""):
    # get contig map
    reference = dict()
    reference_lengths = dict()
    reference_lines = 0
    def save_contig(name, lines):
        reference[name] = "".join(lines)
        reference_lengths[name] = len(reference[name])
    log("Reading true reference: {}".format(true_reference))
    with open(true_reference, 'r') as input:
        current_contig = None
        current_contig_lines = list()
        for line in input:
            reference_lines += 1
            line = line.strip()
            if line.startswith(">"):
                if current_contig is not None:
                    save_contig(current_contig, current_contig_lines)
                current_contig = line[1:].split()[0]
                current_contig_lines = list()
            else:
                current_contig_lines.append(line)
        if current_contig is not None:
            save_contig(current_contig, current_contig_lines)
        else:
            log("ERROR: no contigs found in true reference: {}".format(true_reference))
    total_reference_bases = sum(reference_lengths.values())
    log("Read {} lines from {}, with {} contigs and {} total nucleotides".format(
        reference_lines, true_reference, len(reference), total_reference_bases))

    # read reads
    log("Reading reads from {}".format(polished_alignment))
    alignments = dict()
    all_alignment_summaries = {'=':0,'X':0,'I':0,'D':0,'S':0}
    untracked_alignment_operations = 0
    total_reads = 0
    unaligned_reads = 0
    bamfile = None
    try:
        bamfile = pysam.AlignmentFile(polished_alignment, 'rb' if polished_alignment.endswith("bam") else 'r')
        for read in bamfile.fetch():
            total_reads += 1

            contig = read.reference_name
            aln_start = read.reference_start
            aln_end = read.reference_end
            if aln_end is None or aln_start == -1:
                unaligned_reads += 1
                continue

            if contig not in alignments:
                alignments[contig] = list()
            alignments[contig].append([aln_start, aln_end])

            cigar_summary, _ = read.get_cigar_stats()
            cigar_total = sum(map(lambda x: cigar_summary[x], range(10))) # no nm
            considered_total = sum(map(lambda x: cigar_summary[x], [1,2,4,7,8]))
            untracked_alignment_operations += (cigar_total - considered_total)
            all_alignment_summaries['='] += cigar_summary[7]
            all_alignment_summaries['X'] += cigar_summary[8]
            all_alignment_summaries['I'] += cigar_summary[1]
            all_alignment_summaries['D'] += cigar_summary[2]
            all_alignment_summaries['S'] += cigar_summary[4]

    finally:
        if bamfile is not None: bamfile.close()

    total_summaries = sum(all_alignment_summaries.values())
    total_match = all_alignment_summaries['=']
    total_mismatch = all_alignment_summaries['X']
    total_insert = all_alignment_summaries['I']
    total_delete = all_alignment_summaries['D']
    total_softclip = all_alignment_summaries['S']
    total_aln_nonsc = total_summaries - total_softclip

    global global_outfiles
    if outfile is not None: global_outfiles = [outfile]
    log_result("")
    log_result("############################")
    log_result("Results {}: read_size {}, time {}".format(type_identifier, args.read_size, datetime.datetime.now()))
    log_result("  Total reads:              {:12d}".format(total_reads))
    log_result("  Total unaligned reads:    {:12d}  {:.5f} of total reads".format(unaligned_reads, 1.0 * unaligned_reads / total_reads))
    log_result("")
    log_result("  Total reference bases:    {:12d}".format(total_reference_bases))
    log_result("  Total aligned bases:      {:12d}  {:.5f} of ref bases".format(total_summaries, 1.0 * total_summaries / total_reference_bases))
    log_result("  Total softclip:           {:12d}  {:.5f} of aln bases".format(total_softclip, 1.0 * total_softclip / total_summaries))
    log_result("  Total non-softclip aln:   {:12d}  {:.5f} of ref bases".format(total_aln_nonsc, 1.0 * total_aln_nonsc / total_reference_bases))
    log_result("")
    log_result("  Total exact match:        {:12d}  {:.5f} of non-softclip aln bases".format(total_match, 1.0 * total_match / total_aln_nonsc))
    log_result("                                          {:.5f} of reference bases".format(1.0 * total_match / total_reference_bases))
    log_result("  Total mismatch:           {:12d}  {:.5f} of non-softclip aln bases".format(total_mismatch, 1.0 * total_mismatch / total_aln_nonsc))
    log_result("                                          {:.5f} of reference bases".format(1.0 * total_mismatch / total_reference_bases))
    log_result("  Total insert:             {:12d}  {:.5f} of non-softclip aln bases".format(total_insert, 1.0 * total_insert / total_aln_nonsc))
    log_result("                                          {:.5f} of reference bases".format(1.0 * total_insert / total_reference_bases))
    log_result("  Total delete:             {:12d}  {:.5f} of non-softclip aln bases".format(total_delete, 1.0 * total_delete / total_aln_nonsc))
    log_result("                                          {:.5f} of reference bases".format(1.0 * total_delete / total_reference_bases))
    if untracked_alignment_operations > 0:
        log_result("  Total unconsidered cigop: {:12d}".format(untracked_alignment_operations))
    log_result("")
    global_outfiles = None

    return 1.0 * total_match / total_aln_nonsc, 1.0 * total_match / total_reference_bases



def main(args = None):
    # get our arguments
    args = parse_args() if args is None else parse_args(args)

    # get input file
    true_reference = args.true_reference
    assembly_fasta = args.assembly_fasta
    assembly_bam = args.assembly_bam
    polisher_params = args.polisher_params
    input_files = [true_reference, assembly_fasta, assembly_bam, polisher_params]
    input_files_not_specified = False in list(map(lambda x: x is not None, input_files))
    input_files_specified_but_missing = list(filter(lambda x: x is not None and not os.path.isfile(x), input_files))
    if input_files_not_specified:
        log("Not all input files are specified")
    elif len(input_files_specified_but_missing) > 0:
        log("ERROR: Input files specified, but could not be found:\n\t{}".format(
            "\n\t".join(input_files_specified_but_missing)))
        sys.exit(1)

    ### finding files ###

    # find polisher output files
    polished_reference = "{}.fa".format(args.output_name)
    polished_reference_exists = False
    if not os.path.isfile(polished_reference):
        log("Found no polished reference")
    else:
        log("Found polished reference: {}".format(polished_reference))
        polished_reference_exists = True

    # find polished assembly reads
    polished_reads = "{}.polished_assembly.fq".format(args.output_name)
    polished_reads_exists = False
    if not os.path.isfile(polished_reads):
        log("Found no polished reads")
    else:
        log("Found polished reads: {}".format(polished_reads))
        polished_reads_exists = True
        if (args.force_realign):
            log("  Will overwrite (--force_realign set)")
            polished_reads_exists = False

    # find polished assembly bam
    polished_bam = "{}.polished_assembly.bam".format(args.output_name)
    polished_bam_exists = False
    if not os.path.isfile(polished_bam):
        log("Found no polished bam")
    else:
        log("Found polished bam: {}".format(polished_bam))
        polished_bam_exists = True
        if (args.force_realign):
            log("  Will overwrite (--force_realign set)")
            polished_bam_exists = False

    # find original assembly reads
    orig_assmebly_reads = "{}.original_assembly.fq".format(args.output_name)
    orig_assmebly_reads_exists = False
    if not os.path.isfile(orig_assmebly_reads):
        log("Found no original assembly reads")
    else:
        log("Found original assembly reads: {}".format(orig_assmebly_reads))
        orig_assmebly_reads_exists = True
        if (args.force_realign):
            log("  Will overwrite (--force_realign set)")
            orig_assmebly_reads_exists = False

    # find original assembly bam
    orig_assmebly_bam = "{}.original_assembly.bam".format(args.output_name)
    orig_assmebly_bam_exists = False
    if not os.path.isfile(orig_assmebly_bam):
        log("Found no original assembly bam")
    else:
        log("Found original assembly bam: {}".format(orig_assmebly_bam))
        orig_assmebly_bam_exists = True
        if (args.force_realign):
            log("  Will overwrite (--force_realign set)")
            orig_assmebly_bam_exists = False


    ### running process ###

    # the main process
    log("\nStarting execution")

    # polisher
    log("")
    if polished_bam_exists or polished_reads_exists or polished_reference_exists:
        log("Skipping polishing")
    elif input_files_not_specified:
        log("ERROR: Found no polishing output and missing polishing files!")
        sys.exit(1)
    else:
        log("Running polisher")
        run_polisher(assembly_fasta, assembly_bam, polisher_params, args.output_name, args.marginPolish_invocation)
    if not os.path.exists(polished_reference):
        log("ERROR: Could not find polished reference: {}".format(polished_reference))
        sys.exit(1)

    # post-polishing: generate reads
    log("")
    if polished_bam_exists or polished_reads_exists:
        log("Skipping post-polishing read generation")
    else:
        log("Running post-polshing read generation")
        run_polished_output_to_reads(polished_reference, polished_reads, args.read_size)
    if not os.path.exists(polished_reads):
        log("ERROR: Could not find polished reads: {}".format(polished_reads))
        sys.exit(1)

    # post-polishing: generate alignment
    log("")
    if polished_bam_exists:
        log("Skipping post-polishing alignment")
    else:
        log("Running post-polishing alignment")
        run_polished_read_alignment(true_reference, polished_reads, polished_bam, args.minimap2_invocation)
    if not os.path.exists(polished_bam):
        log("ERROR: Could not find polished read alignment: {}".format(polished_bam))
        sys.exit(1)

    # alignment analysis
    log("")
    log("Analyzing polishing alignment")
    polished_assembly_final_results = "{}.polished_assembly_results.txt".format(args.output_name)
    polish_exact_match_nonscaln, polish_exact_match_ref = run_alignment_analysis(
        true_reference, polished_bam, args, outfile=polished_assembly_final_results, type_identifier="polished assembly")

    ### original assembly alignments ###

    # generate original assembly reads
    log("")
    if orig_assmebly_bam_exists or orig_assmebly_reads_exists:
        log("Skipping original assembly read generation")
    else:
        log("Running original assembly read generation")
        run_polished_output_to_reads(assembly_fasta, orig_assmebly_reads, args.read_size)
    if not os.path.exists(orig_assmebly_reads):
        log("ERROR: Could not find original assembly reads: {}".format(orig_assmebly_reads))
        sys.exit(1)

    # generate original assembly alignment
    log("")
    if orig_assmebly_bam_exists:
        log("Skipping original assembly alignment")
    else:
        log("Running original assembly alignment")
        run_polished_read_alignment(true_reference, orig_assmebly_reads, orig_assmebly_bam, args.minimap2_invocation)
    if not os.path.exists(orig_assmebly_bam):
        log("ERROR: Could not find original assembly alignment: {}".format(orig_assmebly_bam))
        sys.exit(1)

    # alignment analysis
    log("")
    log("Analyzing original assembly alignment")
    orig_assembly_final_results = "{}.original_assembly_results.txt".format(args.output_name)
    orig_assembly_exact_match_nonscaln, orig_assembly_exact_match_ref = run_alignment_analysis(
        true_reference, orig_assmebly_bam, args, outfile=orig_assembly_final_results, type_identifier="original assembly")

    global global_outfiles
    global_outfiles = [orig_assembly_final_results, polished_assembly_final_results]
    log_result("Polishing difference (exact matches):")
    log_result("  Non-softclip aligned bases:")
    log_result("    Polished:    {}".format(polish_exact_match_nonscaln))
    log_result("    Original:    {}".format(orig_assembly_exact_match_nonscaln))
    log_result("    Difference:  {}".format(polish_exact_match_nonscaln - orig_assembly_exact_match_nonscaln))
    log_result("  True reference bases:")
    log_result("    Polished:    {}".format(polish_exact_match_ref))
    log_result("    Original:    {}".format(orig_assembly_exact_match_ref))
    log_result("    Difference:  {}".format(polish_exact_match_ref - orig_assembly_exact_match_ref))
    log_result("")

    ### fin ###
    log("Fin.")



if __name__ == "__main__":
    main()
