
import tarfile
from contextlib import closing
import pysam

from marginphase_core import *


def merge_chunks__determine_chunk_splitting(job, chunk_identifier, left_chunk_location, right_chunk_location, chunk_boundary):

    # get read and phase block info
    l_reads, l_phase_blocks = merge_chunks__read_chunk(job, chunk_identifier, left_chunk_location)
    r_reads, r_phase_blocks = merge_chunks__read_chunk(job, chunk_identifier, right_chunk_location)

    # organize chunk comparison
    all_phase_blocks, perfect_matches, inverted_matches, shared_phase_blocks = merge_chunks__organize_reads_and_blocks(
        job, chunk_identifier, l_reads, l_phase_blocks, r_reads, r_phase_blocks)

    # recommend strategy
    split_position, phase_block, invert_right, decision_summary = merge_chunks__recommend_merge_strategy(
        job, chunk_identifier, chunk_boundary, perfect_matches, inverted_matches, shared_phase_blocks)

    # implement strategy
    left_reads_writing, right_reads_writing, vcf_split_pos, vcf_right_phase_action = merge_chunks__specify_split_action(
        job, chunk_identifier, split_position, phase_block, invert_right,
        l_reads, l_phase_blocks, r_reads, r_phase_blocks)

    # log summarization
    job.fileStore.logToMaster(("{}:merge_chunks:determine_chunk_splitting: "
                              "read merge action: write {} from left, {} from right")
                              .format(chunk_identifier, len(left_reads_writing), len(right_reads_writing)))
    job.fileStore.logToMaster(("{}:merge_chunks:determine_chunk_splitting: call merge action: "
                              "split at {}, right action {}")
                              .format(chunk_identifier, vcf_split_pos, vcf_right_phase_action))

    # return
    return left_reads_writing, right_reads_writing, vcf_split_pos, vcf_right_phase_action, decision_summary


def merge_chunks__organize_reads_and_blocks(job, chunk_identifier, l_reads, l_phase_blocks, r_reads, r_phase_blocks):

    # get all phase blocks with same start pos
    all_phase_blocks = list(set(l_phase_blocks.keys()).union(set(r_phase_blocks.keys())))
    all_phase_block_count = len(all_phase_blocks)

    # get all phase blocks with same start pos
    shared_phase_blocks = list(set(l_phase_blocks.keys()).intersection(set(r_phase_blocks.keys())))
    shared_phase_blocks.sort()

    # get phase blocks by start and end pos
    l_pb_uniq_ids = list(map(lambda x: "{}-{}".format(x[PB_START_POS], x[PB_END_POS]), l_phase_blocks.values()))
    r_pb_uniq_ids = set(map(lambda x: "{}-{}".format(x[PB_START_POS], x[PB_END_POS]), r_phase_blocks.values()))

    # find all matches
    perfect_matches = list()
    inverted_matches = list()
    for l_pb_uniq in l_pb_uniq_ids:
        if l_pb_uniq in r_pb_uniq_ids:
            # we know phase block positions align
            shared_phase_block = int(l_pb_uniq.split("-")[0])
            phase_block_median_pos = int((int(l_pb_uniq.split("-")[0]) + int(l_pb_uniq.split("-")[1]))/2.0)

            # get all reads in phase block (get from both blocks, in case some were skipped)
            reads_in_phase_block = set()
            for reads in [l_phase_blocks[shared_phase_block][PB_HAP1_READS],l_phase_blocks[shared_phase_block][PB_HAP2_READS],
                          r_phase_blocks[shared_phase_block][PB_HAP1_READS],r_phase_blocks[shared_phase_block][PB_HAP2_READS]]:
                for read in reads:
                    reads_in_phase_block.add(read)

            # get reads spannign median position
            reads_spanning_median = list(filter(
                lambda x:   l_reads[x][RD_ALN_START] <= phase_block_median_pos <= l_reads[x][RD_ALN_END]
                            if x in l_reads else
                            r_reads[x][RD_ALN_START] <= phase_block_median_pos <= r_reads[x][RD_ALN_END],
                reads_in_phase_block))

            # perfect match?
            perfect_matched_reads = list(filter(
                lambda x:   l_reads[x][RD_HAPLOTYPE_TAG] == r_reads[x][RD_HAPLOTYPE_TAG]
                            if x in l_reads and x in r_reads else False,
                reads_spanning_median
            ))
            if len(perfect_matched_reads) == len(reads_spanning_median):
                perfect_matches.append(l_pb_uniq)
                continue

            # inverted match?
            inverted_matched_reads = list(filter(
                lambda x: l_reads[x][RD_HAPLOTYPE_TAG].replace("h1","<TMP>").replace("h2","h1").replace("<TMP>","h2")
                          == r_reads[x][RD_HAPLOTYPE_TAG] if x in l_reads and x in r_reads else False,
                reads_spanning_median
            ))
            if len(inverted_matched_reads) == len(reads_spanning_median):
                inverted_matches.append(l_pb_uniq)
                continue

    # loggit
    job.fileStore.logToMaster("{}:merge_chunks:organize_reads_and_blocks: Found {} distinct phase blocks"
        .format(chunk_identifier, all_phase_block_count))
    job.fileStore.logToMaster("{}:merge_chunks:organize_reads_and_blocks: Found {} ({}%) perfect matches"
        .format(chunk_identifier, len(perfect_matches), percent(len(perfect_matches), all_phase_block_count)))
    job.fileStore.logToMaster("{}:merge_chunks:organize_reads_and_blocks: Found {} ({}%) inverted matches"
        .format(chunk_identifier, len(inverted_matches), percent(len(inverted_matches), all_phase_block_count)))
    job.fileStore.logToMaster("{}:merge_chunks:organize_reads_and_blocks: Found {} ({}%) matched phase starts"
        .format(chunk_identifier, len(shared_phase_blocks), percent(len(shared_phase_blocks), all_phase_block_count)))

    # return what we found
    return all_phase_blocks, perfect_matches, inverted_matches, shared_phase_blocks


def merge_chunks__recommend_merge_strategy(job, chunk_identifier, chunk_boundary, perfect_matches,
                                           inverted_matches, shared_phase_blocks):

    # helper function
    def pick_closest_elem_to_chunk_boundary(elems):
        elems.sort(key=lambda x: abs(chunk_boundary - sum(map(int, str(x).split("-"))) / len(str(x).split("-"))))
        return elems[0]

    # description of merge strategy
    split_position, phase_block, invert_right, decision_summary = None, None, None, None

    # case1: perfect match:
    #   READS: left of and spanning split_pos from left chunk, starting after split_pos from right chunk
    #   CALLS: left of split pos from left VCF, right of split pos from right VCF
    if len(perfect_matches) > 0:
        parts = map(int, pick_closest_elem_to_chunk_boundary(perfect_matches).split("-"))
        split_position = int(sum(parts) / len(parts))
        phase_block = parts[0]
        invert_right = False
        decision_summary = "PERFECT_MATCH"
        job.fileStore.logToMaster("{}:merge_chunks:recommend_merge_strategy: Found perfect match at "
                                  "pos {} in phase block {}".format(chunk_identifier, split_position, phase_block))


    # case2: perfect match but inverted haploptyes
    #   READS: left of and spanning split_pos from left chunk, starting after split_pos from right chunk
    #       reverse haplotype of included reads in phase_block from right chunk
    #   CALLS: left of split pos from left VCF, right of split pos from right VCF
    #       reverse phasing of calls in phase block from right chunk
    elif len(inverted_matches) > 0:
        parts = map(int, pick_closest_elem_to_chunk_boundary(inverted_matches).split("-"))
        split_position = int(sum(parts) / len(parts))
        phase_block = parts[0]
        invert_right = True
        decision_summary = "INVERT_MATCH"
        job.fileStore.logToMaster(
            "{}:merge_chunks:recommend_merge_strategy: Found inverted match at "
            "pos {} in phase block {}".format(chunk_identifier, split_position, phase_block))

    # case3: found a phase block starting at the same posistion in each chunk
    #   READS: finishing before split_pos from left chunk, starting after split pos from right chunk
    #       reads spanning split_pos get hap info from left before split_pos, and hap info from right after and including
    #   CALLS: left of split pos from left VCF, right of split pos from right VCF
    elif len(shared_phase_blocks) > 0:
        phase_block = pick_closest_elem_to_chunk_boundary(shared_phase_blocks)
        split_position = phase_block
        invert_right = False
        decision_summary = "PHASE_START_MATCH"
        job.fileStore.logToMaster("{}:merge_chunks:recommend_merge_strategy: Found phase block start match "
                                  "at {}".format(chunk_identifier, phase_block))

    # case4: no matching phase blocks
    #   READS: finishing before split_pos from left chunk, reads spanning split_pos get phasing finishing left of
    #       split_pos from left chunk, phasing in phase blocks spanning split_pos get split in two, phasing in phase
    #       blocks starting after split_pos from right chunk
    #   CALLS: starting left of split_pos from left VCF, starting after split_pos from right VCF, calls from right
    #       phase block spanning split_pos get new phase_block
    else:
        phase_block = None
        split_position = chunk_boundary
        invert_right = False
        decision_summary = "NO_MATCH"
        job.fileStore.logToMaster("{}:merge_chunks:recommend_merge_strategy: Found no match, creating "
                                  "new phase block at {}".format(chunk_identifier, split_position))

    # return data
    return split_position, phase_block, invert_right, decision_summary


def merge_chunks__specify_split_action(job, chunk_identifier, split_position, phase_block, invert_right,
                                       l_reads, l_phase_blocks, r_reads, r_phase_blocks):

    # describes read inclusion and modifications to haplotype string
    left_reads_writing = dict()
    right_reads_writing = dict()
    right_phase_block_conversion = dict()

    for read in l_reads.values():
        # this read belongs wholly to the right chunk
        if read[RD_ALN_START] > split_position:
            continue

        # this read belongs wholly to the left chunk
        elif read[RD_ALN_END] <= split_position:
            left_reads_writing[read[RD_ID]] = None

        # this read spans the split_position (and needs analysis)
        elif read[RD_ALN_START] <= split_position and read[RD_ALN_END] > split_position:

            # case4: new phase block created at split pos
            if phase_block is None:
                l_read = read
                r_read = r_reads[read[RD_ID]]
                new_hap_str, old_right_phase_block = merge_chunks__create_new_phase_block_at_position(
                    job, chunk_identifier, split_position, l_read, r_read)
                left_reads_writing[read[RD_ID]] = [[read[RD_HAPLOTYPE_TAG], new_hap_str]]
                if old_right_phase_block is not None: # None indicates read was not haplotyped (ie, h0)
                    right_phase_block_conversion[old_right_phase_block] = split_position

                # santity check (should never find two phase blocks from right)
                if len(right_phase_block_conversion) > 1:
                    raise UserError("SANITY_CHECK_FAILURE:{}: got inconsistent phase blocks ({}) spanning {} for read {}"
                                    .format(chunk_identifier, right_phase_block_conversion.keys(), split_position, read[RD_ID]))

            # case3: take hap info before split_pos from left, after right.  phase block exists at split_pos
            elif phase_block == split_position:
                l_read = l_reads[read[RD_ID]]
                haps = list(filter(lambda x: x[RPB_BLOCK_ID] < split_position, l_read[RD_PHASE_BLOCKS]))
                if read[RD_ID] in r_reads:
                    r_read = r_reads[read[RD_ID]]
                    haps.extend(list(filter(lambda x: x[RPB_BLOCK_ID] >= split_position, r_read[RD_PHASE_BLOCKS])))
                haps.sort(key=lambda x: x[RPB_BLOCK_ID])
                new_hap_str = ";".join(map(merge_chunks__encode_phase_info, haps))
                left_reads_writing[read[RD_ID]] = [[read[RD_HAPLOTYPE_TAG], new_hap_str]]

            # case2, case1:
            else:
                left_reads_writing[read[RD_ID]] = None

    # get right reads we care about (reads in the spanning-the-split_pos phase chunk)
    # (everthing before split_pos comes from left chunk, everything after is unchanged)
    analysis_read_ids = set()
    if phase_block is None:
        if len(right_phase_block_conversion) == 0:
            job.fileStore.logToMaster("{}:merge_chunks:specify_split_action: No reads spanning {} were found!"
                .format(chunk_identifier, split_position))
        else:
            analysis_phase_block_id = list(right_phase_block_conversion.keys())[0]
            analysis_phase_block = r_phase_blocks[analysis_phase_block_id]
            analysis_read_ids = analysis_phase_block[PB_HAP1_READS].union(analysis_phase_block[PB_HAP2_READS])
    else:
        analysis_phase_block = r_phase_blocks[phase_block]
        analysis_read_ids = analysis_phase_block[PB_HAP1_READS].union(analysis_phase_block[PB_HAP2_READS])

    for read_id in analysis_read_ids:
        read = r_reads[read_id]

        # this read belongs wholly to the left chunk
        if read[RD_ALN_END] <= split_position:
            continue

        # this was analyzed with the left reads
        elif read_id in left_reads_writing:
            continue

        # now we need to analyize - we know these reads start after split_pos
        # case4
        if phase_block is None:
            if len(right_phase_block_conversion) == 0:
                raise UserError("SANITY_CHECK_FAILURE:{}: new phase block determined, but no conversion for read {}"
                                .format(chunk_identifier, read_id))
            pb_from = list(right_phase_block_conversion.keys())[0]
            pb_to = right_phase_block_conversion[pb_from]
            new_hap_str = read[RD_HAPLOTYPE_TAG].replace("p{},".format(pb_from), "p{},".format(pb_to))
            right_reads_writing[read_id] = [[read[RD_HAPLOTYPE_TAG], new_hap_str]]

        # case2
        elif invert_right:
            h1_str = "h1,p{}".format(phase_block)
            h1_tmp = "h1,p<TMP>"
            h2_str = "h2,p{}".format(phase_block)
            new_hap_str = read[RD_HAPLOTYPE_TAG].replace(h1_str, h1_tmp).replace(h2_str, h1_str).replace(h1_tmp, h2_str)
            right_reads_writing[read_id] = [[read[RD_HAPLOTYPE_TAG], new_hap_str]]

        # case1, case3
        else:
            pass

    # summarize vcf
    vcf_split_position = split_position
    vcf_right_phase_conversion = dict()

    if phase_block is None: # case 4
        if len(right_phase_block_conversion) != 0:
            # this also requires a new phase block when started
            vcf_right_phase_conversion = right_phase_block_conversion
        else: # no reads span this, so no action
            pass
    elif invert_right: # case 2
        vcf_right_phase_conversion = {phase_block: ACTION_INVERT}
    else:
        # case1 or case 3: no action
        pass

    # finish
    return left_reads_writing, right_reads_writing, vcf_split_position, vcf_right_phase_conversion


def merge_chunks__parse_phase_info(phase_info):
    parts = phase_info.split(",")
    if len(parts) != 4 or not (parts[0].startswith('h') and parts[1].startswith('p') and parts[2].startswith('r')
                               and parts[3].startswith('l')):
        return None, None, None, None
    parts = map(lambda x: int(x[1:]), parts)
    haplotype, phase_block, read_start, read_length = parts[0], parts[1], parts[2], parts[3]
    return haplotype, phase_block, read_start, read_length


def merge_chunks__encode_phase_info(read_data):
    haplotype = 0 if read_data[RPB_IS_HAP1] is None else (1 if read_data[RPB_IS_HAP1] else 2)
    phase_block = read_data[RPB_BLOCK_ID]
    read_start = read_data[RPB_READ_START]
    read_length = read_data[RPB_READ_LENGTH]
    return "h{},p{},r{},l{}".format(haplotype, phase_block, read_start, read_length)


def merge_chunks__save_read_info(all_reads, read):
    # get info
    read_id = read.query_name
    align_start = read.reference_start
    align_end = read.reference_end

    # save storage data
    read_data = dict()
    all_reads[read_id] = read_data
    read_data[RD_ID] = read_id
    read_data[RD_ALN_START] = align_start
    read_data[RD_ALN_END] = align_end
    read_data[RD_HAPLOTYPE_TAG] = read.get_tag(TAG_HAPLOTYPE)
    read_data[RD_PHASE_BLOCKS] = list()

    # return data
    return read_data


def merge_chunks__save_phase_block_info(job, chunk_identifier, phase_blocks, phase_info, read_id):
    # get info
    haplotype, phase_block, read_start, read_length = merge_chunks__parse_phase_info(phase_info)
    if haplotype is None:
        raise UserError("SANITY_CHECK_FAILURE:{}: malformed phase info for read {}: {}"
                        .format(chunk_identifier, read_id, phase_info))
    if haplotype == 0:
        return None

    # get storage data
    if phase_block in phase_blocks:
        phase_block_data = phase_blocks[phase_block]
    else:
        phase_block_data = dict()
        phase_blocks[phase_block] = phase_block_data
        phase_block_data[PB_BLOCK_ID] = phase_block
        phase_block_data[PB_HAP1_READS] = set()
        phase_block_data[PB_HAP2_READS] = set()
        phase_block_data[PB_END_POS] = None
        phase_blocks[phase_block] = phase_block_data

    # read pb info
    read_phase_block_info = dict()
    read_phase_block_info[RPB_BLOCK_ID] = phase_block
    read_phase_block_info[RPB_READ_START] = read_start
    read_phase_block_info[RPB_READ_LENGTH] = read_length
    read_phase_block_info[RPB_IS_HAP1] = None

    # save read
    if haplotype == 1:
        phase_block_data[PB_HAP1_READS].add(read_id)
        read_phase_block_info[RPB_IS_HAP1] = True
    elif haplotype == 2:
        phase_block_data[PB_HAP2_READS].add(read_id)
        read_phase_block_info[RPB_IS_HAP1] = False
    else:
        raise UserError("SANITY_CHECK_FAILURE:{}: unexpected haplotype in phase_info for read {}: {}"
                        .format(chunk_identifier, read_id, phase_info))

    # return read phase data
    return read_phase_block_info


def merge_chunks__read_chunk(job, chunk_identifier, chunk_location):
    # log
    job.fileStore.logToMaster("{}:merge_chunks:read_chunk: reading chunk {}"
                              .format(chunk_identifier, os.path.basename(chunk_location)))

    # data storage
    phase_blocks = dict()
    reads = dict()
    read_count = 0
    failed_reads = 0

    with closing(pysam.AlignmentFile(chunk_location, 'rb' if chunk_location.endswith("bam") else 'r')) as aln:

        start = time.time()
        for read in aln.fetch():
            # get read data
            read_id = read.query_name
            read_count += 1

            # find haplotype tag
            for tag in [TAG_HAPLOTYPE, TAG_CHUNK_ID]:
                if not read.has_tag(tag):
                    job.fileStore.logToMaster("{}:merge_chunks:read_chunk: read {} had no {} tag"
                                              .format(chunk_identifier, read_id, tag))
                    failed_reads += 1
                    continue

            # save read data
            read_data = merge_chunks__save_read_info(reads, read)

            # save haplotpye data
            haplotype_tags = read.get_tag(TAG_HAPLOTYPE).split(";")
            for pb_tag in haplotype_tags:
                rpb_info = merge_chunks__save_phase_block_info(job, chunk_identifier, phase_blocks, pb_tag, read_id)
                if rpb_info is not None: read_data[RD_PHASE_BLOCKS].append(rpb_info)

        job.fileStore.logToMaster("{}:merge_chunks:read_chunk: read {} reads ({}s)"
                                  .format(chunk_identifier, read_count, int(time.time() - start)))

    # finish phase block analysis
    phase_block_ids = list(phase_blocks.keys())
    phase_block_ids.sort()
    prev_pb = None
    for pb_id in phase_block_ids:
        curr_pb = phase_blocks[pb_id]
        if prev_pb is not None: prev_pb[PB_END_POS] = curr_pb[PB_START_POS]
        prev_pb = curr_pb
    # we aren't going to use this last one anyway
    prev_pb[PB_END_POS] = prev_pb[PB_START_POS]

    # return chunk data
    return reads, phase_blocks


def merge_chunks__create_new_phase_block_at_position(job, chunk_identifier, split_position, l_read, r_read):
    # get all documented haplotpyes
    l_haps = l_read[RD_PHASE_BLOCKS]
    r_haps = r_read[RD_PHASE_BLOCKS]

    # data we want at the end
    haps = list()
    old_right_phase_block = None

    # get desired (and modified) haplotypes from l_read
    for hap in l_haps:
        if hap[RPB_BLOCK_ID] >= split_position:
            continue  # belongs to r_read
        elif hap[RPB_BLOCK_ID] + hap[RPB_READ_LENGTH] < split_position:
            haps.append(hap)  # before split (or h0 read)
        else:  # this read needs to be split
            new_hap = {
                RPB_IS_HAP1: hap[RPB_IS_HAP1],
                RPB_BLOCK_ID: hap[RPB_BLOCK_ID],
                RPB_READ_START: hap[RPB_READ_START],
                RPB_READ_LENGTH: split_position - hap[RPB_BLOCK_ID]
            }
            haps.append(new_hap)

    # get desired (and modified) haplotypes from r_read
    for hap in r_haps:
        if hap[RPB_BLOCK_ID] >= split_position:
            haps.append(hap)  # after split
        elif hap[RPB_BLOCK_ID] + hap[RPB_READ_LENGTH] < split_position:
            continue  # belongs to l_read
        else:  # this read needs to be split
            split_diff = split_position - hap[RPB_BLOCK_ID]
            new_hap = {
                RPB_IS_HAP1: hap[RPB_IS_HAP1],
                RPB_BLOCK_ID: split_position,
                RPB_READ_START: hap[RPB_READ_START] + split_diff,
                RPB_READ_LENGTH: hap[RPB_READ_LENGTH] - split_diff
            }
            haps.append(new_hap)
            # sanity check and save old haplotype
            if old_right_phase_block is not None:
                raise UserError("{}:SANITY_CHECK_FAILURE: found multiple phase_blocks ({}, {}) spanning split_position "
                                "{} for read {}:".format(chunk_identifier, old_right_phase_block, hap[RPB_BLOCK_ID],
                                                         split_position, l_read[RD_ID]))
            old_right_phase_block = hap[RPB_BLOCK_ID]

    # edge case for when l_hap is haplotyped after split pos, but and r_haps is not-haplotyped
    if len(haps) == 0:
        job.fileStore.logToMaster("{}:merge_chunks:create_new_phase_block: output no haplotypes for read {} ({}-{}) "
                                  "spanning split_position {}, with haplotypes (left) {} and (right) {}"
                                  .format(chunk_identifier, l_read[RD_ID], l_read[RD_ALN_START], l_read[RD_ALN_END],
                                          split_position, l_read[RD_HAPLOTYPE_TAG], r_read[RD_HAPLOTYPE_TAG]))
        haps.append({RPB_IS_HAP1: None, RPB_BLOCK_ID: 0, RPB_READ_START: 0, RPB_READ_LENGTH: 0})

    # save haploptyes
    haps.sort(key=lambda x: x[RPB_BLOCK_ID])
    new_hap_str = ";".join(map(merge_chunks__encode_phase_info, haps))


    return new_hap_str, old_right_phase_block


def merge_chunks__append_sam_reads(job, chunk_identifier, input_sam_file, output_sam_file, included_reads,
                                   excluded_reads):
    # data to track
    header_lines = 0
    written_read_ids = set()
    reads_not_written = 0
    explicitly_excluded = 0
    modified_writes = 0

    # include header?
    include_header = not os.path.isfile(output_sam_file)
    if include_header:
        job.fileStore.logToMaster("{}:merge_chunks:append_sam_reads: output sam file does not exist, writing header"
                                  .format(chunk_identifier))

    # read and write
    with open(output_sam_file, 'a') as output, open(input_sam_file, 'r') as input:
        for line in input:
            # header
            if line.startswith("@"):
                if include_header:
                    output.write(line)
                    header_lines += 1
                continue
            read_id = line.split("\t")[0]

            # already written in a previous chunk
            if read_id in excluded_reads:
                reads_not_written += 1
                explicitly_excluded += 1
                continue

            # should be written in this chunk
            elif read_id in included_reads:
                if included_reads[read_id] is not None:
                    for substitution in included_reads[read_id]:
                        # sanity check
                        if type(substitution) != list and len(substitution) != 2:
                            raise UserError("{}:merge_chunks:append_sam_reads: malformed substitution for {} in {}: {}"
                                            .format(chunk_identifier, read_id, included_reads[read_id], substitution))
                        # replace old hap str for new one
                        line = line.replace(substitution[0], substitution[1])
                    modified_writes += 1
                output.write(line)
                written_read_ids.add(read_id)

            # edge case for inlusion of all reads in last chunk
            elif None in included_reads:
                output.write(line)
                written_read_ids.add(read_id)

            # if not explicitly told to write, read is not written
            else:
                reads_not_written += 1

    if include_header:
        job.fileStore.logToMaster("{}:merge_chunks:append_sam_reads: wrote {} header lines from {} to {}"
                                  .format(chunk_identifier, header_lines,
                                          os.path.basename(input_sam_file), os.path.basename(output_sam_file)))
    job.fileStore.logToMaster("{}:merge_chunks:append_sam_reads: wrote {} reads ({} modified) with {} excluded "
                              "({} explicitly so) from {} to {}"
                              .format(chunk_identifier, len(written_read_ids), modified_writes,
                                      reads_not_written, explicitly_excluded,
                                      os.path.basename(input_sam_file), os.path.basename(output_sam_file)))
    return written_read_ids


def merge_chunks__append_vcf_calls(job, chunk_identifier, input_vcf_file, output_vcf_file, start_pos, end_pos,
                                   vcf_conversion, mp_identifier=None):

    # include header?
    include_header = not os.path.isfile(output_vcf_file)
    if include_header:
        job.fileStore.logToMaster("{}:merge_chunks:append_vcf_calls: output vcf file does not exist, writing header"
                                  .format(chunk_identifier))

    # data we track
    written_header_lines = 0
    written_lines = 0
    lines_outside_boundaries = 0

    # read and write
    with open(output_vcf_file, 'a') as output, open(input_vcf_file, 'r') as input:
        first_analyzed_line = True # may need to manage the phase set (only for the first phase of a chunk)
        for line in input:
            if line.startswith("#"):  # header
                if include_header:
                    output.write(line)
                    written_header_lines += 1
                continue

            # break line into parts
            line = line.rstrip().split("\t")
            position = int(line[1])
            if position < start_pos or position >= end_pos: #only include positions in given range
                lines_outside_boundaries += 1
                continue

            # get info and tags (and positions)
            info = line[-1].split(":")
            tags = line[-2].split(":")
            genotype_tag_idx = None
            phase_block_tag_idx = None
            idx = 0
            for tag in tags:
                if tag == TAG_GENOTYPE: genotype_tag_idx = idx
                if tag == TAG_PHASE_SET: phase_block_tag_idx = idx
                idx += 1
            if genotype_tag_idx is None:
                raise UserError("{}:SANITY_CHECK_FAILURE: malformed vcf {} phasing line (no {} tag): {}"
                                .format(chunk_identifier, os.path.basename(input_vcf_file), TAG_GENOTYPE, line))
            if phase_block_tag_idx is None:
                raise UserError("{}:SANITY_CHECK_FAILURE: malformed vcf {} phasing line (no {} tag): {}"
                                .format(chunk_identifier, os.path.basename(input_vcf_file), TAG_PHASE_SET, line))

            # phase block
            initial_phase_block = int(info[phase_block_tag_idx])
            reverse_phasing = False
            updated_phase_block = str(initial_phase_block)
            if initial_phase_block in vcf_conversion:
                if vcf_conversion[initial_phase_block] == ACTION_INVERT:
                    reverse_phasing = True
                else:
                    updated_phase_block = str(vcf_conversion[initial_phase_block])
            # set phase block (may be same as old phase block)
            info[phase_block_tag_idx] = updated_phase_block

            # get genotype
            genotype = info[genotype_tag_idx]
            continued_phase_block = "|" in genotype
            new_phase_block = "/" in genotype
            if not (continued_phase_block ^ new_phase_block):
                raise UserError("{}:SANITY_CHECK_FAILURE: Malformed vcf {} phasing line (unexpected genotype): {}"
                                .format(chunk_identifier, os.path.basename(input_vcf_file), line))
            genotype = genotype.split("|") if continued_phase_block else genotype.split("/")
            # make updates to genotype (potentially)
            if reverse_phasing:
                genotype.reverse()
            if first_analyzed_line and initial_phase_block != updated_phase_block:
                new_phase_block = True
            # save genotype
            genotype = ("/" if new_phase_block else "|").join(genotype)
            info[genotype_tag_idx] = genotype


            # add identifier (if appropriate)
            if mp_identifier is not None:
                info.append(mp_identifier)
                tags.append(TAG_MARGIN_PHASE_IDENTIFIER)
                line[-2] = ":".join(tags)

            # cleanup
            line[-1] = ":".join(map(str, info))
            line = "\t".join(line) + "\n"
            output.write(line)
            written_lines += 1
            first_analyzed_line = False

    # loggit
    job.fileStore.logToMaster(
        "{}:merge_chunks:append_vcf_calls: wrote {} calls ({} skipped from being outside boundaries {}-{}) from {} to {}"
            .format(chunk_identifier, written_lines, lines_outside_boundaries, start_pos, end_pos,
                    os.path.basename(input_vcf_file), os.path.basename(output_vcf_file)))


def merge_chunks__extract_chunk_tarball(job, config, tar_work_dir, chunk):
    # prep
    os.mkdir(tar_work_dir)
    tar_file = os.path.join(tar_work_dir, "chunk.tar.gz")

    # get file
    job.fileStore.readGlobalFile(chunk[CI_OUTPUT_FILE_ID], tar_file, mutable=True)
    with tarfile.open(tar_file, 'r') as tar:
        tar.extractall(tar_work_dir)

    # find desired files
    sam_unified, vcf = None, None
    for name in os.listdir(tar_work_dir):
        if name.endswith(VCF_SUFFIX): vcf = name
        elif name.endswith(SAM_UNIFIED_SUFFIX): sam_unified = name
    sam_unified_file = None if sam_unified is None else os.path.join(tar_work_dir, sam_unified)
    vcf_file = None if vcf is None else os.path.join(tar_work_dir, vcf)

    # return file locations
    return sam_unified_file, vcf_file
