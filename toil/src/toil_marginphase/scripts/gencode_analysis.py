GENCODE = 'hg38_gencode27.bed'
LONGREAD_ONLY = 'gencode.longread_only_1k.txt'
LONGREAD_ONLY_FULL = 'hg38_gencode27.longread_only_1k.bed'


def generate_longread_only():
    longread_only = list()
    with open(LONGREAD_ONLY) as input:
        for line in input:
            longread_only.append(line.strip())

    gencode_map = dict()
    with open(GENCODE) as input:
        for line in input:
            line = line.split()
            id = line[3]
            gencode_map[id] = line

    chroms = dict()
    with open(LONGREAD_ONLY_FULL, 'w') as output:
        for lro in longread_only:
            gene = gencode_map[lro]
            output.write("\t".join(gene) + "\n")
            chrom = gene[0]
            if chrom in chroms: chroms[chrom] += 1
            else: chroms[chrom] = 1

    keys = chroms.keys()
    keys.sort()
    for key in keys:
        print("{}: {}".format(key, chroms[key]))

def get_coverage_stats():
    import glob
    print
    # for file in glob.glob("*.coverage"):
    for file in glob.glob("../workdir_cent/*.coverage"):
        coverage = 0.0
        overlap = 0.0
        with open(file) as input:
            for line in input:
                if "genome" not in line: continue
                line = line.split()
                if int(line[1]) != 0: coverage += float(line[4])
                if int(line[1]) > 1: overlap += float(line[4])
        print("{}\n\tcoverage: {}".format(file, coverage))


RMSK_LONGREAD = "INC_repeat_masker.INC_gencode_longread.bed"
def get_rmsk_percentages(rmsk):
    classes = dict()
    total = 0
    with open(rmsk) as input:
        for line in input:
            line = line.split()
            clazz = line[6]
            cover = int(line[2]) - int(line[1])
            total += cover
            if clazz in classes: classes[clazz] += cover
            else: classes[clazz] = cover

    print(rmsk)
    for k,v in classes.iteritems():
        print("%16s:\t%8d (%2d%%)" % (k, v, int(100.0*v/total)))

get_rmsk_percentages("INC_rmsk.INC_gencode.INC_repeatmask.EXC_gvcf_clbl.bed")
get_rmsk_percentages("../rmsk.bed   ")