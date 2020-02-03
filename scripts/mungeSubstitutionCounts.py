"""
Converts Ryan's substitution matrix to one used by margin
"""

counts = { k:{ i:{ j:0 for j in "ACTG" } for i in "ACTG" } for k in "FR" }
fh = open("subs.txt")
line = fh.readline()
while line != '':
    assert line[0] == '>'
    base, strand = line[1], line[3]
    aCount = float(fh.readline().split(",")[0])
    cCount = float(fh.readline().split(",")[0])
    gCount = float(fh.readline().split(",")[0])
    tCount = float(fh.readline().split(",")[0])
    total = aCount + cCount + gCount + tCount
    counts[strand][base]["A"] = aCount / total
    counts[strand][base]["C"] = cCount / total
    counts[strand][base]["G"] = gCount / total
    counts[strand][base]["T"] = tCount / total
    line = fh.readline()
    
for strand in "FR":
    print("Strand: ", strand)
    for refBase in "ACGT":
        rc = lambda b : { "A":"T", "C":"G", "T":"A", "G":"C" }[b]
        #refBase = refBase if strand == "F" else rc(refBase)
        for readBase in "ACGT":
            readBase = readBase if strand == "F" else rc(readBase)
            print("{:.3f}".format(counts[strand][refBase][readBase]) + ", ", end="")
        print()
    print("")
    
#print(counts)