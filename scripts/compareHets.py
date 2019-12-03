

import sys
 
if len(sys.argv) != 3:
	print("compareHets.py trueHetsFile predictedHetsFile")
	sys.exit(0)

trueHetsFile = sys.argv[1]
predictedHetsFile = sys.argv[2]

def rle(s):
	return [ s[i] for i in range(len(s)) if i == 0 or s[i] != s[i-1] ]

def compareRLEs(str1, str2):
	return rle(str1[5:-5]) != rle(str2[5:-5])

def getMismatches(mismatchFile):
	mismatches = {}
	mismatchLines = []
	with open(mismatchFile) as fh:
		for line in fh:
			tokens = line.split()
			if len(tokens) > 0 and tokens[0] == "Mismatch":
				mismatchLines.append((tokens[-2], tokens[-1], int(tokens[1])))
				if tokens[-1] > tokens[-2]:
					mismatches[tokens[-2]] = (tokens[-1], int(tokens[1]))
				else:
					mismatches[tokens[-1]] = (tokens[-2], int(tokens[1]))	
	return mismatches, mismatchLines

trueMismatches, trueMismatchLines = getMismatches(trueHetsFile)
predictedMismatches, predictedMismatchLines = getMismatches(predictedHetsFile)

totalCommonMismatches = set(trueMismatches.keys()).intersection(set(predictedMismatches.keys()))

print("Total common mismatches %s, total true mismatches %s, total predicted mismatches %s" % \
	(len(totalCommonMismatches), len(trueMismatches), len(predictedMismatches)))

print("Predicted mismatches")
for mismatch1, mismatch2, location in predictedMismatchLines:
	print("Predicted mismatch: %s %s %s %s %s" % (mismatch1, mismatch2, location, mismatch1 in trueMismatches or mismatch2 in trueMismatches, compareRLEs(mismatch1, mismatch2)))

print("True mismatches")
for mismatch1, mismatch2, location in trueMismatchLines:
        print("Predicted mismatch: %s %s %s %s %s" % (mismatch1, mismatch2, location, mismatch1 in predictedMismatches or mismatch2 in predictedMismatches, compareRLEs(mismatch1, mismatch2)))

