#!/usr/bin/env python3
#from __future__ import print_function
import argparse
import subprocess
import sys
import os
import math

def main():
    trueAssembly = sys.argv[1]
    assembly = open(sys.argv[2]).readlines()
    verbose = True if len(sys.argv) == 4 and sys.argv[-1] == "verbose" else False
    
    read = []
    for line in assembly[1:]:
        if ">" in line:
            break
        read.append(line[:-1])
    read = "".join(read)
    
    # This spike in was used to confirm that the x sequence that is printed by the code  is the true assembly 
    # and the y sequence is the assembly
    # so x gaps are insertions with respect to the true assembly
    # and y gaps are deletions with respect to the true assembly
    #read = read[:30000] + "GATTACAGATTACA" + read[30000:]
    
    print("Read len:", len(read))
    #read = read[:70000]
    fh = open("temp.fa", 'w')
    fh.write(">hello\n%s\n" % read)
    fh.close()
    
    subprocess.call("cPecanLastz %s temp.fa --format=axt > temp.axt" % (trueAssembly,), shell=True)
    
    axt = open("temp.axt", "r").readlines()
    axt = [ i for i in axt if i[0] != '#']
    
    identities = []
    totalMatches = 0
    totalXLen = 0
    totalYLen = 0
    totalXGaps = 0
    totalYGaps = 0
    totalXGapEvents = 0
    totalYGapEvents = 0
    totalMismatches = 0
    for j in range(len(axt)):
        if "hello" in axt[j]:
    
            x = axt[j+1]
            y = axt[j+2]
            
            xLen = len([ i for i in x if i != '-'])
            yLen = len([ i for i in y if i != '-' ]) 
 
            if xLen > 50000 and yLen > 50000:
                xGaps = len(x) - xLen
                yGaps = len(y) - yLen
                xGapEvents = len([ i for i in range(1, len(x)) if x[i] == '-' and x[i-1] != '-'  ])
                yGapEvents = len([ i for i in range(1, len(y)) if y[i] == '-' and y[i-1] != '-'  ])
                mismatches = sum([ 1 if (x[i] != y[i] and x[i] != '-' and y[i] != '-') else 0 for i in range(len(x)) ])
                mismatchLocations = [ (i, x[i-10:i+10], y[i-10:i+10]) for i in range(len(x)) if (x[i] != y[i] and x[i] != '-' and y[i] != '-')]
                matches = sum([ 1 if x[i] == y[i] else 0 for i in range(len(x)) ])

                identity = matches * 2 / float(xLen + yLen)
                
                identities.append((xLen, yLen, matches, identity, mismatches, xGaps, yGaps, mismatchLocations))
                
                totalMatches += matches
                totalXLen += xLen
                totalYLen += yLen
                totalXGaps += xGaps
                totalYGaps += yGaps
                totalXGapEvents += xGapEvents
                totalYGapEvents += yGapEvents
                totalMismatches += mismatches
    
    identities.sort()
    
    for xLen, yLen, matches, identity, mismatches, xGaps, yGaps, mismatchLocations in identities:
        print("xLen: %s, yLen: %s, matches: %s, identity: %s, mismatches: %s, xGaps: %s, yGaps: %s" % (xLen, yLen, matches, identity, mismatches, xGaps, yGaps))
        if verbose:
            for location, xString, yString in mismatchLocations:
                print(" Mismatch", location, xString, yString)

    print("total-xLen: %s, total-yLen: %s, total-matches: %s, total-identity: %s, total-mismatches: %s, total-xGaps: %s, total-yGaps: %s, total-xGapEvents: %s, total-yGapEvents: %s" % 
		(totalXLen, totalYLen, totalMatches, 2.0*totalMatches/float(totalXLen + totalYLen + 0.0001), totalMismatches, totalXGaps, totalYGaps, totalXGapEvents, totalYGapEvents))

if __name__ == '__main__':
    main()
