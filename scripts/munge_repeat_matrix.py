#!/usr/bin/env python3
#from __future__ import print_function

import sys
import os

def main():
    assembly = sys.argv[1]
    
    with open(sys.argv[1]) as fh:
        seqIndex = -1
        header = None
        seq = ""
        for line in fh:
            if line[0] == '>':
                def writeSeq():
                    if header != None:
                        # Write out fasta seq
                        with open(sys.argv[1] + ("_%s" % seqIndex), 'w') as fh2:
                            fh2.write(header + "\n" + seq + "\n")
                writeSeq()
                seqIndex += 1
                header = line
                seq = ""
            else:
                seq += line[:-1]
        writeSeq()
        
if __name__ == '__main__':
    main()