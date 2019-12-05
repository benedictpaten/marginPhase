#!/usr/bin/env python3
import sys
import json
import math

"""
Updates a margin parameters file to calculate the repeat length emission probabilities
for hmms.
"""

if len(sys.argv) != 3:
    print("calculate_repeat_emissions.py params.json output_params.json")
    sys.exit(0)

with open(sys.argv[1]) as fh:
    p = json.loads(fh.read())
    
    alphabetSize = 4
    
    ## Calculate repeat counts  probs
     
    repeatCountsAT = p["polish"]["repeatCountSubstitutionMatrix"]["baseLogRepeatCounts_AT"]
    repeatCountsGC = p["polish"]["repeatCountSubstitutionMatrix"]["baseLogRepeatCounts_GC"]
    m = len(repeatCountsAT)
    assert len(repeatCountsGC) == m
    repeatCounts = [ math.exp(repeatCountsAT[i]) + math.exp(repeatCountsGC[i]) for i in range(m) ]
    totalProb = sum(repeatCounts)
    normalizedRepeatCounts = [ i/totalProb for i in repeatCounts ]
    print("Summed normalized repeat counts: {}".format(sum(normalizedRepeatCounts)))
    
    # Update repeat count emission probs
    
    emissionsSize = len(p["polish"]["hmm"]["emissions"])
    assert emissionsSize == alphabetSize * alphabetSize + m*m + 2*(m + alphabetSize)
    assert emissionsSize == len(p["polish"]["hmmConditional"]["emissions"])
    
    f = lambda x : [ 0.0001 if i == 0.0 else i for i in [ float("{:.4f}".format(i)) for i in x ] ]
    
    i = alphabetSize * alphabetSize + m*m + alphabetSize
    p["polish"]["hmm"]["emissions"][i: i + m] = f(normalizedRepeatCounts)
    i = alphabetSize * alphabetSize + m*m + alphabetSize * 2 + m 
    p["polish"]["hmm"]["emissions"][i: i+ m] = f(normalizedRepeatCounts)
    p["polish"]["hmmConditional"]["emissions"][i: i+ m] = f(normalizedRepeatCounts)
    
    ## Now deal with repeat count substitution matrix
    repeatCountMatrices = p["polish"]["repeatCountSubstitutionMatrix"]
    
    # Initialize combined repeat count matrix
    matrixSize = len(list(repeatCountMatrices.values())[0])
    repeatCountMatrix = [ 0.0 ] * matrixSize
    assert m*m == matrixSize
    print("Matrix size is {0}x{0}".format(m))
    
    # Sum the probabilities in non-log space
    for k in repeatCountMatrices:
        for i, prob in enumerate(repeatCountMatrices[k]):
            repeatCountMatrix[i] += math.exp(prob)
            
    # Normalize the probabilities per row 
    rowNormalizedRepeatCountMatrix = [ 0.0 ] * matrixSize
    for i in range(m):
        totalProb = sum([ repeatCountMatrix[i*m + j] for j in range(m) ])
        for j in range(m):
            rowNormalizedRepeatCountMatrix[i*m + j] =  repeatCountMatrix[i*m + j]/totalProb
            
    # Normalize the whole matrix to one
    normalizedRepeatCountMatrix = [ 0.0 ] * matrixSize
    for i in range(m):
        totalProb = sum([ repeatCountMatrix[i*m + j] for j in range(m) ])
        for j in range(m):
            normalizedRepeatCountMatrix[i*m + j] =  normalizedRepeatCounts[i] * repeatCountMatrix[i*m + j]/totalProb
    
    print("Summed normalized matrix: {}".format(sum(f(normalizedRepeatCountMatrix))))
    
    #totalProb = sum(repeatCountMatrix)
    #normalizedRepeatCountMatrix = [ i/totalProb for i in repeatCountMatrix ]
        
    i = alphabetSize*alphabetSize
    p["polish"]["hmm"]["emissions"][i:i+matrixSize] = f(normalizedRepeatCountMatrix)
    p["polish"]["hmmConditional"]["emissions"][i:i+matrixSize] = f(rowNormalizedRepeatCountMatrix)
     
    # Output the updated params file
    with open(sys.argv[2], "w") as fh2:
        fh2.write(json.dumps(p, sort_keys=True, indent=4))
        
    #with open("out.txt", "w") as fh2:
    #    fh2.write(json.dumps({ "emissions":p["polish"]["hmmConditional"]["emissions"] }))
        
    