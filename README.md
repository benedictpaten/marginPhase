# MarginPhase #

MarginPhase is a program for simultaneous haplotyping and genotyping.

## Dependencies ##
cmake version 3.7 (or higher):
```
wget https://cmake.org/files/v3.7/cmake-3.7.2-Linux-x86_64.sh && mkdir /opt/cmake && sh cmake-3.7.2-Linux-x86_64.sh --prefix=/opt/cmake --skip-license && ln -s /opt/cmake/bin/cmake /usr/local/bin/cmake
```

If you're on Ubuntu:
```
apt-get -y install git make gcc g++ autoconf bzip2 lzma-dev zlib1g-dev libcurl4-openssl-dev libcrypto++-dev libpthread-stubs0-dev libbz2-dev liblzma-dev
```

- Check out the repository:
```
git clone git@github.com:benedictpaten/marginPhase.git
```

- Check out submodules:
```
cd marginPhase
git submodule update --init
```

- Make build directory:
```
mkdir build
cd build
```

- Generate Makefile with cmake:
```
cmake ..
 ```

- Build with make:
```
make
 ```

## Running the program ##

- to run marginPhase:
``` marginPhase <BAM_FILE> <REFERENCE_FASTA> <PARAMS> [options] ```

- program OPTIONS:
```
Required arguments:
    BAM is the alignment of reads.  All reads must be aligned to the same contig 
        and be in bam format.
    REFERENCE_FASTA is the reference sequence for the BAM's contig in fasta format.
    PARAMS is the file with marginPhase parameters.


Default options:
    -h --help              : Print this help screen
    -o --outputBase        : Base output identifier [default = "output"]
                               "example" -> "example.sam", "example.vcf"
    -a --logLevel          : Set the log level [default = info]
    -t --tag               : Annotate all output reads with this value for the 
                               'mp' tag

Nucleotide probabilities options:
    -s --singleNuclProbDir : Directory of single nucleotide probabilities files
    -S --onlySNP           : Use only single nucleotide probabilities information,
                               so reads that aren't in SNP dir are discarded

VCF Comparison options:
    -r --referenceVCF      : Reference vcf file for output comparison
    -v --verbose           : Bitmask controlling outputs 
                                 1 - LOG_TRUE_POSITIVES
                                 2 - LOG_FALSE_POSITIVES
                                 4 - LOG_FALSE_NEGATIVES
                               example: 0 -> N/A; 2 -> LFP; 7 -> LTP,LFP,LFN)
```

- Nucleotide Probabilities - this is an alternate input format where reads aligned in a bam have nucleotide alignment posteriors stored in an external location. This option expects the files to be of the form: ${singleNuclProbDir}/${readId}.tsv  If specified, MarginPhase will load the posteriors into its model instead of alignments taken directly from the BAM.  Experimentally, we have found that alignments in this form help ameliorate high error rates found in long reads.

- VCF Comparison - passing in a referenceVCF into the program will do a comparison of the generated VCF and the specified truth set.  Verbosity level will determine what variant classifications will be printed.

- to run tests: ./allTests
(This runs every test. You can comment out ones you don't want to run in allTests.c)

## Notes: ##
- you'll need to have a reference genome available to write out all the files


