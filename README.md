# MarginPolish #

MarginPolish is a graph-based assembly polisher that takes advantage of multiple probable alignment paths.  It takes as input a FASTA assembly and BAM, and produces a polished FASTA assembly.  Secondarily, it produces image files used by the HELEN polisher.

## Dependencies ##
cmake version 3.7 (or higher):
```
wget https://cmake.org/files/v3.7/cmake-3.7.2-Linux-x86_64.sh && mkdir /opt/cmake && sh cmake-3.7.2-Linux-x86_64.sh --prefix=/opt/cmake --skip-license && ln -s /opt/cmake/bin/cmake /usr/local/bin/cmake
```

If you're on Ubuntu:
```
apt-get -y install git make gcc g++ autoconf bzip2 lzma-dev zlib1g-dev libcurl4-openssl-dev libpthread-stubs0-dev libbz2-dev liblzma-dev libhdf5-dev
```

- Check out the repository:
```
git clone https://github.com/UCSC-nanopore-cgl/marginPolish.git
```

- Check out submodules:
```
cd marginPolish
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

## Running ##

### Running MarginPolish ###


- to run marginPolish:
``` marginPolish <BAM_FILE> <ASSEMBLY_FASTA> <PARAMS> [options] ```

- program OPTIONS:
```
Polishes the ASSEMBLY_FASTA using alignments in BAM_FILE.

Required arguments:
    BAM_FILE is the alignment of reads to the assembly (or reference).
    ASSEMBLY_FASTA is the reference sequence BAM file in fasta format.
    PARAMS is the file with marginPolish parameters.

Default options:
    -h --help                : Print this help screen
    -a --logLevel            : Set the log level [default = info]
    -t --threads             : Set number of concurrent threads [default = 1]
    -o --outputBase          : Name to use for output files [default = 'output']
    -r --region              : If set, will only compute for given chromosomal region.
                                 Format: chr:start_pos-end_pos (chr3:2000-3000).

HELEN feature generation options:
    -f --outputFeatureType   : output features of chunks for HELEN.  Valid types:
                                 simpleWeight:    weighted likelihood from POA nodes (non-RLE)
                                 rleWeight:       weighted likelihood from POA nodes (RLE)
                                 nuclAndRlWeight: weighted likelihood, split into nucleotide and run length
                                 splitRleWeight:  weighted likelihood, with run lengths split into chunks
    -L --splitRleWeightMaxRL : max run length (for 'splitRleWeight' type only) [default = 10]
    -u --trueReferenceBam    : true reference aligned to ASSEMBLY_FASTA, for HELEN
                               features.  Setting this parameter will include labels
                               in output.
    -5 --hdf5Only            : only output H5 feature files.  Default behavior is to output
                               h5, tsv, and fa for each chunk.

Miscellaneous supplementary output options:
    -i --outputRepeatCounts  : Output base to write out the repeat counts [default = NULL]
    -j --outputPoaTsv        : Output base to write out the poa as TSV file [default = NULL]

```


#### Sample Execution

```./marginPolish ../tests/NA12878.np.chr3.5kb.bam ../tests/hg19.chr3.9mb.fa ../params/allParams.np.human.json -o example_out```


### MarginPhase ###

MarginPhase is a program for simultaneous haplotyping and genotyping with long reads.  It relies on much of the same infrastructure as marginPolish.  We refer you to the marginPhase repo (from which this is forked) at:

https://github.com/benedictpaten/marginPhase


### Tests ###

After building run the 'allTests' executable in your build directory.  This runs every test. You can comment out ones you don't want to run in allTests.c
