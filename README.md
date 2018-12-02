# MarginPhase #

MarginPhase is a program for simultaneous haplotyping and genotyping.  Included is a second program MarginPolish, which polishes long read assemblies.

## Dependencies ##
cmake version 3.7 (or higher):
```
wget https://cmake.org/files/v3.7/cmake-3.7.2-Linux-x86_64.sh && mkdir /opt/cmake && sh cmake-3.7.2-Linux-x86_64.sh --prefix=/opt/cmake --skip-license && ln -s /opt/cmake/bin/cmake /usr/local/bin/cmake
```

If you're on Ubuntu:
```
apt-get -y install git make gcc g++ autoconf bzip2 lzma-dev zlib1g-dev libcurl4-openssl-dev libpthread-stubs0-dev libbz2-dev liblzma-dev
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

## Running ##

### Running MarginPhase ###

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

- you'll need to have a reference genome available to write out all the files

#### Sample Execution

```./marginPhase ../tests/NA12878.np.chr3.5kb.bam ../tests/hg19.chr3.9mb.fa ../params/params.nanopore.json -o example_out```


### Running MarginPolish ###

MarginPolish is a polishing utility for long read data.  It takes an assembled fasta file and a BAM file of reads aligned to the assembly (and a parameters file).  It outputs a polished FASTA, which will eventually include a BAM with the reads re-alinged to the polished assembly.

- to run marginPolish:
``` marginPhase <BAM_FILE> <ASSEMBLY_FASTA> <PARAMS> [options] ```

- program OPTIONS:
```

Required arguments:
    BAM_FILE is the alignment of reads to the assembly (or reference)
    ASSEMBLY_FASTA is the reference sequence BAM file in fasta format.
    PARAMS is the file with marginPolish parameters.

Default options:
    -h --help              : Print this help screen
    -a --logLevel          : Set the log level [default = info]
    -o --outputBase        : Name to use for output files [default = output]
    -r --region            : If set, will only compute for given chromosomal region.
                               Format: chr:start_pos-end_pos (chr3:2000-3000).

```


#### Sample Execution

```./marginPolish ../tests/NA12878.np.chr3.5kb.bam ../tests/hg19.chr3.9mb.fa ../params/params.nanopore.json -o example_out```

#### Evaluation Workflow

The 'assembly_experiment.py' script in the scripts/ directory is a simple workflow for evaluating the effectiveness of marginPolish on an assembly.  It requires an assembly fasta, alignment (to the assembly), true reference, and set of parameters, as well as an existing installation of minimap2 which accepts the '--eqx' parameter (the latest version does this).  

The script will run marginPolish on the assembly and alignment.  Then it will break the polished assembly into chunks (size is configurable), align them to the true reference, and compute statistics on the alignment.  It also performs this on the original assembly for comparison.  

All intermediate output (as well as the summarized results at the end) are stored, and the program will reuse these if present for analysis.

```
usage: Runs marginPolish on assembly, compares polished assembly to true reference
       [-h] --true_reference TRUE_REFERENCE [--assembly_fasta ASSEMBLY_FASTA]
       [--assembly_bam ASSEMBLY_BAM] [--polisher_params POLISHER_PARAMS]
       [--output_name OUTPUT_NAME] [--read_size READ_SIZE]
       [--marginPolish_invocation MARGINPOLISH_INVOCATION]
       [--minimap2_invocation MINIMAP2_INVOCATION]

```

### Tests ###

After building run the 'allTests' executable in your build directory.  This runs every test. You can comment out ones you don't want to run in allTests.c
