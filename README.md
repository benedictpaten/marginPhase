This is a program for genotyping and haplotyping.

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

- Build htsLib (only required once):
```
cd externalTools/htslib
autoconf
autoheader
./configure
make
cd ../../
```
If the `lzma` and `bz2` packages are not installed, when configuring htslib run
``` ./configure --disable-lzma --disable-bz2 ```
to compile without those packages.

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
``` ./marginPhase <PATH/TO/BAM> <PATH/TO/REFERENCE> [OPTIONS] ```

- program OPTIONS:
```
    -p --params <PATH/TO/JSON>
    -r --referenceVCF <PATH/TO/REFERENCE/VCF>
    -a --logLevel <critical, info, debug [default = info]>
    -o --outputBase <\"example\" -> \"example1.sam\", \"example2.sam\", \"example.vcf\")\n">
    -v --verbose <verbosity bitmask>
```


- to run tests: ./allTests
(This runs every test. You can comment out ones you don't want to run in allTests.c)

## Notes: ##
- you'll need to have a reference genome available to write out all the files


