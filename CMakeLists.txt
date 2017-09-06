cmake_minimum_required(VERSION 3.7)
project(marginPhase)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=gnu99 -D_XOPEN_SOURCE -D_POSIX_C_SOURCE=200112L -mpopcnt -O3")
# add DNDEBUG flag back when done

include(FindOpenMP)
if(OPENMP_FOUND)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

set(HTSLIB_HEADERS
        externalTools/htslib
        externalTools/htslib/htslib
        externalTools/htslib/cram
        )
set(SONLIB_HEADERS
#        externalTools/sonLib/lib   # old header file location for sonlib
        externalTools/sonLib/C/inc
        externalTools/sonLib/C/impl
        externalTools/sonLib/C/tests
        externalTools/sonLib/externalTools/cutest
        externalTools/sonLib/externalTools/quicktree_1.1/include
        )

include_directories(inc)
include_directories(${HTSLIB_HEADERS})
include_directories(${SONLIB_HEADERS})
include_directories(externalTools/jsmn)

file(GLOB htslib_cram_SRC "externalTools/htslib/cram/*.c")
file(GLOB htslib_test_SRC "externalTools/htslib/test/*.c")
file(GLOB sonlib_impl_SRC "externalTools/sonLib/C/impl/*.c")
file(GLOB sonlib_cutest_SRC "externalTools/sonLib/externalTools/cutest/*.c")
file(GLOB sonlib_quicktree_SRC "externalTools/sonLib/externalTools/quicktree_1.1/src/*.c")

set(HTSLIB_SOURCES
        ${htslib_cram_SRC}
        ${htslib_test_SRC}
        externalTools/htslib/bcf_sr_sort.c
        externalTools/htslib/bgzf.c
        externalTools/htslib/bgzip.c
        externalTools/htslib/errmod.c
        externalTools/htslib/faidx.c
        externalTools/htslib/hfile.c
        externalTools/htslib/hfile_gcs.c
        externalTools/htslib/hfile_libcurl.c
        externalTools/htslib/hfile_net.c
        externalTools/htslib/hfile_s3.c
        externalTools/htslib/hts.c
        externalTools/htslib/htsfile.c
        externalTools/htslib/kfunc.c
        externalTools/htslib/knetfile.c
        externalTools/htslib/kstring.c
        externalTools/htslib/md5.c
        externalTools/htslib/multipart.c
        # externalTools/htslib/plugin.c    # build problems including this file
        externalTools/htslib/probaln.c
        externalTools/htslib/realn.c
        externalTools/htslib/regidx.c
        externalTools/htslib/sam.c
        externalTools/htslib/synced_bcf_reader.c
        externalTools/htslib/tabix.c
        externalTools/htslib/tbx.c
        externalTools/htslib/textutils.c
        externalTools/htslib/thread_pool.c
        externalTools/htslib/vcf.c
        externalTools/htslib/vcf_sweep.c
        externalTools/htslib/vcfutils.c
        )

set(SONLIB_SOURCES
        ${sonlib_impl_SRC}
        ${sonlib_cutest_SRC}
        ${sonlib_quicktree_SRC}
        externalTools/sonLib/C/tests/allTests.c
        externalTools/sonLib/C/tests/cigarsTest.c
        externalTools/sonLib/C/tests/fastaCTest.c
        # externalTools/sonLib/C/tests/kt_connect_test.cpp    # build problems including this file
        externalTools/sonLib/C/tests/kvDatabaseTest.c
        externalTools/sonLib/C/tests/kvDatabaseTestCommon.c
        externalTools/sonLib/C/tests/sonLibAlignTest.c
        externalTools/sonLib/C/tests/sonLibCacheTest.c
        externalTools/sonLib/C/tests/sonLibCommonTest.c
        externalTools/sonLib/C/tests/sonLibCompressionTest.c
        externalTools/sonLib/C/tests/sonLibConnectivityTests.c
        externalTools/sonLib/C/tests/sonLibDoubleTuplesTest.c
        externalTools/sonLib/C/tests/sonLibEdgeContainerTests.c
        externalTools/sonLib/C/tests/sonLibEulerTest.c
        externalTools/sonLib/C/tests/sonLibExceptTest.c
        externalTools/sonLib/C/tests/sonLibFileTest.c
        externalTools/sonLib/C/tests/sonLibGraphTest.c
        externalTools/sonLib/C/tests/sonLibHashTest.c
        externalTools/sonLib/C/tests/sonLibIntTuplesTest.c
        externalTools/sonLib/C/tests/sonLibListTest.c
        externalTools/sonLib/C/tests/sonLibPosetAlignmentTest.c
        externalTools/sonLib/C/tests/sonLibRandomTest.c
        externalTools/sonLib/C/tests/sonLibSetTest.c
        externalTools/sonLib/C/tests/sonLibSortedSetTest.c
        externalTools/sonLib/C/tests/sonLibStringTest.c
        externalTools/sonLib/C/tests/sonLibTreapTest.c
        externalTools/sonLib/C/tests/sonLibTreeTest.c
        externalTools/sonLib/C/tests/sonLibTuplesTest.c
        externalTools/sonLib/C/tests/stMatrixTest.c
        externalTools/sonLib/C/tests/stPhylogenyTest.c
        externalTools/sonLib/C/tests/stThreadPoolTest.c
        externalTools/sonLib/C/tests/stUnionFindTest.c
        )

set(SOURCE_FILES
        impl/column.c
        impl/coordination.c
        impl/emissions.c
        impl/genomeFragment.c
        impl/hmm.c
        impl/mergeColumn.c
        impl/partitions.c
        impl/profileSeq.c
        impl/referencePriorProbs.c
        impl/parser.c
        impl/outputWriter.c
        impl/vcfTools.c
        )


# Dependencies for htslib
find_library(ZLIB NAMES z)
if(${ZLIB} STREQUAL "ZLIB-NOTFOUND")
    message(FATAL_ERROR "Couldn't find the 'z' library")
endif()
find_library(BZ2LIB bz2)
if(${BZ2LIB} STREQUAL "BZ2LIB-NOTFOUND")
    message(FATAL_ERROR "Couldn't find the 'bz2' library")
endif()
find_library(LZMALIB lzma)
if(${LZMALIB} STREQUAL "LZMALIB-NOTFOUND")
    message(FATAL_ERROR "Couldn't find the 'lzma' library")
endif()
find_library(CURLLIB curl)
if(${CURLLIB} STREQUAL "CURLLIB-NOTFOUND")
    message(FATAL_ERROR "Couldn't find the 'curl' library")
endif()
find_library(CRYPTOLIB crypto)
if(${CRYPTOLIB} STREQUAL "CRYPTOLIB-NOTFOUND")
    message(FATAL_ERROR "Couldn't find the 'crypto' library")
endif()
find_library(PTHREADLIB pthread)
if(${PTHREADLIB} STREQUAL "PTHREADLIB-NOTFOUND")
    message(FATAL_ERROR "Couldn't find the 'pthread' library")
endif()

add_library(son ${SONLIB_SOURCES})
add_library(hts ${HTSLIB_SOURCES})
add_library(jsmn externalTools/jsmn/jsmn.c)

target_link_libraries(hts pthread crypto bz2 lzma z curl)

add_executable(marginPhase marginPhase.c ${SOURCE_FILES})
target_link_libraries(marginPhase son hts jsmn)

add_executable(vcfCompare vcfCompare.c ${SOURCE_FILES})
target_link_libraries(vcfCompare son hts jsmn)


enable_testing()
add_executable(allTests tests/allTests.c tests/marginPhaseTest.c tests/stRPHmmTest.c tests/parserTest.c ${SOURCE_FILES})
target_link_libraries(allTests hts son jsmn)
add_test(stRPHmm allTests)