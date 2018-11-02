cmake_minimum_required(VERSION 3.7)
project(marginPhase)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3 -DNDEBUG -std=gnu99 -D_XOPEN_SOURCE=500 -D_POSIX_C_SOURCE=200112L -mpopcnt ")
#set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3 -DNDEBUG")

#include(FindOpenMP)
#if(OPENMP_FOUND)
#    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
#    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
#    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
#endif()

set(HTSLIB_HEADERS
        externalTools/htslib
        externalTools/htslib/htslib
        externalTools/htslib/cram
        externalTools/htslib/os
        )
set(SONLIB_HEADERS
        externalTools/sonLib/C/inc
        externalTools/sonLib/C/impl
        externalTools/sonLib/C/tests
        externalTools/sonLib/externalTools/cutest
        externalTools/sonLib/externalTools/quicktree_1.1/include
        )
set(CPECAN_HEADERS
        externalTools/cPecan/inc
        externalTools/cPecan/tests
        externalTools/cPecan/externalTools/lastz-distribu-1.03.54/src
        )

include_directories(inc)
include_directories(${HTSLIB_HEADERS})
include_directories(${SONLIB_HEADERS})
include_directories(${CPECAN_HEADERS})
include_directories(externalTools/jsmn)

set(HTSLIB_SOURCES
        externalTools/htslib/cram/cram_codecs.c
        externalTools/htslib/cram/cram_decode.c
        externalTools/htslib/cram/cram_encode.c
        externalTools/htslib/cram/cram_external.c
        externalTools/htslib/cram/cram_index.c
        externalTools/htslib/cram/cram_io.c
        externalTools/htslib/cram/cram_samtools.c
        externalTools/htslib/cram/cram_stats.c
        externalTools/htslib/cram/files.c
        externalTools/htslib/cram/mFILE.c
        externalTools/htslib/cram/open_trace_file.c
        externalTools/htslib/cram/pooled_alloc.c
        externalTools/htslib/cram/rANS_static.c
        externalTools/htslib/cram/sam_header.c
        externalTools/htslib/cram/string_alloc.c
        externalTools/htslib/os/rand.c
        externalTools/htslib/test/fieldarith.c
        externalTools/htslib/test/hfile.c
        externalTools/htslib/test/hts_endian.c
        externalTools/htslib/test/sam.c
        externalTools/htslib/test/test-bcf-sr.c
        externalTools/htslib/test/test-regidx.c
        externalTools/htslib/test/test-vcf-api.c
        externalTools/htslib/test/test-vcf-sweep.c
        externalTools/htslib/test/test_bgzf.c
        externalTools/htslib/test/test_view.c
        externalTools/htslib/test/thrash_threads1.c
        externalTools/htslib/test/thrash_threads2.c
        externalTools/htslib/test/thrash_threads3.c
        externalTools/htslib/test/thrash_threads4.c
        externalTools/htslib/test/thrash_threads5.c
        externalTools/htslib/test/thrash_threads6.c
        externalTools/htslib/bcf_sr_sort.c
        externalTools/htslib/bgzf.c
        externalTools/htslib/bgzip.c
        externalTools/htslib/errmod.c
        externalTools/htslib/faidx.c
        externalTools/htslib/hfile.c
        externalTools/htslib/hfile_gcs.c
        externalTools/htslib/hfile_libcurl.c
        externalTools/htslib/hfile_net.c
        #externalTools/htslib/hfile_s3.c
        externalTools/htslib/hts.c
        externalTools/htslib/htsfile.c
        externalTools/htslib/kfunc.c
        externalTools/htslib/knetfile.c
        externalTools/htslib/kstring.c
        externalTools/htslib/md5.c
        externalTools/htslib/multipart.c
        #externalTools/htslib/plugin.c
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
        externalTools/sonLib/C/impl/avl.c
        externalTools/sonLib/C/impl/bioioC.c
        externalTools/sonLib/C/impl/commonC.c
        externalTools/sonLib/C/impl/fastCMaths.c
        externalTools/sonLib/C/impl/hashTableC.c
        externalTools/sonLib/C/impl/hashTableC_itr.c
        externalTools/sonLib/C/impl/lz4.c
        externalTools/sonLib/C/impl/lz4hc.c
        externalTools/sonLib/C/impl/pairwiseAlignment.c
        externalTools/sonLib/C/impl/sonLibAlign.c
        externalTools/sonLib/C/impl/sonLibCache.c
        externalTools/sonLib/C/impl/sonLibCommon.c
        externalTools/sonLib/C/impl/sonLibCompression.c
        externalTools/sonLib/C/impl/sonLibConnectivity.c
        externalTools/sonLib/C/impl/sonLibEulerTour.c
        externalTools/sonLib/C/impl/sonLibExcept.c
        externalTools/sonLib/C/impl/sonLibFile.c
        externalTools/sonLib/C/impl/sonLibHash.c
        externalTools/sonLib/C/impl/sonLibKVDatabase.c
        externalTools/sonLib/C/impl/sonLibKVDatabase_BigRecordFile.c
        externalTools/sonLib/C/impl/sonLibKVDatabase_KyotoTycoon.cpp
        externalTools/sonLib/C/impl/sonLibKVDatabase_MySql.c
        externalTools/sonLib/C/impl/sonLibKVDatabase_TokyoCabinet.c
        externalTools/sonLib/C/impl/sonLibKVDatabaseConf.c
        externalTools/sonLib/C/impl/sonLibList.c
        externalTools/sonLib/C/impl/sonLibMath.c
        externalTools/sonLib/C/impl/sonLibNaiveConnectivity.c
        externalTools/sonLib/C/impl/sonLibRandom.c
        externalTools/sonLib/C/impl/sonLibSet.c
        externalTools/sonLib/C/impl/sonLibSortedSet.c
        externalTools/sonLib/C/impl/sonLibString.c
        externalTools/sonLib/C/impl/sonLibTreap.c
        externalTools/sonLib/C/impl/sonLibTree.c
        externalTools/sonLib/C/impl/sonLibTuples.c
        externalTools/sonLib/C/impl/stGraph.c
        externalTools/sonLib/C/impl/stMatrix.c
        externalTools/sonLib/C/impl/stPhylogeny.c
        externalTools/sonLib/C/impl/stPosetAlignment.c
        externalTools/sonLib/C/impl/stSafeC.c
        externalTools/sonLib/C/impl/stThreadPool.c
        externalTools/sonLib/C/impl/stUnionFind.c
        externalTools/sonLib/C/tests/allTests.c
        externalTools/sonLib/C/tests/cigarsTest.c
        externalTools/sonLib/C/tests/fastaCTest.c
        #externalTools/sonLib/C/tests/kt_connect_test.cpp
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
        externalTools/sonLib/externalTools/cutest/AllTests.c
        externalTools/sonLib/externalTools/cutest/CuTest.c
        externalTools/sonLib/externalTools/cutest/CuTestTest.c
        externalTools/sonLib/externalTools/quicktree_1.1/src/align.c
        externalTools/sonLib/externalTools/quicktree_1.1/src/aln2dist.c
        externalTools/sonLib/externalTools/quicktree_1.1/src/aln2tree.c
        externalTools/sonLib/externalTools/quicktree_1.1/src/buildtree.c
        externalTools/sonLib/externalTools/quicktree_1.1/src/cluster.c
        externalTools/sonLib/externalTools/quicktree_1.1/src/dist2tree.c
        externalTools/sonLib/externalTools/quicktree_1.1/src/distancemat.c
        externalTools/sonLib/externalTools/quicktree_1.1/src/distancemat_merops.c
        externalTools/sonLib/externalTools/quicktree_1.1/src/options.c
        externalTools/sonLib/externalTools/quicktree_1.1/src/quicktree.c
        externalTools/sonLib/externalTools/quicktree_1.1/src/sequence.c
        externalTools/sonLib/externalTools/quicktree_1.1/src/tree.c
        externalTools/sonLib/externalTools/quicktree_1.1/src/util.c
        )

set(LASTZ_SOURCES
        externalTools/cPecan/externalTools/lastz-distrib-1.03.54/src/align_diffs.c
        externalTools/cPecan/externalTools/lastz-distrib-1.03.54/src/axt.c
        externalTools/cPecan/externalTools/lastz-distrib-1.03.54/src/capsule.c
        externalTools/cPecan/externalTools/lastz-distrib-1.03.54/src/chain.c
        externalTools/cPecan/externalTools/lastz-distrib-1.03.54/src/cigar.c
        externalTools/cPecan/externalTools/lastz-distrib-1.03.54/src/continuity_dist.c
        externalTools/cPecan/externalTools/lastz-distrib-1.03.54/src/coverage_dist.c
        externalTools/cPecan/externalTools/lastz-distrib-1.03.54/src/diag_hash.c
        externalTools/cPecan/externalTools/lastz-distrib-1.03.54/src/dna_utilities.c
        externalTools/cPecan/externalTools/lastz-distrib-1.03.54/src/edit_script.c
        externalTools/cPecan/externalTools/lastz-distrib-1.03.54/src/gapped_extend.c
        externalTools/cPecan/externalTools/lastz-distrib-1.03.54/src/genpaf.c
        externalTools/cPecan/externalTools/lastz-distrib-1.03.54/src/gfa.c
        externalTools/cPecan/externalTools/lastz-distrib-1.03.54/src/identity_dist.c
        externalTools/cPecan/externalTools/lastz-distrib-1.03.54/src/infer_scores.c
        externalTools/cPecan/externalTools/lastz-distrib-1.03.54/src/lastz.c
        externalTools/cPecan/externalTools/lastz-distrib-1.03.54/src/lav.c
        externalTools/cPecan/externalTools/lastz-distrib-1.03.54/src/maf.c
        externalTools/cPecan/externalTools/lastz-distrib-1.03.54/src/masking.c
        externalTools/cPecan/externalTools/lastz-distrib-1.03.54/src/output.c
        externalTools/cPecan/externalTools/lastz-distrib-1.03.54/src/pos_table.c
        externalTools/cPecan/externalTools/lastz-distrib-1.03.54/src/quantum.c
        externalTools/cPecan/externalTools/lastz-distrib-1.03.54/src/sam.c
        externalTools/cPecan/externalTools/lastz-distrib-1.03.54/src/seed_search.c
        externalTools/cPecan/externalTools/lastz-distrib-1.03.54/src/seeds.c
        externalTools/cPecan/externalTools/lastz-distrib-1.03.54/src/segment.c
        externalTools/cPecan/externalTools/lastz-distrib-1.03.54/src/sequences.c
        externalTools/cPecan/externalTools/lastz-distrib-1.03.54/src/text_align.c
        externalTools/cPecan/externalTools/lastz-distrib-1.03.54/src/tweener.c
        externalTools/cPecan/externalTools/lastz-distrib-1.03.54/src/utilities.c
        )
       
set(CPECAN_SOURCES
        externalTools/cPecan/impl/multipleAligner.c
        externalTools/cPecan/impl/pairwiseAligner.c
        externalTools/cPecan/impl/stateMachine.c
        externalTools/cPecan/impl/randomSequences.c
        externalTools/cPecan/tests/allTests.c
        externalTools/cPecan/tests/multipleAlignerTest.c
        externalTools/cPecan/tests/pairwiseAlignerLongTest.c
        externalTools/cPecan/tests/pairwiseAlignerTest.c
        )

set(SOURCE_FILES
        impl/chunker.c
        impl/column.c
        impl/coordination.c
        impl/emissions.c
        impl/genomeFragment.c
        impl/hmm.c
        impl/mergeColumn.c
        impl/outputWriter.c
        impl/parser.c
        impl/partitions.c
        impl/polisher.c
        impl/profileSeq.c
        impl/referencePriorProbs.c
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

find_library(PTHREADLIB pthread)
if(${PTHREADLIB} STREQUAL "PTHREADLIB-NOTFOUND")
    message(FATAL_ERROR "Couldn't find the 'pthread' library")
endif()

find_library(CURLLIB curl)
if(${CURLLIB} STREQUAL "CURLLIB-NOTFOUND")
    message(FATAL_ERROR "Couldn't find the 'curl' library")
endif()

add_custom_command(
        OUTPUT "${PROJECT_SOURCE_DIR}/externalTools/htslib/config.h"
        COMMAND autoconf
        COMMAND autoheader
        COMMAND ./configure --disable-lzma --disable-s3 --disable-plugins
        COMMAND make
        WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/externalTools/htslib/
)
add_custom_target(HTSLIB_CONFIGURED DEPENDS "${PROJECT_SOURCE_DIR}/externalTools/htslib/config.h")

add_library(son ${SONLIB_SOURCES})
add_library(lastz ${LASTZ_SOURCES})
add_library(cpecan ${CPECAN_SOURCES})
add_dependencies(cpecan lastz)
add_library(hts ${HTSLIB_SOURCES})
add_dependencies(hts HTSLIB_CONFIGURED)
add_library(jsmn externalTools/jsmn/jsmn.c)

target_compile_options(lastz PUBLIC -c -O3 -Wall -Wextra -D_FILE_OFFSET_BITS=64 -D_LARGEFILE_SOURCE -DVERSION_MAJOR=\"1\" -DVERSION_MINOR=\"03\" -DVERSION_SUBMINOR=\"54\" -DREVISION_DATE=\"20140128\" -DSUBVERSION_REV=\"1827:1830\" -Dscore_type='D')

target_link_libraries(hts pthread bz2 z curl)

add_executable(marginPhase marginPhase.c ${SOURCE_FILES})
target_link_libraries(marginPhase cpecan son hts jsmn)

#TODO uncomment after compilation errors are fixed
add_executable(marginPolish marginPolish.c ${SOURCE_FILES})
target_link_libraries(marginPolish cpecan son hts jsmn)

add_executable(vcfCompare vcfCompare.c ${SOURCE_FILES})
target_link_libraries(vcfCompare cpecan son hts jsmn)

enable_testing()
add_executable(allTests tests/allTests.c tests/marginPhaseTest.c tests/stRPHmmTest.c tests/parserTest.c tests/polisherTest.c ${SOURCE_FILES})
target_link_libraries(allTests cpecan hts cpecan son jsmn)
add_test(stRPHmm allTests)
