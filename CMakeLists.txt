cmake_minimum_required(VERSION 3.7)
project(marginPhase)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=gnu99 -O0 --pedantic -g -fno-inline -fpic")

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
include_directories(${SONLIB_HEADERS})
include_directories(${CPECAN_HEADERS})
include_directories(externalTools/jsmn)

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
#        externalTools/sonLib/C/impl/sonLibCompression.c
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
#        externalTools/sonLib/C/tests/allTests.c
        externalTools/sonLib/C/tests/cigarsTest.c
        externalTools/sonLib/C/tests/fastaCTest.c
        #externalTools/sonLib/C/tests/kt_connect_test.cpp
        externalTools/sonLib/C/tests/kvDatabaseTest.c
        externalTools/sonLib/C/tests/kvDatabaseTestCommon.c
        externalTools/sonLib/C/tests/sonLibAlignTest.c
        externalTools/sonLib/C/tests/sonLibCacheTest.c
        externalTools/sonLib/C/tests/sonLibCommonTest.c
#        externalTools/sonLib/C/tests/sonLibCompressionTest.c
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
#        impl/outputWriter.c
#        impl/vcfTools.c
        impl/polisher.c
        impl/callConsensus.c
        )

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
add_library(cpecan ${CPECAN_SOURCES})
add_library(jsmn externalTools/jsmn/jsmn.c)

enable_testing()
add_executable(allTests tests/stRPHmmTest.c tests/polisherTest.c ${SOURCE_FILES})
target_link_libraries(allTests cpecan son jsmn)
add_test(stRPHmm allTests)

add_executable(simpleConsensusTest tests/simpleConsensusTest.c ${SOURCE_FILES})
target_link_libraries(simpleConsensusTest cpecan son jsmn)

add_library(callConsensus SHARED impl/callConsensus.c ${SOURCE_FILES})
target_link_libraries(callConsensus cpecan son jsmn)
