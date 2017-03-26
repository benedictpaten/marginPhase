rootPath = ./
include ./include.mk

libSources = impl/*.c
libHeaders = inc/*.h
libTests = tests/*.c
testBin = tests/testBin

all : ${binPath}/stMarginPhaseTest ${binPath}/marginPhase

${binPath}/marginPhase : ${libTests} ${libSources} ${libHeaders} ${basicLibsDependencies} 
	${cxx} ${cflags} -mpopcnt -I inc -I impl -I${libPath} -o ${binPath}/marginPhase marginPhase.c ${libSources} ${basicLibs} 

${binPath}/stMarginPhaseTest : ${libTests} ${libSources} ${libHeaders} ${basicLibsDependencies} 
	${cxx} ${cflags} -mpopcnt -I inc -I impl -I${libPath} -o ${binPath}/stMarginPhaseTest ${libTests} ${libSources} ${basicLibs}  

clean : 
	rm -f *.o
	rm -f ${binPath}/stMarginPhaseTest ${binPath}/marginPhase

test : all
	${binPath}/stMarginPhaseTest
