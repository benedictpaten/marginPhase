rootPath = ./
include ./include.mk

libSources = impl/*.c
libHeaders = inc/*.h
libTests = tests/*.c
testBin = tests/testBin

#all : ${binPath}/stMarginPhaseTest ${binPath}/marginPhase
all : ${libPath}/marginPhaseLib.a

${binPath}/marginPhase : ${libTests} ${libSources} ${libHeaders} ${basicLibsDependencies} 
	${cxx} ${cflags} -I inc -I impl -I${libPath} -o ${binPath}/marginPhase marginPhase.c ${libTests} ${libSources} ${basicLibs} 

${libPath}/marginPhaseLib.a : 
	${cxx} ${cflags} -I inc -I ${libPath}/ -c ${libSources}
	ar rc marginPhaseLib.a *.o
	ranlib marginPhaseLib.a
	rm *.o
	mv marginPhaseLib.a ${libPath}/
	cp ${libHeaders} ${libPath}

${binPath}/stMarginPhaseTest : ${libTests} ${libSources} ${libHeaders} ${basicLibsDependencies} 
	${cxx} ${cflags} -I inc -I impl -I${libPath} -o ${binPath}/stMarginPhaseTest ${libTests} ${libSources} ${basicLibs}  

clean : 
	rm -f *.o
	rm -f ${binPath}/stMarginPhaseTest ${binPath}/marginPhase

test : all
	${binPath}/stMarginPhaseTest
