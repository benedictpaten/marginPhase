#Modify this variable to set the location of sonLib
sonLibRootPath ?= ${rootPath}/../sonLib
sonLibPath=${sonLibRootPath}/lib
#Use sonLib bin and lib dirs
binPath=${sonLibRootPath}/bin
libPath=${sonLibPath}

include  ${sonLibRootPath}/include.mk

# Set the compiler here / use of OpenMP
cxx = gcc-6 -fopenmp -DHAVE_OPEN_MP=1
# use export OMP_NUM_THREADS=FOO in shell to define number of threads to use

basicLibs = ${sonLibPath}/sonLib.a ${sonLibPath}/cuTest.a ${dblibs}
basicLibsDependencies = ${sonLibPath}/sonLib.a ${sonLibPath}/cuTest.a 

