#----------------------------------------------------------
# Parameters to be configured by the user

# root for the blas library, for ex. /home/foo/ATLAS/lib/Linux_P4SSE2
#BLASROOT = /home/cpernet/Logiciels/ATLAS_Opteron64/lib/
BLASROOT = /home/cpernet/Logiciels/GotoBLAS

# ATLAS BLAS users : uncomment these lines:
#CXXFLAGS+=-D__LINBOX_HAVE_CBLAS
#LOADLIBES+=-L${BLASROOT} -lcblas -latlas

# GotoBlas BLAS users : uncomment this line
LOADLIBES+=-L${BLASROOT} -lgoto

# Other BLAS users, uncomment this line
#LOADLIBES+=-L${BLASROOT} -lcblas

# architecture parameter for gcc (-march option) for ex. pentium4, athlon ...
#ARCH = -march=pentium3 
#ARCH = -march=pentium4
#ARCH = -march=athlon
#ARCH = -march=opteron 
ARCH = -march=opteron -m64 -mtune=k8

#----------------------------------------------------------

OPTFLAGS+=-O3
OPTFLAGS+=-g -Wall
OPTFLAGS+= ${ARCH}

CXXFLAGS+=${OPTFLAGS} 

INCLUDES+=-I. -I../include -I${BLASROOT}/../../include


#ATLAS users, uncomment this line
#LOADLIBES+= -latlas

CXX=g++ ${INCLUDES}

all: test-fgemm test-invert test-det test-rank test-charpoly test-lqup dense_generator

clean:
	rm -f test-fgemm test-fgemv  test-invert test-det test-rank test-charpoly test-lqup dense_generator

