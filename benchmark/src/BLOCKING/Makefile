OPTFLAGS= -O7
OPTFLAGS+=-funroll-all-loops
OPTFLAGS+=-frerun-loop-opt
OPTFLAGS+=-fexpensive-optimizations -felide-constructors
OPTFLAGS+=-ffast-math 
OPTFLAGS+=-fforce-addr -fforce-mem -fstrength-reduce
OPTFLAGS+=-fstrict-aliasing
#OPTFLAGS+=-mtune=pentium -march=pentium

tblockmat-%:tblockmat.C
	${CXX} ${OPTFLAGS} ${OPT} -o $@ -DNDIM=$* $< ${CXXFLAGS} ${LDFLAGS} ${LOADLIBES}
tblockmatdouble-%:tblockmatdouble.C
	${CXX} ${OPTFLAGS} ${OPT} -o $@ -DNDIM=$* $< ${CXXFLAGS} ${LDFLAGS} ${LOADLIBES}
