# Copyright (c) FFLAS-FFPACK
# Written by T. Gautier, jgdumas
# adapted from Givaro.
# ========LICENCE========
# This file is part of the library FFLAS-FFPACK.
#
# FFLAS-FFPACK is free software: you can redistribute it and/or modify
# it under the terms of the  GNU Lesser General Public
# License as published by the Free Software Foundation; either
# version 2.1 of the License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with this library; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
# ========LICENCE========
#/



include ../../Makefile.Rule

SRC=${wildcard *.C}
OBJ=${SRC:.C=.o}

exec: ${OBJ}

%.o: %.C 
	@echo -n -e "   compiling  $* \t..."
	@echo -e "\nCompiling  $* ..." >> ${TARGET_PATH}/compilation.log
	@echo "${CXX} ${CXXFLAGS} $< ${LDLIBS} -o ${BIN_PATH}/$*" >> ${TARGET_PATH}/compilation.log
	@${CXX} ${CXXFLAGS} $< ${LDLIBS} -o ${BIN_PATH}/$* 2>> ${TARGET_PATH}/compilation.log; if test "$$?" = "0"; then echo "done"; else echo "error "; fi

clean:
	rm -rf *.o
