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