# Makefile for TestMain
# By Jeff Tracey, 12/1/02; Copyright Jeff A. Tracey 2002 All rights reserved.
# DO NOT SEE WHY basicFunctions.* SHOULD BE NEEDED AT ALL (BUT IS INCLUDED IN GeometricPrimatives.h).
# mainTODO.* IS A PLACE HOLDER FOR A TESTING PROGRAM

CC = g++ #icc -ansi # g++
OBJS = GeometricPrimatives.o file_io.o raster.o 
HDRS = GeometricPrimatives.h file_io.h raster.h
DEBUG_FLAG = #but assign -g for debugging with gdb
LCHECK_FLAG = #assign -lefence to use Electric Fence memory debugger
OPT = -O2 # assign -Wall for warnings; assign - O0 for no optim; -O2 for optim
PROF = #assign -pg for profiling with gprof

RasterTestMain : mainTODO.o ${OBJS}
	${CC} ${OPT} ${DEBUG_FLAG} ${LCHECK_FLAG} -o DataX DataXmain.o ${OBJS} ${PROF}
mainTODO.o : mainTODO.cpp ${HDRS}
	${CC} -c ${OPT} ${DEBUG_FLAG} DataXmain.cpp ${PROF}
raster.o : raster.cpp raster.h
	${CC} -c ${OPT} ${DEBUG_FLAG} raster_io.cpp ${PROF}
file_io.o : file_io.cpp file_io.h
	${CC} -c ${OPT} ${DEBUG_FLAG} file_io.cpp ${PROF}
GeometricPrimatives.o : GeometricPrimatives.cpp GeometricPrimatives.h
	${CC} -c ${OPT} ${DEBUG_FLAG} GeometricPrimatives.cpp ${PROF}

.PHONY : clean
clean : 
	-rm core.* *.o

