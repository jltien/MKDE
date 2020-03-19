1# Makefile for RasterDemo
# By Jeff Tracey, 23 Mar 2018

CC = g++ # icc -ansi # g++
DEBUG_FLAG = #but assign -g for debugging with gdb
LCHECK_FLAG = #assign -lefence to use Electric Fence memory debugger
OPT = -O2 # assign -Wall for warnings; assign - O0 for no optim; -O2 for optim
PROF = #assign -pg for profiling with gprof

RasterDemo : RasterDemo.o raster.o basicFunctions.o file_io.o
	${CC} ${OPT} ${DEBUG_FLAG} ${LCHECK_FLAG} -o RasterDemo RasterDemo.o raster.o basicFunctions.o file_io.o ${PROF}
RasterDemo.o : RasterDemo.cpp raster.h basicFunctions.h file_io.h
	${CC} -c ${OPT} ${DEBUG_FLAG} RasterDemo.cpp ${PROF}
raster.o : raster.cpp raster.h  file_io.h basicDefinitions.h basicFunctions.h
	${CC} -c ${OPT} ${DEBUG_FLAG} raster.cpp ${PROF}
file_io.o : file_io.cpp file_io.h basicDefinitions.h
	${CC} -c ${OPT} ${DEBUG_FLAG} file_io.cpp ${PROF}
basicFunctions.o : basicFunctions.cpp basicFunctions.h basicDefinitions.h
	${CC} -c ${OPT} ${DEBUG_FLAG} basicFunctions.cpp ${PROF}

.PHONY : clean
clean : 
	-rm core.* *.o

