# Makefile for RasterDemo
# By Jeff Tracey, 23 Mar 2018

CC = g++ # icc -ansi # g++
DEBUG_FLAG = #but assign -g for debugging with gdb
LCHECK_FLAG = #assign -lefence to use Electric Fence memory debugger
OPT = -O2 # assign -Wall for warnings; assign - O0 for no optim; -O2 for optim
PROF = #assign -pg for profiling with gprof

mkde : mkde.o raster.o basicFunctions.o file_io.o
	${CC} ${OPT} ${DEBUG_FLAG} ${LCHECK_FLAG} -o mkde.o RasterDemo.o raster.o basicFunctions.o file_io.o ${PROF}
mkde.o : mkde.cpp AnimalData.h ./Raster/raster.h ./Raster/basicFunctions.h ./Raster/file_io.h
	${CC} -c ${OPT} ${DEBUG_FLAG} mkde.cpp ${PROF}
raster.o : ./Raster/raster.cpp ./Raster/raster.h  ./Raster/file_io.h ./Raster/basicDefinitions.h ./Raster/basicFunctions.h
	${CC} -c ${OPT} ${DEBUG_FLAG} ./Raster/raster.cpp ${PROF}
file_io.o : file_io.cpp file_io.h ./Raster/basicDefinitions.h
	${CC} -c ${OPT} ${DEBUG_FLAG} ./Raster/file_io.cpp ${PROF}
basicFunctions.o : basicFunctions.cpp ./Raster/basicFunctions.h ./Raster/basicDefinitions.h
	${CC} -c ${OPT} ${DEBUG_FLAG} ./Raster/basicFunctions.cpp ${PROF}

.PHONY : clean
clean : 
	-rm core.* *.o

