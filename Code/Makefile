CC = g++ -std=c++11

mkde: mkde.o pairwiseFunctions.o raster.o file_io.o ioFunctions.o basicFunctions.o
	${CC} mkde.o pairwiseFunctions.o raster.o file_io.o ioFunctions.o basicFunctions.o -o mkde

mkde.o: mkde.cpp mkde.h ioFunctions.o
	${CC} -c mkde.cpp


pairwiseFunctions.o: ./Raster/pairwiseFunctions.cpp ./Raster/pairwiseFunctions.h
	${CC} -c ./Raster/pairwiseFunctions.cpp

raster.o: ./Raster/raster.cpp ./Raster/raster.h
	${CC} -c ./Raster/raster.cpp

file_io.o: ./Raster/file_io.cpp ./Raster/file_io.h
	${CC} -c ./Raster/file_io.cpp

ioFunctions.o: ./Raster/ioFunctions.cpp ./Raster/ioFunctions.h
	${CC} -c ./Raster/ioFunctions.cpp

basicFunctions.o: ./Raster/basicFunctions.cpp ./Raster/basicFunctions.h
	${CC} -c ./Raster/basicFunctions.cpp

clean:
	rm *.o mkde
