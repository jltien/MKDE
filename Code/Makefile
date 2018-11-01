CC = g++ -std=c++11

mkde: mkde.o pairwiseFunctions.o raster.o file_io.o basicFunctions.o
	${CC} mkde.o pairwiseFunctions.o raster.o file_io.o basicFunctions.o -o mkde

mkde.o: mkde.cpp mkde.h
	${CC} -c mkde.cpp

pairwiseFunctions.o: ./Raster/pairwiseFunctions.cpp ./Raster/pairwiseFunctions.h
	${CC} -c ./Raster/pairwiseFunctions.cpp

raster.o: ./Raster/raster.cpp ./Raster/raster.o
	${CC} -c ./Raster/raster.cpp

file_io.o: ./Raster/file_io.cpp ./Raster/file_io.h
	${CC} -c ./Raster/file_io.cpp

basicFunctions.o: ./Raster/basicFunctions.cpp ./Raster/basicFunctions.h
	${CC} -c ./Raster/basicFunctions.cpp

clean:
	rm *.o mkde
