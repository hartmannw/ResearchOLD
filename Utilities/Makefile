CC=g++
CFLAGS=-Wall -O0
#INCLUDE=-I/people/hartmann/research/src/FileUtilities/

all: StringFunctions.o MatrixFunctions.o

StringFunctions.o: StringFunctions.cc StringFunctions.h
	$(CC) $(CFLAGS) $(INCLUDE) -c StringFunctions.cc

MatrixFunctions.o: MatrixFunctions.cc MatrixFunctions.h
	$(CC) $(CFLAGS) $(INCLUDE) -c MatrixFunctions.cc

clean:
	rm *.o
