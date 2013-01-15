CC=g++
CFLAGS=-Wall -O0 -g

all: TestMatrix Matrix.o StringFunctions.o MatrixFunctions.o

Matrix.o: Matrix.h Matrix.cc
	$(CC) $(CFLAGS) $(INCLUDE) -c Matrix.cc

StringFunctions.o: StringFunctions.cc StringFunctions.h
	$(CC) $(CFLAGS) $(INCLUDE) -c StringFunctions.cc

MatrixFunctions.o: MatrixFunctions.cc MatrixFunctions.h Matrix.o
	$(CC) $(CFLAGS) $(INCLUDE) -c MatrixFunctions.cc

TestMatrix: testmatrix.cc Matrix.o
	$(CC) $(CFLAGS) $(INCLUDE) -o testmatrix testmatrix.cc Matrix.o

clean:
	rm *.o testmatrix
