CC=g++
CFLAGS=-Wall -O2
INCLUDE=-I../FileUtilities/ \
				-I../Utilities/

all: TestDTW TestMultiBest

TestDTW: DynamicTimeWarp.o
	$(CC) $(CFLAGS) $(INCLUDE) -o testdtw DynamicTimeWarp.o testdtw.cc \
		../FileUtilities/SpeechFeatures.o ../FileUtilities/ImageIO.o

TestMultiBest: MultiBestPath.o
	$(CC) $(CFLAGS) $(INCLUDE) -o testmulti MultiBestPath.o testmultibest.cc \
  	../FileUtilities/SpeechFeatures.o 

DynamicTimeWarp.o: DynamicTimeWarp.cc DynamicTimeWarp.h
	$(CC) $(CFLAGS) $(INCLUDE) -c DynamicTimeWarp.cc

MultiBestPath.o: MultiBestPath.cc MultiBestPath.h
	$(CC) $(CFLAGS) $(INCLUDE) -c MultiBestPath.cc

clean:
	rm *.o testdtw testmulti
