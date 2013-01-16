CC=g++
CFLAGS=-Wall -O2
INCLUDE=-I../FileUtilities/ \
				-I../Utilities/ \
				-I../Statistics/ 

all: TestDTW TestMultiBest GeneratePronunciations

GeneratePronunciations: GeneratePronunciations.cc
	$(CC) $(CFLAGS) $(INCLUDE) -o GeneratePronunciations \
	GeneratePronunciations.cc ../FileUtilities/SpeechFeatures.o \
	../Utilities/Matrix.o ../Utilities/StringFunctions.o MultiBestPath.o \
	../Statistics/HmmSet.o ../FileUtilities/ImageIO.o \
	../Statistics/DiagonalGaussian.o ../Statistics/MixtureOfDiagonalGaussians.o \
	../Statistics/PosteriorgramGenerator.o

TestDTW: DynamicTimeWarp.o testdtw.cc
	$(CC) $(CFLAGS) $(INCLUDE) -o testdtw DynamicTimeWarp.o testdtw.cc \
		../FileUtilities/SpeechFeatures.o ../FileUtilities/ImageIO.o

TestMultiBest: MultiBestPath.o
	$(CC) $(CFLAGS) $(INCLUDE) -o testmulti MultiBestPath.o testmultibest.cc \
  	../FileUtilities/SpeechFeatures.o ../Utilities/Matrix.o \
		../FileUtilities/ImageIO.o

DynamicTimeWarp.o: DynamicTimeWarp.cc DynamicTimeWarp.h
	$(CC) $(CFLAGS) $(INCLUDE) -c DynamicTimeWarp.cc

MultiBestPath.o: MultiBestPath.cc MultiBestPath.h
	$(CC) $(CFLAGS) $(INCLUDE) -c MultiBestPath.cc

clean:
	rm *.o testdtw testmulti GeneratePronunciations
