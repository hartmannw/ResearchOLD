CC=g++
CFLAGS=-Wall -O0
INCLUDE=-I/people/hartmann/research/src/Utilities/ \
    -I/people/hartmann/research/src/FileUtilities/

all: Test GaussianPosteriorgram AveragePosteriorgram

Test: DiagonalGaussian.o MixtureOfDiagonalGaussians.o HmmSet.o SpectralCluster.o KMeans.o
	$(CC) $(CFLAGS) $(INCLUDE) -o test DiagonalGaussian.o \
    MixtureOfDiagonalGaussians.o HmmSet.o ../Utilities/StringFunctions.o \
    ../FileUtilities/ImageIO.o ../Utilities/MatrixFunctions.o \
    SpectralCluster.o KMeans.o test.cc

GaussianPosteriorgram: HmmSet.o MixtureOfDiagonalGaussians.o GaussianPosteriorgram.cc
	$(CC) $(CFLAGS) $(INCLUDE) -o GaussianPosteriorgram \
    DiagonalGaussian.o MixtureOfDiagonalGaussians.o HmmSet.o \
    ../Utilities/StringFunctions.o \
    ../FileUtilities/ImageIO.o ../FileUtilities/SpeechFeatures.o \
    GaussianPosteriorgram.cc

AveragePosteriorgram: HmmSet.o MixtureOfDiagonalGaussians.o AveragePosteriorgram.cc
	$(CC) $(CFLAGS) $(INCLUDE) -o AveragePosteriorgram \
    DiagonalGaussian.o MixtureOfDiagonalGaussians.o HmmSet.o \
    ../Utilities/StringFunctions.o \
    ../FileUtilities/ImageIO.o ../FileUtilities/SpeechFeatures.o \
    AveragePosteriorgram.cc

DiagonalGaussian.o: DiagonalGaussian.cc DiagonalGaussian.h
	$(CC) $(CFLAGS) $(INCLUDE) -c DiagonalGaussian.cc

MixtureOfDiagonalGaussians.o: MixtureOfDiagonalGaussians.cc \
    MixtureOfDiagonalGaussians.h DiagonalGaussian.o
	$(CC) $(CFLAGS) $(INCLUDE) -c MixtureOfDiagonalGaussians.cc

HiddenMarkovModel.o: HiddenMarkovModel.cc HiddenMarkovModel.h
	$(CC) $(CFLAGS) $(INCLUDE) -c HiddenMarkovModel.cc

HmmSet.o: HmmSet.cc HmmSet.h HiddenMarkovModel.o MixtureOfDiagonalGaussians.o 
	$(CC) $(CFLAGS) $(INCLUDE) -c HmmSet.cc

SpectralCluster.o: SpectralCluster.cc SpectralCluster.h
	$(CC) $(CFLAGS) $(INCLUDE) -c SpectralCluster.cc

KMeans.o: KMeans.cc KMeans.h
	$(CC) $(CFLAGS) $(INCLUDE) -c KMeans.cc

clean:
	rm *.o test
