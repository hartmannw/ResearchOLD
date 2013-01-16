CC=g++
CFLAGS=-Wall -O2
INCLUDE=-I../Utilities/ \
    -I../FileUtilities/

all: test_posteriorgram

test_posteriorgram: HmmSet.o MixtureOfDiagonalGaussians.o \
	../Utilities/Matrix.o PosteriorgramGenerator.o test_posteriorgram.cc
	$(CC) $(CFLAGS) $(INCLUDE) -o test_posteriorgram \
    DiagonalGaussian.o MixtureOfDiagonalGaussians.o HmmSet.o \
    ../Utilities/StringFunctions.o PosteriorgramGenerator.o \
    ../FileUtilities/ImageIO.o ../FileUtilities/SpeechFeatures.o \
    test_posteriorgram.cc

DiagonalGaussian.o: DiagonalGaussian.cc DiagonalGaussian.h
	$(CC) $(CFLAGS) $(INCLUDE) -c DiagonalGaussian.cc

MixtureOfDiagonalGaussians.o: MixtureOfDiagonalGaussians.cc \
  MixtureOfDiagonalGaussians.h DiagonalGaussian.o
	$(CC) $(CFLAGS) $(INCLUDE) -c MixtureOfDiagonalGaussians.cc

HiddenMarkovModel.o: HiddenMarkovModel.cc HiddenMarkovModel.h
	$(CC) $(CFLAGS) $(INCLUDE) -c HiddenMarkovModel.cc

HmmSet.o: HmmSet.cc HmmSet.h HiddenMarkovModel.o MixtureOfDiagonalGaussians.o 
	$(CC) $(CFLAGS) $(INCLUDE) -c HmmSet.cc

PosteriorgramGenerator.o: MixtureOfDiagonalGaussians.o ../Utilities/Matrix.o \
	PosteriorgramGenerator.h PosteriorgramGenerator.cc
	$(CC) $(CFLAGS) $(INCLUDE) -c PosteriorgramGenerator.cc

clean:
	rm *.o test_posteriorgram
