CC=g++
CFLAGS=-Wall -O0
INCLUDE=-I../../Utilities/ \
    -I../../FileUtilities/ \
		-I../

all: GaussianPosteriorgram AveragePosteriorgram

GaussianPosteriorgram: ../HmmSet.o ../MixtureOfDiagonalGaussians.o GaussianPosteriorgram.cc
	$(CC) $(CFLAGS) $(INCLUDE) -o GaussianPosteriorgram \
    ../DiagonalGaussian.o ../MixtureOfDiagonalGaussians.o ../HmmSet.o \
    ../../Utilities/StringFunctions.o \
    ../../FileUtilities/ImageIO.o ../../FileUtilities/SpeechFeatures.o \
    GaussianPosteriorgram.cc

AveragePosteriorgram: ../HmmSet.o ../MixtureOfDiagonalGaussians.o AveragePosteriorgram.cc
	$(CC) $(CFLAGS) $(INCLUDE) -o AveragePosteriorgram \
    ../DiagonalGaussian.o ../MixtureOfDiagonalGaussians.o ../HmmSet.o \
    ../../Utilities/StringFunctions.o \
    ../../FileUtilities/ImageIO.o ../../FileUtilities/SpeechFeatures.o \
    AveragePosteriorgram.cc

clean:
	rm *.o GaussianPosteriorgram AveragePosteriorgram
