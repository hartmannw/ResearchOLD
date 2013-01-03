CC=g++
CFLAGS=-Wall -O2
INCLUDE=-I../FileUtilities/

all: TestDTW IdentifyWordAlignments CreateBestSegments FindSilenceExemplars \
		IdentifySilenceLocations

FindSilenceExemplars: DynamicTimeWarp.o FindSilenceExemplars.cc
	$(CC) $(CFLAGS) $(INCLUDE) -o FindSilenceExemplars DynamicTimeWarp.o \
		FindSilenceExemplars.cc ../FileUtilities/SpeechFeatures.o \
		../FileUtilities/ImageIO.o

IdentifySilenceLocations: DynamicTimeWarp.o IdentifySilenceLocations.cc
	$(CC) $(CFLAGS) $(INCLUDE) -o IdentifySilenceLocations DynamicTimeWarp.o \
		IdentifySilenceLocations.cc ../FileUtilities/SpeechFeatures.o \
		../FileUtilities/ImageIO.o

CreateBestSegments: DynamicTimeWarp.o CreateBestSegments.cc
	$(CC) $(CFLAGS) $(INCLUDE) -o CreateBestSegments DynamicTimeWarp.o \
		CreateBestSegments.cc ../FileUtilities/SpeechFeatures.o \
		../FileUtilities/ImageIO.o

IdentifyWordAlignments: DynamicTimeWarp.o IdentifyWordAlignments.cc
	$(CC) $(CFLAGS) $(INCLUDE) -o IdentifyWordAlignments DynamicTimeWarp.o \
		IdentifyWordAlignments.cc ../FileUtilities/SpeechFeatures.o \
		../FileUtilities/ImageIO.o

TestDTW: DynamicTimeWarp.o
	$(CC) $(CFLAGS) $(INCLUDE) -o testdtw DynamicTimeWarp.o testdtw.cc \
		../FileUtilities/SpeechFeatures.o ../FileUtilities/ImageIO.o

DynamicTimeWarp.o: DynamicTimeWarp.cc DynamicTimeWarp.h
	$(CC) $(CFLAGS) $(INCLUDE) -c DynamicTimeWarp.cc \
	  ../FileUtilities/ImageIO.o

clean:
	rm *.o testdtw IdentifyWordAlignments CreateBestSegments FindSilenceExemplars
