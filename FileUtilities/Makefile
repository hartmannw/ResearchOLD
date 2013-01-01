CC=g++
CFLAGS=-Wall

all: Test ConvertCep2Ascii ConvertCep2Htk ImageIO.o

Test: SpeechFeatures.o ImageIO.o
	$(CC) $(CFLAGS) -o test SpeechFeatures.o ImageIO.o test.cc

ConvertCep2Ascii:  SpeechFeatures.o
	$(CC) $(CFLAGS) -o ConvertCep2Ascii SpeechFeatures.o ConvertCep2Ascii.cc

ConvertCep2Htk: SpeechFeatures.o
	$(CC) $(CFLAGS) -o ConvertCep2Htk SpeechFeatures.o ConvertCep2Htk.cc

SpeechFeatures.o: SpeechFeatures.cc SpeechFeatures.h
	$(CC) $(CFLAGS) -c SpeechFeatures.cc

ImageIO.o: ImageIO.cc ImageIO.h
	$(CC) $(CFLAGS) -c ImageIO.cc

clean:
	rm *.o test ConvertCep2Htk ConvertCep2Ascii
