CC=g++
CFLAGS=-Wall -O2
INCLUDE=-I../Utilities/

all: Test ModelSolarSystem

ModelSolarSystem: OrbitalBody.o GravitationalSystem.o ModelSolarSystem.cc
	$(CC) $(CFLAGS) $(INCLUDE) -o ModelSolarSystem OrbitalBody.o \
		GravitationalSystem.o ../Utilities/StringFunctions.o ModelSolarSystem.cc

Test: OrbitalBody.o GravitationalSystem.o test.cc
	$(CC) $(CFLAGS) $(INCLUDE) -o test OrbitalBody.o GravitationalSystem.o test.cc

GravitationalSystem.o: GravitationalSystem.h GravitationalSystem.cc
	$(CC) $(CFLAGS) $(INCLUDE) -c GravitationalSystem.cc

OrbitalBody.o: OrbitalBody.cc OrbitalBody.h
	$(CC) $(CFLAGS) $(INCLUDE) -c OrbitalBody.cc

clean:
	rm *.o test ModelSolarSystem
