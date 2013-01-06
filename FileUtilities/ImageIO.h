// William Hartmann (hartmannw@gmail.com)
// This is free and unencumbered software released into the public domain.
// See the UNLICENSE file for more information.

// Contains basic functions for dealing with some images. Currently it only
// supports saving a matrix of double data as a binary PGM.

#ifndef FILEUTILITIES_IMAGEIO_
#define FILEUTILITIES_IMAGEIO_

#include<vector>
#include<string>
#include<fstream>
#include<limits>

namespace fileutilities
{

// Normalizes the data to be between 0 and max. The data is then written to an 
// ASCII file in the PGM format.
bool WritePlainPGM(const std::vector< std::vector<double> > &image,
    const unsigned int &max, const std::string &filename);

// Normalizes the data to be between 0 and 255. Data is then written to a binary
// PGM file.
bool WriteBinaryPGM(const std::vector< std::vector<double> > &image,
    const std::string &filename);

}

#endif
