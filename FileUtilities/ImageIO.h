#ifndef FILEUTILITIES_IMAGEIO_
#define FILEUTILITIES_IMAGEIO_

#include<vector>
#include<string>
#include<fstream>
#include<limits>

namespace fileutilities
{

bool WritePlainPGM(const std::vector< std::vector<double> > &image,
    const unsigned int &max, const std::string &filename);

bool WriteBinaryPGM(const std::vector< std::vector<double> > &image,
    const std::string &filename);

}

#endif
