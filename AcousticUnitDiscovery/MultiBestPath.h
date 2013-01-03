// William Hartmann (hartmannw@gmail.com)
// This is free and unencumbered software released into the public domain.
// See the UNLICENSE file for more information.

#ifndef ACOUSTICUNITDISCOVERY_MULTIBESTPATH_H
#define ACOUSTICUNITDISCOVERY_MULTIBESTPATH_H

#include<vector>

namespace acousticunitdiscovery
{

std::vector<int> FindBestPath(const std::vector<std::vector<double> > &pgram, 
    const std::vector<std::vector<double> > &transition, int min_frames);

}

#endif
