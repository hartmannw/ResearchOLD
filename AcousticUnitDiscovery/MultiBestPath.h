// William Hartmann (hartmannw@gmail.com)
// This is free and unencumbered software released into the public domain.
// See the UNLICENSE file for more information.

#ifndef ACOUSTICUNITDISCOVERY_MULTIBESTPATH_H
#define ACOUSTICUNITDISCOVERY_MULTIBESTPATH_H

#include<vector>
#include<cmath>

namespace acousticunitdiscovery
{

// ViterbiInfo stores the information needed for the viterbi algorithm so that
// the best path can be determined and backtracked through.
typedef struct
{
  int parent;   // Index to the parent state.
  double score; // Score for this particular point along the path.
} ViterbiInfo;

std::vector<std::vector<double> > GenerateTransitionMatrix(
    unsigned int states, double self_loop_prob);

std::vector<int> FindBestPath(const std::vector<std::vector<double> > &pgram, 
    const std::vector<std::vector<double> > &transition, int min_frames);

}

#endif
