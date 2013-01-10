// William Hartmann (hartmannw@gmail.com)
// This is free and unencumbered software released into the public domain.
// See the UNLICENSE file for more information.

#ifndef ACOUSTICUNITDISCOVERY_MULTIBESTPATH_H
#define ACOUSTICUNITDISCOVERY_MULTIBESTPATH_H

#include<vector>
#include<cmath>
#include<algorithm>
#include "Matrix.h"
#include "ImageIO.h"

namespace acousticunitdiscovery
{

// ViterbiInfo stores the information needed for the viterbi algorithm so that
// the best path can be determined and backtracked through.
typedef struct
{
  int parent;   // Index to the parent state.
  double score; // Score for this particular point along the path.
} ViterbiInfo;

utilities::Matrix<double> GenerateTransitionMatrix(
    unsigned int states, double self_loop_prob);

std::vector<int> FindBestPath(const utilities::Matrix<double> &pgram, 
    const utilities::Matrix<double> &transition, int min_frames);

std::vector<int> FindViterbiPath(const utilities::Matrix<double> &pgram,         
    const utilities::Matrix<double> &transition, int min_frames, 
    double &final_score);

std::vector<int> FindRestrictedViterbiPath(                                      
    const utilities::Matrix<double> &pgram,                                      
    const utilities::Matrix<double> &transition, int min_frames,                 
    std::vector<int> initial_path, bool force_align, double &final_score);

std::vector<int> ApproximateViterbiSet(
    const std::vector<utilities::Matrix<double> > &pgram_set,
    const utilities::Matrix<double> &transition, int min_frames);

double GetStateScore(const utilities::Matrix<double> &pgram,                     
    const std::vector<int> &initial_path, int min_frames, int state, int frame);

double GetTransitionScore(const utilities::Matrix<double> &transition,           
    const std::vector<int> &initial_path, int min_frames, int from_state,     
    int to_state);

std::vector<int> BestSubPathInViterbiSet(                                        
    const utilities::Matrix<ViterbiInfo> &dp_matrix, int state, int frame);

std::vector<int> BestPathInDpMatrix(
    const utilities::Matrix<ViterbiInfo> &dp_matrix, unsigned int min_frames,
    const std::vector<int> &initial_path, bool force_align, 
    double &final_score);

}

#endif
