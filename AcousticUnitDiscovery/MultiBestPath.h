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

// MultiBestPath contains a set of functions for finding certain types of best
// paths. In general, the functions expect something like a posteriorgram. Where
// the rows are the state posteriors and the columns are the frames. A 
// transition matrix must also be supplied containing the probabilities of
// transitioning between any two states. A minimum number of frames can also be
// specified. Each state in the final path must exist for at least as many
// consecutive frames as specified by min_frames. Note that the algorithms work
// in the log domain and expect all probabilities to be log probabilities.
//
// Types of Best Paths:
//
// FindViterbiPath (standard viterbi): Finds the best sequence of states as in 
// any other implementation of the viterbi algorithm.
//
// FindRestrictedViterbiPath: Same as the previous path, except an initial
// sequence of states can be specified. Finds the best path that begins with the
// initial sequence of states. Can also perform a force alignment such that the
// resulting path is identical to the given path.
//
// ApproximateViterbiSet: Given a set of posteriorgrams, finds the path that 
// maximizes the likelihood for the entire set. The implementation is an 
// approximation and does not guarentee the true path.
//
// BestPathInSet: First finds the best path for each posteriorgram in the set.
// Then it selects from that set of paths the path that maximizes the 
// likelihood for the entire set.

namespace acousticunitdiscovery
{

// ViterbiInfo stores the information needed for the viterbi algorithm so that
// the best path can be determined and backtracked through.
typedef struct
{
  int parent;   // Index to the parent state.
  double score; // Score for this particular point along the path.
} ViterbiInfo;

// Generates a transition matrix where the diagonal elements are self_loop_prob
// and the off diagonal elements are (1-self_loop_prob). Assumes the 
// probabilities are not in the log domain.
utilities::Matrix<double> GenerateTransitionMatrix(
    unsigned int states, double self_loop_prob);

// Initial depracted version of the best path algorithm should be DELETED.
std::vector<int> FindBestPath(const utilities::Matrix<double> &pgram, 
    const utilities::Matrix<double> &transition, int min_frames);

// Standard viterbi best path algorithm operating on a posteriorgram.
std::vector<int> FindViterbiPath(const utilities::Matrix<double> &pgram,         
    const utilities::Matrix<double> &transition, int min_frames, 
    double &final_score);

// Returns the best path that starts with the best path in initial_path. If 
// force_align is set to true, then the final path is equal to initial_path.
std::vector<int> FindRestrictedViterbiPath(                                      
    const utilities::Matrix<double> &pgram,                                      
    const utilities::Matrix<double> &transition, int min_frames,                 
    std::vector<int> initial_path, bool force_align, double &final_score);

// Returns the one best path for a particular posteriorgram in the set that also
// maximizes the likelihood for the entire set.
std::vector<int> BestPathInSet(
    const std::vector<utilities::Matrix<double> > &pgram_set,
    const utilities::Matrix<double> &transition, int min_frames);

// Returns the best single path for an entire set of posteriorgrams. The 
// implementation is approximate, so the best path is not guaranteed.
std::vector<int> ApproximateViterbiSet(
    const std::vector<utilities::Matrix<double> > &pgram_set,
    const utilities::Matrix<double> &transition, int min_frames);

// Handles the logic of determining the state likelihood. Should only be used 
// by functions internal to MultiBestPath.
double GetStateScore(const utilities::Matrix<double> &pgram,                     
    const std::vector<int> &initial_path, int min_frames, int state, int frame);

// Handles the logic of determining the transition likelihood. Should only be 
// used by functions internal to MultiBestPath.
double GetTransitionScore(const utilities::Matrix<double> &transition,           
    const std::vector<int> &initial_path, int min_frames, int from_state,     
    int to_state);

// Used by ViterbiSet. Returns a sub path in the dp_matrix.
std::vector<int> BestSubPathInViterbiSet(                                        
    const utilities::Matrix<ViterbiInfo> &dp_matrix, int state, int frame);

// Once the dynamic programming matrix has been filled, this function is used to
// find the single best path in the matrix. Should only be used by internal
// functions.
std::vector<int> BestPathInDpMatrix(
    const utilities::Matrix<ViterbiInfo> &dp_matrix, unsigned int min_frames,
    const std::vector<int> &initial_path, bool force_align, 
    double &final_score);

}

#endif
