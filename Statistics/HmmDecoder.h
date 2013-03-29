// William Hartmann (hartmannw@gmail.com)
// This is free and unencumbered software released into the public domain.
// See the UNLICENSE file for more information.

#ifndef STATISTICS_HMMDECODER_H_
#define STATISTICS_HMMDECODER_H_

#include <vector>
#include <algorithm>

#include "HmmSet.h"
#include "Matrix.h"
#include "HiddenMarkovModel.h"

namespace statistics
{

// ViterbiInfo stores the information needed for the viterbi algorithm so that
// the best path can be determined and backtracked through.
typedef struct
{
  int parent;   // Index to the parent state.
  double score; // Score for this particular point along the path.
} ViterbiInfo;

class HmmDecoder
{
 private:
  HmmSet *htk_;
  utilities::Matrix<double> *lm_;
  double min_log_;
  double zero_log_;
  std::vector<utilities::Matrix<ViterbiInfo> > forward_prob;
  std::vector<utilities::Matrix<ViterbiInfo> > back_prob;

  ViterbiInfo SetMaximumScore(double a, double b, int a_index, int b_index);
  double SafeLog(double value){ return std::max(std::log(value), zero_log_); }

 public:
  HmmDecoder() : min_log_(-1000000), zero_log_(-30){}
  ~HmmDecoder(){}

  void Initialize(HmmSet *htk, utilities::Matrix<double> *lm);
  bool FeaturesStateScores(const utilities::Matrix<double> &data, 
      utilities::Matrix<double> &scores);

  bool ForwardPass(const utilities::Matrix<double> &scores);
  std::vector<int> BestForwardPath();
  std::vector<int> ViterbiPath(const utilities::Matrix<double> &scores);

  double ForceAlignScore( const utilities::Matrix<double> &scores,
      const std::vector<int> &sequence );
  std::vector<double> ForceAlignExits( const utilities::Matrix<double> &scores,
      const std::vector<int> &sequence );

};

}

#endif
