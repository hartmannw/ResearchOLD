// William Hartmann (hartmannw@gmail.com)
// This is free and unencumbered software released into the public domain.
// See the UNLICENSE file for more information.

#ifndef STATISTICS_POSTERIORGRAMGENERATOR_H_
#define STATISTICS_POSTERIORGRAMGENERATOR_H_

#include<vector>
#include<utility>
#include<algorithm>
#include "MixtureOfDiagonalGaussians.h"
#include "Matrix.h"

// Stores a vector of mixtures of diagonal Gaussians. Each mixture is
// associated with an index so multiple mixtures can contribute to the same
// posterior. Once the GMMs have been set, the class accepts a Matrix of data
// (feature x frame). The data can be used to produce a posteriorgram, the best
// index per frame, or a ranked list of outputs for the entire utterance.

namespace statistics
{

class PosteriorgramGenerator
{
 private:
  std::vector<MixtureOfDiagonalGaussians> mog_;
  std::vector<int> posterior_index_;

 public:
  PosteriorgramGenerator() {}
  ~PosteriorgramGenerator() {}

  bool SetGaussians(std::vector<MixtureOfDiagonalGaussians> mog);
  bool SetGaussians(std::vector<MixtureOfDiagonalGaussians> mog, 
      std::vector<int> indices);

  utilities::Matrix<double> ComputeSimilarityMatrix();
  utilities::Matrix<double> ComputePosteriorgram(
      const utilities::Matrix<double> &data);

  std::vector<int> BestIndexPerFrame(const utilities::Matrix<double> &pgram);

  // Returns a sorted list of pairs. where first is the posterior index and 
  // second is the number of frames where that posterior had the highest score.
  // List is sorted in descending order based on second, so first element in the
  // array is the posterior index that scored highest in the most number of
  // frames.
  std::vector< std::pair<int, int> > BestIndexCount(
      const utilities::Matrix<double> &pgram);

  // Internal function used to sort a vector of pairs properly.
  static bool CompareIndexCountPair(const std::pair<int, int> &a, 
      const std::pair<int, int> &b){ return a.second > b.second; }
};

} // end namespace statistics
#endif
