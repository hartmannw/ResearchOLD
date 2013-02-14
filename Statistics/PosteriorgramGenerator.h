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

  // Set the Gaussians used for calculating posteriors. If a vector of indices
  // are supplied, it is used to group the gaussians for the posterior 
  // calculation. It is assumed the vector is properly formed where every index
  // between 0 and N is present at least once. Also, N < mog_.size(). If no 
  // vector is specificied, each Gaussian is assumed to be its own unique label.
  bool SetGaussians(std::vector<MixtureOfDiagonalGaussians> mog);
  bool SetGaussians(std::vector<MixtureOfDiagonalGaussians> mog, 
      std::vector<int> indices);

  // Compute the similarity matrix based on the Cauchy-Schwarz Divergance
  // measure. (See the MixtureOfDiagonalGaussian code for more information about
  // the technique.) Each element (i,j) computes the divergence between mog_[i]
  // and mog_[j]. The matrix contains value between 0 and 1, where the most
  // similar are 1. Each diagonal element should have a value of 1.
  utilities::Matrix<double> ComputeSimilarityMatrix();

  // Computes the posteriorgram where the rows are posteriors and columns are 
  // frames. Assumes the columns in data are also frames. 
  utilities::Matrix<double> ComputePosteriorgram(
      const utilities::Matrix<double> &data) const;

  // Given the posteriorgram, returns the index of the highest scoring posterior
  // in each frame.
  std::vector<int> BestIndexPerFrame(const utilities::Matrix<double> &pgram);

  // Returns a sorted list of pairs. where first is the posterior index and 
  // second is the number of frames where that posterior had the highest score.
  // List is sorted in descending order based on second, so first element in the
  // array is the posterior index that scored highest in the most number of
  // frames.
  std::vector< std::pair<int, int> > BestIndexCount(
      const utilities::Matrix<double> &pgram);

  // Similar to BestIndexCount except that it uses the total mass in the 
  // posteriorgram to weight the best posterior instead of the number of frames
  // it was the highest index.
  std::vector<std::pair<int, double> > BestIndexMass(      
       const utilities::Matrix<double> &pgram);

  // Internal function used to sort a vector of pairs properly.
  static bool CompareIndexCountPair(const std::pair<int, int> &a, 
      const std::pair<int, int> &b){ return a.second > b.second; }
  
  // Internal function used to sort a vector of pairs properly.
  static bool CompareIndexMassPair(const std::pair<int, double> &a, 
      const std::pair<int, double> &b){ return a.second > b.second; }
};

} // end namespace statistics
#endif
