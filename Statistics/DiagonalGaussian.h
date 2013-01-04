#ifndef STATISTICS_DIAGONALGUASSIAN_
#define STATISTICS_DIAGONALGUASSIAN_

#include <vector>
#include <cmath>

namespace statistics
{

class DiagonalGaussian
{
 public:
  DiagonalGaussian(){}
  ~DiagonalGaussian(){}

  void Initialize(std::vector<double> mean, std::vector<double> variance);
  void CalculateConstant();
  void AddMean(std::vector<double> mean);
  void AddVariance(std::vector<double> variance);
  double Likelihood(std::vector<double> point);
  double LogLikelihood(std::vector<double> point);
  double KLDivergence(DiagonalGaussian g);          // I believe these functions
  double SymmetricKLDivergence(DiagonalGaussian g); // are correct, but they
                                                    // should be tested more.

  unsigned int dimension(){ return mean_.size();}
  double mean(unsigned int i){ return mean_[i];}
  double variance(unsigned int i){return variance_[i];}
  std::vector<double> mean(){return mean_;}
  std::vector<double> variance(){return variance_;}
  double determinant();

 private:
  std::vector<double> mean_;
  std::vector<double> variance_;
  double constant_;
};

}

#endif
