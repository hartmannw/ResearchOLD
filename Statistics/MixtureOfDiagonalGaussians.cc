// William Hartmann (hartmannw@gmail.com)
// This is free and unencumbered software released into the public domain.
// See the UNLICENSE file for more information.

#include "MixtureOfDiagonalGaussians.h"

namespace statistics
{

void MixtureOfDiagonalGaussians::SetAllGaussians(
    std::vector<DiagonalGaussian> gaussian,std::vector<double> weight)
{
  gaussian_ = gaussian;
  weight_ = weight;
}

double MixtureOfDiagonalGaussians::Likelihood(std::vector<double> &point)
{
  double ret = 0;
  for(unsigned int i = 0; i < gaussian_.size(); ++i)
    ret+= (weight_[i] * gaussian_[i].Likelihood(point));
  return ret;
}

double MixtureOfDiagonalGaussians::LogLikelihood(std::vector<double> &point)
{
  return log(Likelihood(point));
}

// Implementation matches the results on the testset provided in the orginal
// author's paper, however, this implementation does not exactly match my
// reading of the paper. The author does have a Matlab implementation and this
// implementation matches that. I tried to contact the author to resolve the
// discrepency, but was unsuccessful. Since this matches the author's
// implementation and evaluates their testset correctly, I am willing to trust
// this implementation is the correct one.
double MixtureOfDiagonalGaussians::CSDivergence(MixtureOfDiagonalGaussians mog)
{
  // Divergence measure is split into three terms.
  double first_term = 0, second_term = 0, third_term = 0;
  //First and second terms
  double inner_term = 0;
  for(unsigned int i = 0; i < components(); ++i)
  {
    for(unsigned int j = 0; j < mog.components(); ++j)
    {
      DiagonalGaussian g = mog.gaussian(j);
      g.AddVariance(gaussian_[i].variance());
      first_term += (weight_[i] * mog.weight(j) * 
          g.Likelihood(gaussian_[i].mean()) );
    }
    for(unsigned int k = 0; k < components(); ++k)
    {
      DiagonalGaussian g = gaussian_[k];
      g.AddVariance(gaussian_[i].variance());
      inner_term += (weight_[i] * weight_[k] *
          g.Likelihood(gaussian_[i].mean()));
    }
  }
  second_term += (1 * inner_term);

  //Third term
  inner_term = 0;
  for(unsigned int j = 0; j < mog.components(); ++j)
  {
    for(unsigned int k = 0; k < mog.components(); ++k)
    {
      DiagonalGaussian g = mog.gaussian(j);
      g.AddVariance(mog.gaussian(k).variance());
      inner_term += (mog.weight(j) * mog.weight(k) *
          g.Likelihood(mog.gaussian(k).mean()));
    }
  }
  third_term += (1 * inner_term);

  first_term =  -std::log(first_term);
  second_term = 0.5 * std::log(second_term);
  third_term = 0.5 * std::log(third_term);
  return (first_term + second_term + third_term);
}

std::vector<double> MixtureOfDiagonalGaussians::WeightedMean()
{
  unsigned int dimension = gaussian_[0].dimension();
  std::vector<double> ret(dimension, 0);
  /*
  double max_weight = 0;
  for(unsigned int c = 0; c < gaussian_.size(); ++c)
  {
    if( weight_[c] > max_weight)
    {
      max_weight = weight_[c];
      for(unsigned int m = 0; m < dimension; ++m)
        ret[m] = gaussian_[c].mean(m);
    }
  }
  */ 
  for(unsigned int c = 0; c < gaussian_.size(); ++c)
    for(unsigned int m = 0; m < dimension; ++m)
      ret[m] += (weight_[c] * gaussian_[c].mean(m));
  return ret;
}

void MixtureOfDiagonalGaussians::NormalizeWeights()
{
  double total = 0;
  for(unsigned int i = 0; i < weight_.size(); ++i)
    total += weight_[i];
  for(unsigned int i = 0; i < weight_.size(); ++i)
    weight_[i] = weight_[i] / total;
}

std::vector<double> MixtureOfDiagonalGaussians::sample(
    std::default_random_engine &generator)
{
  unsigned int dimension = gaussian_[0].dimension();
  std::vector<double> ret(dimension, 0);
  for(unsigned int i = 0; i < weight_.size(); ++i)
  {
    std::vector<double> data = gaussian_[i].sample(generator);
    for(unsigned int m = 0; m < dimension; ++m)
      ret[m] += (weight_[i] * data[m]);
  }
  return ret;
}

}
