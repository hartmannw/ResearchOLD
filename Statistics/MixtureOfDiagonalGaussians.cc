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

double MixtureOfDiagonalGaussians::CSDivergence(MixtureOfDiagonalGaussians mog)
{
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
    //std::cout<<gaussian_[i].dimension()<<" "<<gaussian_[i].dimension() / 2.0<<std::endl;
    /*second_term += ( std::pow(weight_[i],2)  * 
        std::sqrt(1 / gaussian_[i].determinant()) ) /
        std::pow(2*pi, gaussian_[i].dimension() / 2.0);*/
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
    /*third_term += ( std::pow(mog.weight(j),2) * 
        std::sqrt(1 / mog.gaussian(j).determinant()) ) /
        std::pow(2*pi, mog.gaussian(j).dimension() / 2.0);*/
    for(unsigned int k = 0; k < mog.components(); ++k)
    {
      DiagonalGaussian g = mog.gaussian(j);
      g.AddVariance(mog.gaussian(k).variance());
      //std::cout<<g.mean(0)<<" "<<g.mean(1)<<" "<<g.variance(0)<<" "<<g.variance(1)<<std::endl;
      inner_term += (mog.weight(j) * mog.weight(k) *
          g.Likelihood(mog.gaussian(k).mean()));
    }
  }
  third_term += (1 * inner_term);

  first_term =  -std::log(first_term);
  second_term = 0.5 * std::log(second_term);
  third_term = 0.5 * std::log(third_term);
//  return second_term;
  return (first_term + second_term + third_term);
}

}
