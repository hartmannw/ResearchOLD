// William Hartmann (hartmannw@gmail.com)
// This is free and unencumbered software released into the public domain.
// See the UNLICENSE file for more information.

#include "PosteriorgramGenerator.h"

namespace statistics
{

bool PosteriorgramGenerator::SetGaussians(
    std::vector<MixtureOfDiagonalGaussians> mog)
{
  mog_ = mog;
  posterior_index_.resize(mog.size());
  // Since no index is specified, we assume each MOG is its own distinct label.
  for(unsigned int i = 0; i < mog_.size(); ++i)
    posterior_index_[i] = i;
  return true;
}

bool PosteriorgramGenerator::SetGaussians(
    std::vector<MixtureOfDiagonalGaussians> mog, std::vector<int> indices)
{
  mog_ = mog;
  posterior_index_ = indices;
  return true;
}

// The similarity matrix is always symmetric so we only compute one half and 
// fill in the rest of the values based on the result.
utilities::Matrix<double> PosteriorgramGenerator::ComputeSimilarityMatrix()
{
  utilities::Matrix<double> ret;
  double max_divergence = 0;
  ret.Initialize(mog_.size(), mog_.size());
  for(unsigned int r = 0; r < mog_.size(); ++r)
    for(unsigned int c = r; c < mog_.size(); ++c)
    {
      ret(r,c) = mog_[r].CSDivergence(mog_[c]);
      if( ret(r,c) > max_divergence)
        max_divergence = ret(r,c);
    }
  // Normalize the matrix so that the values are between 0 and 1.
  for(unsigned int r = 0; r < ret.NumRows(); ++r)
    for(unsigned int c = r; c < ret.NumCols(); ++c)
    {
      ret(r,c) = 1 - ( ret(r,c) / max_divergence );
      ret(c,r) = ret(r,c);
    }

  return ret;
}

// Assume data is in a (feature x frame) format.
utilities::Matrix<double> PosteriorgramGenerator::ComputePosteriorgram(
    const utilities::Matrix<double> &data)
{
  utilities::Matrix<double> ret;

  // First identify the total number of posteriors.
  int posteriors = 0;
  int frames = data.NumCols();
  for(unsigned int i = 0; i < posterior_index_.size(); ++i)
    if( posterior_index_[i] > posteriors)
      posteriors = posterior_index_[i]; 
  posteriors++; // Number of posteriors is the highest index + 1.
  
  ret.Initialize(posteriors, frames, 0);
  for(int f = 0; f < frames; ++f)
  {
    std::vector<double> frame_data = data.GetCol(f);
    double total_value = 0;
    for(unsigned int g = 0; g < mog_.size(); ++g)
    {
      double like =  mog_[g].Likelihood(frame_data);
      ret(posterior_index_[g], f) += like;
      total_value +=  like;
    }
    // Normalize the frame by total mass in frame.
    for(int p = 0; p < posteriors; ++p)
      ret(p, f) = ( ret(p, f) / total_value );

  }
  return ret;
}

std::vector<int> PosteriorgramGenerator::BestIndexPerFrame(
    const utilities::Matrix<double> &pgram)
{
  std::vector<int> ret;
  ret.resize( pgram.NumCols() );

  for(unsigned int f = 0; f < ret.size(); ++f)
  {
    double bestval = -1;
    for(unsigned p = 0; p < pgram.NumRows(); ++p)
    {
      if( pgram(p, f) > bestval )
      {
        bestval = pgram(p, f);
        ret[f] = p;
      }
    }
  }
  return ret;
}

std::vector< std::pair<int, int> > PosteriorgramGenerator::BestIndexCount(                             
      const utilities::Matrix<double> &pgram)
{
  std::vector< std::pair<int, int> > ret;
  std::vector<int> one_best = BestIndexPerFrame(pgram);
  
  ret.resize(pgram.NumRows());
  // Initialize the vector of pairs
  for(unsigned int i = 0; i < ret.size(); ++i)
    ret[i] = std::make_pair(i, 0);
  for(unsigned int i = 0; i < one_best.size(); ++i)
    ret[ one_best[i] ].second++;

  std::sort(ret.begin(), ret.end(), CompareIndexCountPair);
  return ret;
}

std::vector<std::pair<int, double> > PosteriorgramGenerator::BestIndexMass(
    const utilities::Matrix<double> &pgram)
{
  std::vector<std::pair<int, double> > ret;
  ret.resize(pgram.NumRows());
  double total_mass = 0;
  // Initialize the vector of pairs.
  for(unsigned int i = 0; i < ret.size(); ++i)
    ret[i] = std::make_pair(i,0);

  for(unsigned int f = 0; f < pgram.NumCols(); ++f)
    for(unsigned int p = 0; p < pgram.NumRows(); ++p)
    {
      ret[p].second += pgram(p,f);
      total_mass += pgram(p,f);
    }
 
  // Normalize by total mass to produce posterior.
  for(unsigned int i = 0; i < ret.size(); ++i)
    ret[i].second = ret[i].second / total_mass;

  std::sort(ret.begin(), ret.end(), CompareIndexMassPair);
  return ret;
}

} // end namespace statistics
