// William Hartmann (hartmannw@gmail.com)
// This is free and unencumbered software released into the public domain.
// See the UNLICENSE file for more information.

// Creates a Similarity Matrix for the HMMs in a set of Hidden Markov Models.

#include<cmath>
#include<iostream>
#include<fstream>
#include<vector>
#include<string>

#include "Matrix.h"
#include "PosteriorgramGenerator.h"
#include "StringFunctions.h"
#include "HmmSet.h"

int main(int argc, char* argv[])
{
  if( argc < 2)
  {
    std::cout<<"Usage is <HMM File>"<<std::endl;
    exit(0);
  }

  std::string hmmfile = std::string(argv[1]);

  statistics::HmmSet htk;
  htk.LoadHtkHmmSet(hmmfile);
  std::vector<statistics::MixtureOfDiagonalGaussians> mog = htk.states();
  std::vector<statistics::HiddenMarkovModel> hmmset = htk.hmms();
  statistics::PosteriorgramGenerator pg;                                         

  pg.SetGaussians(mog);
  utilities::Matrix<double> sm = pg.ComputeSimilarityMatrix();

  utilities::Matrix<double> hmm_sm;
  hmm_sm.Initialize(hmmset.size(), hmmset.size(), 0);

  for(unsigned int i = 0; i < hmmset.size(); ++i)
    for(unsigned int j = i+1; j < hmmset.size(); ++j)
    {
      unsigned int states = 
          std::min(hmmset[i].NumberOfStates(), hmmset[j].NumberOfStates());
      for(unsigned int s = 0; s < states; ++s)
        hmm_sm(i,j) += (1 - sm(hmmset[i].state(s), hmmset[j].state(s)));
      hmm_sm(j,i) = hmm_sm(i,j);
    }

  for(unsigned int r = 0; r < hmm_sm.NumRows(); ++r)
  {
    for(unsigned int c = 0; c < hmm_sm.NumCols(); ++c)
      std::cout<<hmm_sm(r,c)<<" ";
    std::cout<<"\n";
  }

  return 0;
}
