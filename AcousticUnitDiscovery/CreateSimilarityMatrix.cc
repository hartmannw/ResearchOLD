// William Hartmann (hartmannw@gmail.com)
// This is free and unencumbered software released into the public domain.
// See the UNLICENSE file for more information.

// Creates a Similarity Matrix for the GMMs in Hidden Markov Model.

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
  if( argc < 3)
  {
    std::cout<<"Usage is <HMM File> <alpha update>"<<std::endl;
    exit(0);
  }

  std::string hmmfile = std::string(argv[1]);
  double alpha_update = utilities::ToNumber<double>( std::string(argv[2]) );

  statistics::HmmSet htk;
  htk.LoadHtkHmmSet(hmmfile);
  std::vector<statistics::MixtureOfDiagonalGaussians> mog = htk.states();
  std::vector<statistics::HiddenMarkovModel> hmmset = htk.hmms();
  statistics::PosteriorgramGenerator pg;                                         

  pg.SetGaussians(mog);
  utilities::Matrix<double> sm = pg.ComputeSimilarityMatrix();

  for(unsigned int h = 0; h < hmmset.size(); ++h)
  {
    for(unsigned int i = 0; i < hmmset[h].NumberOfStates(); ++i)
      for(unsigned int j = (i+1); j < hmmset[h].NumberOfStates(); ++j)
      {
        unsigned int state_i = hmmset[h].state(i);
        unsigned int state_j = hmmset[h].state(j);
        sm(state_i, state_j) = ( alpha_update * sm(state_i, state_j) ) +
          ( 1 - alpha_update );
        sm(state_j, state_i) = sm(state_i, state_j);
      }
  }

  for(unsigned int r = 0; r < sm.NumRows(); ++r)
  {
    for(unsigned int c = 0; c < sm.NumCols(); ++c)
      std::cout<<sm(r,c)<<" ";
    std::cout<<"\n";
  }

  return 0;
}
