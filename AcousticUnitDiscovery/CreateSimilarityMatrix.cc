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

double GiniIndex(std::vector<double> u, bool normalize)
{
  double ret = 0;
  double sumu = 0;
  double N = static_cast<double>(u.size());
  for(unsigned int i = 0; i < u.size(); ++i)
    sumu+= u[i];
  std::sort(u.begin(), u.end());
  for(unsigned int i = 0; i < u.size(); ++i)
  {
    double k = static_cast<double>(i);
    ret += (u[i] / sumu) * ( (N - (k+1) +0.5) / N);
  }
  ret = 1 - (2 * ret);
  if(normalize)
    ret = (N * ret) / (N-1);
  return ret;
}

std::vector<std::string> SplitTriphone(std::string s)
{
  std::vector<std::string> ret, tokens;
  utilities::TokenizeString(s, '-', tokens);
  if(tokens.size() > 1)
  {
    ret.push_back(tokens[0]);
    s = tokens[1];
  }
  else
    ret.push_back("");
  tokens.clear();

  utilities::TokenizeString(s, '+', tokens);
  ret.push_back(tokens[0]);
  if(tokens.size() > 1)
    ret.push_back(tokens[1]);
  else
    ret.push_back("");


  //std::cout<<ret[0]<<"-"<<ret[1]<<"+"<<ret[2]<<": "<<origs<<std::endl;

  return ret;

}

utilities::Matrix<double> GenerateLocationMatrix(
    statistics::HiddenMarkovModel hmm, unsigned int steps, 
    double &expected_length)
{
  utilities::Matrix<double> ret;
  unsigned int states = hmm.NumberOfStates();
  ret.Initialize(states, steps, 0);
  ret(0,0) = 1;
  for(unsigned int t = 1; t < steps; ++t)
  {
    ret(0,t) = ret(0,t-1) * hmm.transition(1,1);
    for(unsigned int s=1; s < states; ++s)
      ret(s,t) = (ret(s-1,t-1) * hmm.transition(s,s+1)) + 
        (ret(s,t-1) * hmm.transition(s+1,s+1));
  }

  expected_length = 0;
  for(unsigned int r = 0; r < states; ++r)
  {
    for(unsigned int c = 0; c < steps; ++c)
    {
      expected_length += ret(r,c);
    }
  }
  return ret;
}

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

  // Remove sp and sil from the hmm list.
  for (std::vector<statistics::HiddenMarkovModel>::iterator it = hmmset.begin();
      it != hmmset.end(); ++it)
  {
    if( it->name() == "sp" || it->name() == "sil")
    {
      hmmset.erase(it);
      it = hmmset.begin();
    }
    std::vector<std::string> tri = SplitTriphone(it->name());
    if( tri[0].length() < 1 || tri[2].length() < 1)
    {
      hmmset.erase(it);
      it = hmmset.begin();
    }
  }

  utilities::Matrix<double> hmm_sm, occupancy;
  hmm_sm.Initialize(hmmset.size(), hmmset.size(), 0);
  occupancy.Initialize(hmmset.size(), 3, 0);

  // Calculate the stationary distribution / state occupancy for each HMM.
  for(unsigned int i = 0; i < hmmset.size(); ++i)
  {
    double b,c;
    //a = hmmset[i].transition(3,4) / (1 - hmmset[i].transition(1,1));
    b = hmmset[i].transition(1,2) / (1 - hmmset[i].transition(2,2));
    c = hmmset[i].transition(2,3) / (1 - hmmset[i].transition(3,3));
    occupancy(i,0) = 1 / (1 + b + (c*b));
    occupancy(i,1) = b * occupancy(i,0);
    occupancy(i,2) = c * occupancy(i,1);
  }

  // Calculate the probability of being in a particular state at any given time
  std::vector<utilities::Matrix<double> > location;
  std::vector<double> expected_length;
  for(unsigned int i = 0; i < hmmset.size(); ++i)
  {
    double elength = 0;
    location.push_back(GenerateLocationMatrix(hmmset[i], 100, elength));
    expected_length.push_back(elength);
  }

  for(unsigned int i = 0; i < hmmset.size(); ++i)
    for(unsigned int j = i; j < hmmset.size(); ++j)
    {
      unsigned int statesi = hmmset[i].NumberOfStates();
      unsigned int statesj = hmmset[j].NumberOfStates();
      utilities::Matrix<double> correspondence;
      correspondence.Initialize(statesi, statesj, 0);
      for(unsigned int si = 0; si < statesi; ++si)
      {
        for(unsigned int sj = 0; sj < statesj; ++sj)
        {
          correspondence(si,sj) = occupancy(i, si) * occupancy(j,sj);// *  
              //sm(hmmset[i].state(si), hmmset[j].state(sj)); 
              //hmmset[i].transition(si+1, si+1) *
              //hmmset[j].transition(sj+1, sj+1);
          //hmm_sm(i,j) += correspondence(si,sj);
        }
      }

      utilities::Matrix<double> correspondence2;
      correspondence2.Initialize(statesi, statesj, 0);
      for(unsigned int t = 0; t < location[i].NumCols(); ++t)
      {
        for(unsigned int si = 0; si < statesi; ++si)
          for(unsigned int sj = 0; sj < statesj; ++sj)
          {
            correspondence2(si,sj) += location[i](si,t) * location[j](sj,t);
          }
      }
     
      for(unsigned int r = 0; r < correspondence2.NumRows(); ++r)
      {
        for(unsigned int c = 0; c < correspondence2.NumCols(); ++c)
        {
          correspondence2(r,c) = correspondence2(r,c) / 
            std::max(expected_length[i], expected_length[j]);
          hmm_sm(i,j) += (correspondence2(r,c) * 
              sm(hmmset[i].state(r), hmmset[j].state(c)) );
        }
      }

      // Normalize correspondence matrix
      for(unsigned int r = 0; r < correspondence.NumRows(); ++r)
        for(unsigned int c = 0; c < correspondence.NumCols(); ++c)
          correspondence(r,c) = correspondence(r,c) / hmm_sm(i,j);

      // Now compute similarity
      /*
      double lvalue = 0, rvalue = 0;
      double Mr = static_cast<double>(statesi);
      double Mc = static_cast<double>(statesj);
      for(unsigned int r = 0; r < correspondence.NumRows(); ++r)
        lvalue += GiniIndex(correspondence.GetRow(r), true);
      for(unsigned int c = 0; c < correspondence.NumCols(); ++c)
        rvalue += GiniIndex(correspondence.GetCol(c), true);
      lvalue = lvalue * (1 / Mr);
      rvalue = rvalue * (1 / Mc);
      //hmm_sm(i,j) = 0.5 * (lvalue + rvalue);
      */
      std::vector<std::string> namei;
      std::vector<std::string> namej;
      namei = SplitTriphone(hmmset[i].name());
      namej = SplitTriphone(hmmset[j].name());

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
