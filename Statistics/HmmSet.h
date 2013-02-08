// William Hartmann (hartmannw@gmail.com)
// This is free and unencumbered software released into the public domain.
// See the UNLICENSE file for more information.

#ifndef STATISTICS_HMMSET_
#define STATISTICS_HMMSET_

#include<map>
#include<vector>
#include<string>
#include<fstream>
#include<iostream>

#include "DiagonalGaussian.h"
#include "MixtureOfDiagonalGaussians.h"
#include "HiddenMarkovModel.h"
#include "StringFunctions.h"

namespace statistics
{

// Stores a set of HMMs. Intended for use with a set of HMMs trained using the
// hidden markov model toolkit (HTK).
class HmmSet
{
 private:
  // Allows a MOG to be accessed by name instead of ID.
  std::map<std::string, unsigned int> mixture_index_;
  // Allows an HMM to be accessed by name instead of ID.
  std::map<std::string, unsigned int> hmm_index_;

  // Since it is common to share MOGs across different states in different HMMs,
  // we store them all together.
  std::vector<MixtureOfDiagonalGaussians> states_;

  // The HMM maintains information about which state uses which MOG.
  std::vector<HiddenMarkovModel> hmms_;

  // Set of functions used for reading in a set of HMMs from an HTK file.
  bool FindHtkModelHeader(std::ifstream &fin, std::string &line);
  unsigned int LoadSingleHtkState(std::ifstream &fin, std::string &line, 
      std::string name);
  bool LoadSingleHtkHmm(std::ifstream &fin, std::string &line);
  DiagonalGaussian LoadSingleHtkGaussian(std::ifstream &fin, std::string &line);

 public:
  HmmSet(){}
  ~HmmSet(){}

  // Loads the HMMs from a given HTK file. Assumes the HMMs uses mixures of
  // Gaussians with diagonal covariance matrices. Other types of models are
  // possible with HTK, but are not supported by this function.
  void LoadHtkHmmSet(std::string filename);

  // Standard accessor functions.
  std::vector<HiddenMarkovModel> hmms() { return hmms_; }
  std::vector<MixtureOfDiagonalGaussians> states(){return states_;}
  std::vector<std::vector<std::string> > mixture_names();
  HiddenMarkovModel Hmm(std::string name) { return hmms_[ hmm_index_[name] ];}

};

}

#endif
