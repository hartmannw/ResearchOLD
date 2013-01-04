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

class HmmSet
{
 private:
  std::map<std::string, unsigned int> mixture_index_;
  std::vector<MixtureOfDiagonalGaussians> states_;
  std::vector<HiddenMarkovModel> hmms_;

  bool FindHtkModelHeader(std::ifstream &fin, std::string &line);
  unsigned int LoadSingleHtkState(std::ifstream &fin, std::string &line, 
      std::string name);
  bool LoadSingleHtkHmm(std::ifstream &fin, std::string &line);
  DiagonalGaussian LoadSingleHtkGaussian(std::ifstream &fin, std::string &line);

 public:
  HmmSet(){}
  ~HmmSet(){}

  void LoadHtkHmmSet(std::string filename);
  std::vector<MixtureOfDiagonalGaussians> states(){return states_;}
  std::vector<std::vector<std::string> > mixture_names();

};

}

#endif
