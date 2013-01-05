#include "HmmSet.h"

namespace statistics
{

void HmmSet::LoadHtkHmmSet(std::string filename)
{
  std::ifstream fin(filename.c_str(), std::ios::in);
  std::string line;
  std::cout<<filename<<std::endl;

  while(FindHtkModelHeader(fin, line))
  {
    std::cout<<line<<std::endl;
    if(line.substr(0,4) == std::string("~s \""))
    {
      std::string name = line.substr(4, line.length()-5);
      getline(fin, line);
      LoadSingleHtkState(fin, line, name);
    }
    else
    {
      LoadSingleHtkHmm(fin, line);
    }
  }

  fin.close();
}

bool HmmSet::FindHtkModelHeader(std::ifstream &fin, std::string &line)
{
  while( fin.good() )
  {
    if(line.length() > 4 && (line.substr(0,4) == std::string("~s \"") || 
        line.substr(0,4) == std::string("~h \"")))
    {
      return true;
    }
    std::getline(fin, line);
  }
  return false;
}

unsigned int HmmSet::LoadSingleHtkState(std::ifstream &fin, std::string &line, 
    std::string name)
{
  MixtureOfDiagonalGaussians mog;
  // Reference to previously loaded state.
  if(line.length() > 4 && line.substr(0,4) == std::string("~s \""))
  {
    name = line.substr(4, line.length()-5);
    getline(fin, line);
    std::cout<<mixture_index_[name]<<std::endl;
    return mixture_index_[name];
  }
  else
  {
    if(line.find( std::string("<NUMMIXES>") ) == 0)
    {
      std::vector<std::string> tokens;
      line = utilities::TrimString(line);
      utilities::TokenizeString(line, ' ', tokens);
      unsigned int num_mixes = utilities::ToNumber<unsigned int>(tokens[1]);
      getline(fin, line);
      for(unsigned int i = 0; i < num_mixes; ++i)
      {
        tokens.clear();
        line = utilities::TrimString(line);
        utilities::TokenizeString(line, ' ', tokens);
        double mixture_weight = utilities::ToNumber<float>(tokens[2]);
        getline(fin, line);
        mog.AddGaussian(LoadSingleHtkGaussian(fin, line), mixture_weight);
      }
    }
    else
    {
      DiagonalGaussian g;
      g = LoadSingleHtkGaussian(fin, line);
      mog.AddGaussian(g, 1.0);
    }
  }
  unsigned int next_index = states_.size();
  states_.push_back(mog);
  mixture_index_[name] = next_index;
  return next_index;
}

bool HmmSet::LoadSingleHtkHmm(std::ifstream &fin, std::string &line)
{
  // Assume we begin by looking at the ~h line
  std::string name = line.substr(4, line.length()-5);
  std::cout<<name<<std::endl;
  HiddenMarkovModel hmm;
  hmm.SetName(name);
  getline(fin, line); // BEGINHMM
  getline(fin, line); // NUMSTATES
  std::vector<std::string> tokens;
  line = utilities::TrimString(line);
  utilities::TokenizeString(line, ' ', tokens);
  unsigned int state_count = utilities::ToNumber<unsigned int>(tokens[1]) - 2;
  getline(fin, line); // STATE
  for(unsigned int i = 0; i < state_count; ++i)
  {
    if( line.find("<STATE>") > line.length())
      return false;
    std::string state_name = name + std::string("_") + utilities::ToString(i);
    std::cout<<state_name<<std::endl;
    getline(fin, line);
    hmm.AddState( LoadSingleHtkState(fin, line, state_name));
  }
  while( line.find( std::string("<ENDHMM>")) != 0)
    getline(fin, line);

  return true;

}

DiagonalGaussian HmmSet::LoadSingleHtkGaussian(std::ifstream &fin, 
    std::string &line)
{
  DiagonalGaussian g;
  // Assuming first line is mean
  getline(fin, line);
  line = utilities::TrimString(line);
  std::vector<std::string> tokens;
  utilities::TokenizeString(line, ' ', tokens);
  std::vector<double> mean, variance;
  for(unsigned int i =0; i < tokens.size(); ++i)
    mean.push_back(utilities::ToNumber<double>(tokens[i]));
  getline(fin, line);
  getline(fin, line);
  line = utilities::TrimString(line);
  tokens.clear();
  utilities::TokenizeString(line, ' ', tokens);
  for(unsigned int i = 0; i < tokens.size(); ++i)
    variance.push_back(utilities::ToNumber<double>(tokens[i]));
  g.Initialize(mean, variance);
  getline(fin, line); // Skip GCONST line
  getline(fin, line);

  return g;
}

std::vector<std::vector<std::string> > HmmSet::mixture_names()
{
  std::vector<std::vector<std::string> > ret;
  ret.resize(states_.size());
  for(std::map<std::string, unsigned int>::iterator it = mixture_index_.begin();
      it != mixture_index_.end(); ++it)
    ret[it->second].push_back(it->first);

  return ret;
}

}
