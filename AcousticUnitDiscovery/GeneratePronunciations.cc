// William Hartmann (hartmannw@gmail.com)
// This is free and unencumbered software released into the public domain.
// See the UNLICENSE file for more information.

// Used for generating pronunciations for a specific word. At the moment it 
// works from feature files. In the future I would like to change it so that it
// works from Posteriorgrams and they are generatedin a separate, earlier step.

#include<cmath>
#include<iostream>
#include<fstream>
#include<vector>
#include<string>

#include "Matrix.h"
#include "SpeechFeatures.h"
#include "PosteriorgramGenerator.h"
#include "StringFunctions.h"
#include "HmmSet.h"
#include "MultiBestPath.h"

std::vector<int> ReadVector(std::string filename)
{
  std::vector<int> ret;
  std::ifstream fin;
  fin.open(filename.c_str(), std::ios::in);
  while(fin.good())
  {
    double n;
    fin >> n;
    ret.push_back(static_cast<int>(n)-1);
  }
  ret.resize( ret.size() - 1 );
  return ret;
}

int GroupStatesByIndex(const std::vector<std::vector<std::string> > &original_names,
  std::vector<int> &index, std::vector<std::string> &names, std::string idxfile)
{
  int total = 0;
  index = ReadVector(idxfile);
  for(unsigned int i = 0; i < index.size(); ++i)
    if( index[i] > total)
      total = index[i];
  total++;
  names.resize(total);

  for(unsigned int i = 0; i < index.size(); ++i)
    names[index[i]] += original_names[i][0];
  return total;
}

int GroupStatesByName(const std::vector<std::vector<std::string> > &original_names,
  std::vector<int> &index, std::vector<std::string> &names)
{
  int total = 0;
  index.resize(original_names.size());
  names.resize(0);
  for(unsigned int i = 0; i < index.size(); ++i)
  {
    std::vector<std::string> tokens;
    utilities::TokenizeString(original_names[i][0], '_', tokens);
    int token_id = -1;
    for(unsigned int j = 0; j < names.size(); ++j)
      if(names[j] == tokens[0])
        token_id = j;

    if(token_id < 0)
    {
      index[i] = total;
      names.push_back(tokens[0]);
      ++total;
    }
    else
    {
      index[i] = token_id;
    }
  }
  return total;
}

std::vector<statistics::MixtureOfDiagonalGaussians> GroupMixturesByIndex( 
    const std::vector<statistics::MixtureOfDiagonalGaussians> &gmm,
    const std::vector<int> &index)
{
  std::vector<statistics::MixtureOfDiagonalGaussians> ret;

  // Find the total number of groups. 
  int total_groups = 0;
  for(unsigned int i = 0; i < index.size(); ++i)
    if(index[i] > total_groups)
      total_groups = index[i];
  total_groups++;
  ret.resize(total_groups);

  for(unsigned int i = 0; i < index.size(); ++i)
    for(unsigned int g = 0; g < gmm[i].components(); ++g)
      ret[ index[i] ].AddGaussian(gmm[i].gaussian(g), gmm[i].weight(g));

  for(unsigned int i = 0; i < ret.size(); ++i)
    ret[i].NormalizeWeights();

  return ret;
}

int main()
{
  std::string listfile("/people/hartmann/research/SegmentalModel/tidigit_exp/word_locations/eight.train.list");
  std::string suffix(".htk");
  std::string maindir("/people/hartmann/research/SegmentalModel/tidigit_exp/");
  std::string hmmfile("/people/hartmann/research/SegmentalModel/tidigit_exp/hmm_plp_grapheme11/hmm29/hmmdefs");
  int maximum_examples = 20;
  int min_frames = 3;
  int examples = 0;
  std::ifstream fin;

  fin.open(listfile.c_str());

  if( !fin.good())
  {
    std::cout<<"File "<<listfile<<" could not be opened.\n";
    exit(1);
  }

  statistics::HmmSet htk;
  htk.LoadHtkHmmSet(hmmfile);
  std::vector<statistics::MixtureOfDiagonalGaussians> mog = htk.states();        
  std::vector<std::vector<std::string> > names = htk.mixture_names();            
  std::vector<utilities::Matrix<double> > pgram_set;                            
  statistics::PosteriorgramGenerator pg;                                         
  utilities::Matrix<double> transition;
  std::vector<int> group_index;
  std::vector<std::string> group_name;

  int total_groups = GroupStatesByName(names, group_index, group_name);
  mog = GroupMixturesByIndex(mog, group_index);
  //int total_groups = GroupStatesByIndex(names, group_index, group_name, 
  //    std::string("idx_sm.pgm"));

  transition = 
      acousticunitdiscovery::GenerateTransitionMatrix(total_groups, 0.9);
  //pg.SetGaussians(mog, group_index);
  pg.SetGaussians(mog);
  utilities::Matrix<double> sm = pg.ComputeSimilarityMatrix();
  for(unsigned int r = 0; r < sm.NumRows(); ++r)
  {
    for(unsigned int c = 0; c < sm.NumCols(); ++c)
      std::cout<<sm(r,c)<<" ";
    std::cout<<"\n";
  }

  std::cout<<total_groups<<std::endl;
  for(unsigned int i = 0; i < group_name.size(); ++i)
    std::cout<<group_name[i]<<std::endl;


  while(fin.good() && examples < maximum_examples)
  {
    std::string line;
    std::getline(fin, line);
    utilities::TrimString(line);
    std::cout<<examples<<" "<<line<<"\n";

    std::vector<std::string> tokens;
    utilities::TokenizeString(line, ' ', tokens);
    
    if( tokens.size() == 3)
    {
      std::string feature_file = maindir + tokens[0] + suffix;
      fileutilities::SpeechFeatures sf;
      int start_frame = utilities::ToNumber<int>(tokens[1]);
      int end_frame = utilities::ToNumber<int>(tokens[2]);
      sf.ReadHtkFile(feature_file);
      utilities::Matrix<double> data = sf.frames(0, start_frame, end_frame);
      data.Transpose();
      utilities::Matrix<double> pgram = pg.ComputePosteriorgram(data);
      for(unsigned int r = 0; r < pgram.NumRows(); ++r)
        for(unsigned c = 0; c < pgram.NumCols(); ++c)
          pgram(r,c) = std::log(pgram(r,c));
      pgram_set.push_back(pgram);
      std::vector<int> path =                                                        
          acousticunitdiscovery::FindBestPath(pgram, transition, min_frames);   
      for(unsigned int i = 0; i < path.size(); ++i)
        std::cout<<group_name[path[i]]<<" ";
      std::cout<<"\n";
    }
    ++examples;
  }
  fin.close();

  std::vector<int> path = 
      acousticunitdiscovery::ApproximateViterbiSet(pgram_set, transition,     
      min_frames);
  for(unsigned int i = 0; i < path.size(); ++i)
    std::cout<<group_name[path[i]]<<" ";
  std::cout<<"\n";

  return 0;
}
