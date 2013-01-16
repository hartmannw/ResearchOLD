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

int main()
{
  std::string listfile("/people/hartmann/research/SegmentalModel/tidigit_exp/word_locations/one.train.list");
  std::string suffix(".htk");
  std::string maindir("/people/hartmann/research/SegmentalModel/tidigit_exp/");
  std::string hmmfile("/people/hartmann/research/SegmentalModel/tidigit_exp/hmm_plp_phone/hmm29/hmmdefs");
  int maximum_examples = 10;
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
  pg.SetGaussians(mog);
  utilities::Matrix<double> transition;                                          
  transition = acousticunitdiscovery::GenerateTransitionMatrix(mog.size(), 0.9);

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
        std::cout<<names[path[i]][0]<<" ";
      std::cout<<"\n";
    }
    ++examples;
  }
  fin.close();

  std::vector<int> path = 
      acousticunitdiscovery::ApproximateViterbiSet(pgram_set, transition,     
      min_frames);
  for(unsigned int i = 0; i < path.size(); ++i)
    std::cout<<names[path[i]][0]<<" ";
  std::cout<<"\n";

  return 0;
}
