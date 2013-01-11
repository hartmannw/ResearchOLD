#include "MultiBestPath.h"
#include "SpeechFeatures.h"
#include "Matrix.h"
#include<string>
#include<vector>
#include<iostream>

void PrintVector(std::vector<double> v)
{
  for(unsigned int i = 0; i < v.size(); ++i)
    std::cout<<i<<": "<<v[i]<<std::endl;
}

int main()
{
  std::string fname;
  fileutilities::SpeechFeatures sf; 
  fname = "example1.pgram";
  int min_frames = 5;

  sf.ReadHtkFile(fname);
  utilities::Matrix<double> transition;
  transition = acousticunitdiscovery::GenerateTransitionMatrix(100, 0.5);
  utilities::Matrix<double> pgram = sf.frames(0, 201, 300);
  //utilities::Matrix<double> pgram = sf.record(0);
  pgram.Initialize(pgram.NumRows(), 20);
  pgram.Transpose();
  std::vector<utilities::Matrix<double> > pgram_set;

  //PrintVector(pgram.GetCol(4));
  //PrintVector(pgram.GetCol(5));
  //PrintVector(pgram.GetCol(6));

  std::cout<<pgram.NumRows()<<" "<<pgram.NumCols()<<std::endl;
  for(unsigned int r = 0; r < pgram.NumRows(); ++r)
    for(unsigned int c = 0; c < pgram.NumCols(); ++c)
      pgram(r,c) = std::log(pgram(r,c));
  pgram_set.push_back(pgram);

  std::vector<int> path = 
    acousticunitdiscovery::FindBestPath(pgram, transition, min_frames);

  for(unsigned int i = 0; i < path.size(); ++i)
    std::cout<<path[i]<<" ";
  std::cout<<std::endl;
  std::cout<<path.size()<<std::endl;
  
  std::vector<int> initial_path;
  //initial_path.push_back(3);
  //initial_path.push_back(15);
  //initial_path.push_back(16);

  double score;
  path = acousticunitdiscovery::FindRestrictedViterbiPath(pgram, transition, 
      min_frames, initial_path, false, score);

  for(unsigned int i = 0; i < path.size(); ++i)
    std::cout<<path[i]<<" ";
  std::cout<<std::endl;
  std::cout<<score<<std::endl;

  path = acousticunitdiscovery::ApproximateViterbiSet(pgram_set, transition, 
      min_frames);

  for(unsigned int i = 0; i < path.size(); ++i)
    std::cout<<path[i]<<" ";
  std::cout<<std::endl;
  
  return 0;
}
