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
  fileutilities::SpeechFeatures sf; 
  std::vector<std::string> fname;
  fname.push_back("example1.pgram");
  fname.push_back("example2.pgram");
  fname.push_back("example3.pgram");
  fname.push_back("example4.pgram");
  fname.push_back("example5.pgram");
  int min_frames = 3;
  utilities::Matrix<double> transition;
  transition = acousticunitdiscovery::GenerateTransitionMatrix(25, 0.9);
  std::vector<utilities::Matrix<double> > pgram_set;

  for(unsigned int i = 0; i < fname.size(); ++i)
  {
    sf.ReadHtkFile(fname[i]);
    utilities::Matrix<double> pgram = sf.frames(0, 201, 300);
    pgram.Initialize(pgram.NumRows(), 25);
    pgram.Transpose();

    for(unsigned int r = 0; r < pgram.NumRows(); ++r)
      for(unsigned int c = 0; c < pgram.NumCols(); ++c)
        pgram(r,c) = std::log(pgram(r,c));
    pgram_set.push_back(pgram);
  }

  std::vector<int> path = 
    acousticunitdiscovery::FindBestPath(pgram_set[0], transition, min_frames);

  for(unsigned int i = 0; i < path.size(); ++i)
    std::cout<<path[i]<<" ";
  std::cout<<std::endl;
  std::cout<<path.size()<<std::endl;
  
  std::vector<int> initial_path;
  //initial_path.push_back(3);
  //initial_path.push_back(15);
  //initial_path.push_back(16);

  double score;
  path = acousticunitdiscovery::FindRestrictedViterbiPath(pgram_set[0], transition, 
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
