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

  sf.ReadHtkFile(fname);
  utilities::Matrix<double> transition;
  transition = acousticunitdiscovery::GenerateTransitionMatrix(100, 0.5);
  utilities::Matrix<double> pgram = sf.record(0);
  pgram.Transpose();
  std::vector<utilities::Matrix<double> > pgram_set;
  pgram_set.push_back(pgram);
  pgram_set.push_back(pgram);

  //PrintVector(pgram.GetCol(4));
  //PrintVector(pgram.GetCol(5));
  //PrintVector(pgram.GetCol(6));

  std::cout<<pgram.NumRows()<<" "<<pgram.NumCols()<<std::endl;
  for(unsigned int r = 0; r < pgram.NumRows(); ++r)
    for(unsigned int c = 0; c < pgram.NumCols(); ++c)
      pgram(r,c) = std::log(pgram(r,c));

  std::vector<int> path = 
    acousticunitdiscovery::FindBestPath(pgram, transition, 10);

  for(unsigned int i = 0; i < path.size(); ++i)
    std::cout<<path[i]<<" ";
  std::cout<<std::endl;
  std::cout<<path.size()<<std::endl;
  
  std::vector<int> initial_path;
  initial_path.push_back(51);
  initial_path.push_back(17);
  initial_path.push_back(19);

  double score;
  path = acousticunitdiscovery::FindRestrictedViterbiPath(pgram, transition, 10,
      initial_path, true, score);

  for(unsigned int i = 0; i < path.size(); ++i)
    std::cout<<path[i]<<" ";
  std::cout<<std::endl;
  std::cout<<score<<std::endl;

  path = acousticunitdiscovery::ApproximateViterbiSet(pgram_set, transition, 
      10);

  for(unsigned int i = 0; i < path.size(); ++i)
    std::cout<<path[i]<<" ";
  std::cout<<std::endl;
  
  return 0;
}
