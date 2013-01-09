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

  PrintVector(pgram.GetCol(6));

  std::cout<<pgram.NumRows()<<" "<<pgram.NumCols()<<std::endl;
  for(unsigned int r = 0; r < pgram.NumRows(); ++r)
    for(unsigned int c = 0; c < pgram.NumCols(); ++c)
      pgram(r,c) = std::log(pgram(r,c));

  std::vector<int> path = 
    acousticunitdiscovery::FindBestPath(pgram, transition, 3);

  for(unsigned int i = 0; i < path.size(); ++i)
    std::cout<<path[i]<<" ";
  std::cout<<std::endl;
  std::cout<<path.size()<<std::endl;
  
  std::vector<int> initial_path;
  initial_path.push_back(62);
  initial_path.push_back(12);
  initial_path.push_back(14);

  path = acousticunitdiscovery::FindRestrictedViterbiPath(pgram, transition, 3,
      initial_path);

  for(unsigned int i = 0; i < path.size(); ++i)
    std::cout<<path[i]<<" ";
  std::cout<<std::endl;

  return 0;
}
