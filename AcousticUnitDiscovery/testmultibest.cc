#include "MultiBestPath.h"
#include "SpeechFeatures.h"
#include "Matrix.h"
#include<string>
#include<vector>
#include<iostream>

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
  std::vector<int> path = 
    acousticunitdiscovery::FindBestPath(pgram, transition, 5);

  for(unsigned int i = 0; i < path.size(); ++i)
    std::cout<<path[i]<<" ";
  std::cout<<std::endl;
  return 0;
}
