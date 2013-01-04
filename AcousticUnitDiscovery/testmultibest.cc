#include "MultiBestPath.h"
#include "SpeechFeatures.h"
#include<string>
#include<vector>
#include<iostream>

int main()
{
  std::string fname;
  fileutilities::SpeechFeatures sf; 
  fname = "../../babel/features/BP_101/scripted/training/mog_pgram/BABEL_BP_101_11422_20111019_143344_D2_scripted.pgram";
  fname = "../../babel/features/BP_101/scripted/training/plp/BABEL_BP_101_11422_20111019_143344_D2_scripted.htk.mvn";


  sf.ReadHtkFile(fname);
  std::vector<std::vector<double> > transition;
  transition = acousticunitdiscovery::GenerateTransitionMatrix(0.5, 100);
  return 0;
}
