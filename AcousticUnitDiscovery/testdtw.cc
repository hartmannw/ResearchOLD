#include "SpeechFeatures.h"
#include "DynamicTimeWarp.h"
#include "Matrix.h"
#include <vector>
#include <iostream>
#include <string>

std::vector<bool> ReadSilence(std::string filename)                             
{                                                                               
  std::ifstream fin(filename.c_str());                                          
  std::string line, one = "1";                                                  
  getline(fin, line);                                                           
  std::vector<bool> ret;                                                        
                                                                                
  for(unsigned int i = 0; i < line.length(); ++i)                               
  {                                                                             
    ret.push_back( line[i] == '1');                                             
  }                                                                             
  return ret;                                                                   
}

int main()
{
  fileutilities::SpeechFeatures sf1, sf2;
  std::string fname1, fname2, fname_sil;
  //acousticunitdiscovery::DynamicTimeWarp dtw;

  fname1 = "example1.pgram";
  //fname1 = "../FileUtilities/example.htk";
  fname2 = "example2.pgram";
  //fname_sil = "../../babel/features/BP_101/scripted/training/silence/BABEL_BP_101_11422_20111019_143344_D2_scripted.sil";
  std::vector<bool> silence;                                                
  //silence = ReadSilence(fname_sil); 

  sf1.ReadHtkFile(fname1);
  sf2.ReadHtkFile(fname2);

  //utilities::PrintMatrix(sf1.record(0), std::string("%0.2f"));
for(int i = 0; i < 1; i++)
{
  std::cout<<"Iteration "<<i<<"\n";
  acousticunitdiscovery::DynamicTimeWarp dtw;
  utilities::Matrix<double> r1, r2;
  //r1.Initialize(sf1.record(0));
  //r2.Initialize(sf2.record(0));
  r1 = sf1.record(0);
  r2 = sf2.record(0);
  dtw.set_utterance_one(r1);
  dtw.set_utterance_two(r2);


  dtw.ComputeSimilarityMatrix();
  //dtw.ComputeStandardDTW();
  dtw.ComputeSegmentalDTW(20);
  //dtw.SaveResultAsPGM( std::string("result.pgm"));
  //dtw.IncreaseSilenceCost(silence);
  dtw.PrunePathsByLCMA(50, 0.1);
  dtw.SaveResultAsPGM( std::string("result_pruned.pgm"));
}
  return 0;
}
