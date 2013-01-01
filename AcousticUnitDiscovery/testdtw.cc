#include "SpeechFeatures.h"
#include "DynamicTimeWarp.h"
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

  //fname1 = "utt1.cep.mvn";
  //fname2 = "utt2.cep.mvn";
  fname1 = "../../babel/features/BP_101/scripted/training/plp/BABEL_BP_101_11422_20111019_143344_D2_scripted.htk.mvn";
  fname2 = "../../babel/features/BP_101/scripted/training/plp/BABEL_BP_101_12003_20111116_130201_D2_scripted.htk.mvn";
  fname2 = "../../babel/features/BP_101/scripted/training/plp/BABEL_BP_101_32890_20111130_215110_S8_scripted.htk.mvn";
  fname_sil = "../../babel/features/BP_101/scripted/training/silence/BABEL_BP_101_11422_20111019_143344_D2_scripted.sil";
  std::vector<bool> silence;                                                
  silence = ReadSilence(fname_sil); 


  sf1.ReadHtkFile(fname1);
  sf2.ReadHtkFile(fname2);

for(int i = 0; i < 1; i++)
{
  std::cout<<"Iteration "<<i<<"\n";
  acousticunitdiscovery::DynamicTimeWarp dtw;
  dtw.set_utterance_one(sf1.record(0));
  dtw.set_utterance_two(sf2.record(0));


  dtw.ComputeSimilarityMatrix();
  //dtw.ComputeStandardDTW();
  dtw.ComputeSegmentalDTW(20);
  //dtw.SaveResultAsPGM( std::string("result.pgm"));
  dtw.IncreaseSilenceCost(silence);
  dtw.PrunePathsByLCMA(50, 0.1);
  dtw.SaveResultAsPGM( std::string("result_pruned.pgm"));
}
  return 0;
}
