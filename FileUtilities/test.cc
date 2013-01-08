#include "SpeechFeatures.h"
#include "ImageIO.h"
#include<iostream>


int main()
{
  fileutilities::SpeechFeatures sf;
  std::string filename = "test_data/example.htk2";
  std::string fileout = "test_data/example.htk2";
  sf.ReadHtkFile(filename);

  std::vector<double> mean = sf.CalculateRecordMean(0);
  for(int i = 0; i < mean.size(); ++i)
    std::cout<<mean[i]<<std::endl;
  sf.WriteHtkFile(fileout);

  fileutilities::WriteBinaryPGM(sf.record(0).GetVectorOfVectors(), 
      std::string("test_data/example.pgm"));

  return 0;
}

