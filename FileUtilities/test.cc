#include "SpeechFeatures.h"
#include "ImageIO.h"
#include<iostream>


int main()
{
  fileutilities::SpeechFeatures sf;
 // std::string filename = "/people/hartmann/research/SegmentalModel/tidigit_exp/data/htkplp/train/man/ae/12a.mlp";
  std::string filename = "/people/hartmann/research/src/AcousticUnitDiscovery/utt1.cep.mvn";
  std::string fileout = "example.htk2";
  sf.ReadCepFile(filename);

  std::vector<double> mean = sf.CalculateRecordMean(0);
  for(int i = 0; i < mean.size(); ++i)
    std::cout<<mean[i]<<std::endl;
  sf.WriteHtkFile(fileout);

  fileutilities::WriteBinaryPGM(sf.record(0), std::string("example.pgm"));

  return 0;
}

