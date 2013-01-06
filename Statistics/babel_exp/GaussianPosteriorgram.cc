#include "DiagonalGaussian.h"
#include "MixtureOfDiagonalGaussians.h"
#include "HmmSet.h"
#include "ImageIO.h"
#include "SpeechFeatures.h"
#include <iostream>
#include <vector>
#include <string>
#include <fstream>

std::vector<int> ReadVector(std::string filename)
{
  std::vector<int> ret;
  std::ifstream fin;
  fin.open(filename.c_str(), std::ios::in);
  while(fin.good())
  {
    double n;
    fin >> n;
    ret.push_back(static_cast<int>(n)-1);
  }
  ret.resize( ret.size() - 1 );
  return ret;
}

fileutilities::SpeechFeatures ComputePosteriorgram(
    fileutilities::SpeechFeatures &sf, 
    std::vector<statistics::MixtureOfDiagonalGaussians> &mog,
    std::vector<int> &index,
    unsigned int width)
{
  fileutilities::SpeechFeatures ret;
  std::vector<std::vector<double> > pgram, record;
  record = sf.record(0);
  pgram.resize(record.size());
  for(unsigned int i = 0; i < record.size(); ++i)
  {
    pgram[i].resize(width, 0);
    double maxval = 0;
    for(unsigned x = 0; x < mog.size(); ++x)
    {
      pgram[i][ index[x] ] += mog[x].Likelihood(record[i]);
      if(pgram[i][ index[x] ] > maxval)
        maxval = pgram[i][ index[x] ];
    }
    if(maxval > 0)
      for(unsigned int x = 0; x < width; ++x)
        pgram[i][x] = pgram[i][x] / maxval;

  }
  ret.Initialize(pgram);

  return ret;
}

int main()
{

  statistics::HmmSet htk;
  std::string filename("/people/hartmann/research/babel/FlatStartGrapheme/hmm_training/cantonese_plp_grapheme00/hmm39/hmmdefs");
  htk.LoadHtkHmmSet(filename);
  std::vector<statistics::MixtureOfDiagonalGaussians> mog = htk.states();
  std::vector<std::vector<std::string> > names = htk.mixture_names();

  std::vector<int> cluster;
  cluster = ReadVector(std::string("idx.txt"));
  std::cout<<cluster.size()<<std::endl;

  fileutilities::SpeechFeatures sf, pgram;
  std::string indir = "/people/hartmann/research/babel/features/BP_101/scripted/training/plp/";
  std::string outdir = "/people/hartmann/research/babel/features/BP_101/scripted/training/mog_pgram/";
  filename = "/people/hartmann/research/babel/features/BP_101/scripted/training/plp/BABEL_BP_101_27619_20111102_195109_L1_scripted.htk.mvn";
  sf.ReadHtkFile(filename);

  std::vector<double> mean = sf.CalculateRecordMean(0);
  for(unsigned int i = 0; i < mean.size(); ++i)
    std::cout<<mean[i]<<std::endl;

  fileutilities::WriteBinaryPGM(sf.record(0), std::string("example.pgm"));
  pgram = ComputePosteriorgram(sf, mog, cluster, 100);
  fileutilities::WriteBinaryPGM(pgram.record(0), std::string("posteriorgram.pgm"));
  pgram.WriteHtkFile(std::string("test.htk"));

  std::string filelist = "/people/hartmann/research/babel/filelist/bp101.scripted.training.list";
  std::ifstream fin;
  fin.open(filelist.c_str(), std::ios::in);
  int count = 0;
  while(fin.good())
  {
    std::string line;
    fin >> line;
    if(line.length() > 0)
    {
     std::cout<<count<<": "<<line<<std::endl;
     filename = indir + line + std::string(".htk.mvn");
     sf.ReadHtkFile(filename);
     pgram = ComputePosteriorgram(sf, mog, cluster, 100);
     filename = outdir + line + std::string(".pgram");
     pgram.WriteHtkFile(filename);
     count++;
    }
  }
  fin.close();

  return 0;
}
