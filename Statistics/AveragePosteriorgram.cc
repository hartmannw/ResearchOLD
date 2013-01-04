#include "DiagonalGaussian.h"
#include "MixtureOfDiagonalGaussians.h"
#include "HmmSet.h"
#include "ImageIO.h"
#include "SpeechFeatures.h"
#include "StringFunctions.h"
#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <cmath>

struct WordLocation
{
  std::string file;
  int start;
  int end;
};
typedef struct WordLocation WordLocation;

struct ViterbiInfo
{
  double score;
  int from;
};
typedef struct ViterbiInfo ViterbiInfo;

std::vector<double> RescaleVector(std::vector<double> v, int scale);
std::vector<int> FindPronunciation(std::vector<std::vector<double> > pgram);

std::vector<std::vector<double> > LoadPosteriorgram(std::string filename, 
    int start, int end)
{
  fileutilities::SpeechFeatures pgram;
  pgram.ReadHtkFile(filename);
  return pgram.frames(0, start, end);
}

std::vector<std::vector<double> > RescalePosteriorgram(
    std::vector<std::vector<double> > pgram, int size)
{
  std::vector<std::vector<double> > ret;
  ret.resize(size);
  for(int f = 0; f < size; ++f)
    ret[f].resize(pgram[0].size());
  for(unsigned int d = 0; d < pgram[0].size(); ++d)
  {
    std::vector<double> dimension;
    for(unsigned int f = 0; f < pgram.size(); ++f)
      dimension.push_back(pgram[f][d]);
    dimension = RescaleVector(dimension, size);
    for(unsigned int i = 0; i < dimension.size(); ++i)
      ret[i][d] = dimension[i];
  }
  return ret;
}

std::vector<std::vector<double> > AveragePosteriorgram(
    std::vector<WordLocation> locations, std::string directory, 
    std::string suffix, int size)
{
  std::vector<std::vector<double> > ret;
  for(unsigned int i = 0; i < locations.size(); ++i)
  {
    std::string filename = directory + locations[i].file + suffix;
    std::cout<<filename<<std::endl;
    std::vector<std::vector<double> > pgram;
    pgram= LoadPosteriorgram(filename, locations[i].start, locations[i].end);
    //FindPronunciation(pgram);
    pgram = RescalePosteriorgram(pgram, size);
    if(i == 0)
    {
      ret = pgram;
    }
    else
    {
      for(unsigned int f = 0; f < pgram.size(); ++f)
        for(unsigned int d = 0; d < pgram[f].size(); ++d)
          ret[f][d] += pgram[f][d];
    }
  }

  for(unsigned int f = 0; f < ret.size(); ++f)
    for(unsigned int d = 0; d < ret[f].size(); ++d)
      ret[f][d] = ret[f][d] /  locations.size();
  return ret;
}

std::vector<int> FindPronunciation(std::vector<std::vector<double> > pgram)
{
  double trans = 0.99;
  std::vector<int> ret;
  std::vector<std::vector<ViterbiInfo> > dp;
  dp.resize(pgram.size());
  dp[0].resize(pgram[0].size());

  for(unsigned int f = 0; f < pgram.size(); ++f)
    for(unsigned int d = 0; d < pgram[f].size(); ++d)
    {
      pgram[f][d] = std::log(pgram[f][d]);
      if(pgram[f][d] < -15)
        pgram[f][d] = -15;
    }


  for(unsigned int d = 0; d < dp[0].size(); ++d)
    dp[0][d].score = pgram[0][d];
  for(unsigned int f = 1; f < pgram.size(); ++f)
  {
    dp[f].resize(pgram[f].size());
    for(unsigned int d = 0; d < pgram[f].size(); ++d)
    {
      ViterbiInfo best = {-1000000, 0};
      for(unsigned int i = 0; i < pgram[f-1].size(); ++i)
      {
        double score = dp[f-1][i].score;
        if( i == d )
          score += log(trans);
        else
          score += log(1 - trans);
        if(score > best.score)
        {
          best.score = score;
          best.from = i;
        }
      }
      best.score += pgram[f][d];
      dp[f][d] = best;
    }
  }
  std::cout<<"Beginning backtrack"<<std::endl;
  // Find our starting point
  ViterbiInfo best = dp[dp.size()-1][0];
  std::cout<<best.score<<" "<<best.from<<std::endl;
  int index = 0;
  for(unsigned int i = 1; i < dp[dp.size()-1].size(); ++i)
  {
    if(dp[dp.size()-1][i].score > best.score)
    {
      best = dp[dp.size()-1][i];
      index = i;
    }
  }
  ret.push_back(index);
  for(int f = dp.size() -2; f >= 0; --f)
  {
    ret.push_back(best.from);
    std::cout<<best.from<<" ";
    best = dp[f][best.from];
  }
  std::cout<<std::endl;
  return ret;
}

std::vector<double> RescaleVector(std::vector<double> v, int scale)
{
  std::vector<double> ret;
  double step = static_cast<double>(v.size() - 1) / (scale - 1);
  for(int i = 0; i < (scale - 1); ++i)
  {
    double point = i * step;
    int lind = static_cast<int>( std::floor(point) );
    int rind = static_cast<int>( std::ceil(point) );
    if( lind == rind)
    {
      ret.push_back( v[lind] );
    }
    else
    {
      double m = (v[lind] - v[rind]) / (lind - rind);
      double b = v[lind] - (m * lind);
      ret.push_back(m * point + b);
    }
  }
  ret.push_back(v[ v.size() - 1]);
  return ret;
}

std::vector<WordLocation> FindWordLocations(std::string filename, 
    std::string word)
{
  std::vector<WordLocation> ret;
  std::ifstream fin;
  fin.open(filename.c_str(), std::ios::in);
  int filecount = 0;
  std::string line;
  getline(fin, line);
  while(fin.good())
  {
    filecount++;
    getline(fin, line);
    if(line.length() > 1)
    {
     
     line = utilities::TrimString(line);
     std::string name = line.substr(3, line.size() - 12);
     getline(fin, line);
     line = utilities::TrimString(line);
     while(line != std::string(".") )
     {
      std::vector<std::string> tokens;
      utilities::TokenizeString(line, ' ', tokens);
      if(tokens.size() > 2 && tokens[2] == word)
      {
        int start = std::floor(utilities::stoi(tokens[0]) / 100000);
        int end = std::floor(utilities::stoi(tokens[1]) / 100000) - 1;
        WordLocation wl = {name, start, end};
        ret.push_back(wl);
      }
      getline(fin, line);
      line = utilities::TrimString(line);
     }
    }
  }
  return ret;
}

int main()
{
  std::vector<double> a;
  a.resize(5, 0);
  a[0] = 0; 
  a[1] = 1; 
  a[2] = 2;
  a[3] = 4;
  a[4] = 8;

  for(unsigned int i = 0; i < a.size(); ++i)
    std::cout<<a[i]<<std::endl;

  a = RescaleVector(a, 17);
  std::cout<<std::endl;

  for(unsigned int i = 0; i < a.size(); ++i)
    std::cout<<a[i]<<std::endl;

  std::string word("BOU2WU6");
  std::string filename("/people/hartmann/research/babel/FlatStartGrapheme/testing/align_grapheme.mlf");
  std::string directory("/people/hartmann/research/babel/features/BP_101/scripted/training/mog_pgram/");
  std::string suffix(".pgram");
  std::vector<WordLocation> locations;
  locations = FindWordLocations(filename, word);
  std::vector<std::vector<double> > pgram;
  pgram = AveragePosteriorgram(locations, directory, suffix, 100);
  fileutilities::WriteBinaryPGM(pgram, std::string("average_bou2wu6.pgm"));

  return 0;
}
