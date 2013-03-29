// William Hartmann (hartmannw@gmail.com)
// This is free and unencumbered software released into the public domain.
// See the UNLICENSE file for more information.

// Used for generating pronunciations for a list of words.

#include<cmath>
#include<iostream>
#include<fstream>
#include<vector>
#include<string>
#include<assert.h>

#include "Matrix.h"
#include "SpeechFeatures.h"
#include "PosteriorgramGenerator.h"
#include "StringFunctions.h"
#include "HmmSet.h"
#include "MultiBestPath.h"
#include "HmmDecoder.h"


// Stores all the basic information that is needed for generating 
// pronunciations. Most of the settings can be set by the command line. By
// keeping them together in a structure, parameter lists are shorter and code is
// cleaner.
typedef struct
{
  std::string word_information;
  std::string locationdir;
  std::string hmmfile;
  std::string monohmmfile;
  std::string cluster_file;
  std::string datadir;
  std::string suffix;
  unsigned int max_examples;
  unsigned int min_examples;
  double length_variance;
  unsigned int pronunciation_type;

  // The following parameters determined by the data.
  unsigned int dimension;
  unsigned int total_clusters;
} PronunciationParameters;

// Stores basic information about the location of a segment in a speech file.
typedef struct
{
  std::string file; // File prefix.
  unsigned int start; // Frame where the word begins.
  unsigned int end; // Frame where the word ends.
} WordLocation;


std::vector<std::string> SplitTriphone(std::string s)
{
  std::vector<std::string> ret, tokens;
  utilities::TokenizeString(s, '-', tokens);
  if(tokens.size() > 1)
  {
    ret.push_back(tokens[0]);
    s = tokens[1];
  }
  else
    ret.push_back("");
  tokens.clear();

  utilities::TokenizeString(s, '+', tokens);
  ret.push_back(tokens[0]);
  if(tokens.size() > 1)
    ret.push_back(tokens[1]);
  else
    ret.push_back("");
  return ret;

}

std::vector<int> ReadTriphoneLabels(std::string filename, 
    statistics::HmmSet &htk)
{
  std::vector<int> ret;
  std::ifstream fin;
  fin.open(filename.c_str(), std::ios::in);
  std::vector<statistics::HiddenMarkovModel> hmmset = htk.hmms();
  for(unsigned int i = 0; i < hmmset.size(); ++i)
  {
    std::vector<std::string> tri = SplitTriphone(htk.Hmm(i).name());
    if( tri[0].length() < 1 || tri[2].length() < 1 )
      ret.push_back(0);
    else
    {
      double n;
      fin >> n;
      ret.push_back(static_cast<int>(n) - 1);
    }
  }
  return ret;
}

// Note this only works if we are dealing with clustered graphemes.
void CorrectIndependentClusters( std::vector<int> &cluster, statistics::HmmSet &htk)
{
  std::vector<utilities::Matrix<int> > lcounts, rcounts;
  std::vector<std::vector<int> > ccounts;
  int labels = 0;
  //Find number of labels
  for(unsigned int i = 0; i < cluster.size(); ++i)
    if(cluster[i] > labels)
      labels = cluster[i];
  labels++;
  lcounts.resize(labels);
  rcounts.resize(labels);
  ccounts.resize(labels);
  for(int i = 0; i < labels; ++i)
  {
    lcounts[i].Initialize(26,26,0);
    rcounts[i].Initialize(26,26,0);
    ccounts[i].resize(26,0);
  }
  for(unsigned int i = 0; i < cluster.size(); ++i)
  {
    std::vector<std::string> tri = SplitTriphone(htk.Hmm(i).name());
    if(tri[0].length() > 0 && tri[2].length() > 0)
    {
      rcounts[cluster[i]](tri[1][0] - 'a', tri[2][0] - 'a')++;
      lcounts[cluster[i]](tri[1][0] - 'a', tri[0][0] - 'a')++;
      ccounts[cluster[i]][ tri[1][0] - 'a' ]++;
    }
  }
  // Now we modify the actual clusters
  for(unsigned int i = 0; i < cluster.size(); ++i)
  {
    std::vector<std::string> tri = SplitTriphone(htk.Hmm(i).name());
    // First find best context independent label
    int ind = 0;
    unsigned int c = tri[1][0] - 'a';
    for(unsigned int j = 0; j < rcounts.size(); ++j)
      if(ccounts[j][c] > ccounts[ind][c])
        ind=j;
    if(tri[1] == "sil" || tri[1] == "sp")
    {
      //Do nothing
    }
    else if(tri[0].length() < 1 && tri[2].length() < 1)
    {
      cluster[i] = ind;
    }
    else if(tri[0].length() < 1)
    {
      assert(tri[1][0] - 'a' >= 0);
      assert(tri[1][0] - 'a' < 26);
      int maxind = cluster[i];
      unsigned int r = tri[1][0] - 'a';
      unsigned int c = tri[2][0] - 'a';
      for(unsigned int j = 0; j < rcounts.size(); ++j)
        if(rcounts[j](r,c) > rcounts[maxind](r,c))
          maxind = j;
      cluster[i] = maxind;
      if(rcounts[maxind](r,c) == 0)
        cluster[i] = ind;
    }
    else if(tri[2].length() < 1)
    {
      int maxind = cluster[i];
      unsigned int r = tri[1][0] - 'a';
      unsigned int c = tri[0][0] - 'a';
      for(unsigned int j = 0; j < lcounts.size(); ++j)
        if(lcounts[j](r,c) > lcounts[maxind](r,c))
          maxind = j;
      cluster[i] = maxind;
      if(lcounts[maxind](r,c) == 0)
        cluster[i] = ind;
    }
  }

  return;
}

void RemoveSilenceFromLM(utilities::Matrix<double> &lm, statistics::HmmSet &htk)
{
  for(unsigned int h = 0; h < htk.HmmCount(); ++h)
  {
    if(htk.Hmm(h).name() == "sp" || htk.Hmm(h).name() == "sil")
      for(unsigned int i = 0; i < lm.NumRows(); ++i)
      {
        lm(h,i) = -1000000;
        lm(i,h) = -1000000;
      }
  }
}

std::vector<WordLocation> LoadLocationList( const std::string &location_file )
{
  std::ifstream location_fin;
  location_fin.open(location_file.c_str());
  std::vector<WordLocation> locations;
  while( location_fin.good() ) // Load the information about each word 
  {                            // location into a list.
    std::string loc_line;
    std::getline(location_fin, loc_line);
    if(loc_line.length() > 0)
    {
      utilities::TrimString(loc_line);
      std::vector<std::string> loc_tokens;
      utilities::TokenizeString(loc_line, ' ', loc_tokens);
      WordLocation loc;
      loc.file = loc_tokens[0];
      loc.start = utilities::ToNumber<unsigned int>(loc_tokens[1]);
      loc.end = utilities::ToNumber<unsigned int>(loc_tokens[2]);
      locations.push_back(loc);
    } // end if loc_line.length() > 0
  } // end while location_fin.good()
  location_fin.close();

  return locations;
}

std::vector<WordLocation> PruneLocationList(
    const std::vector<WordLocation> &locations, double length_variance)
{
  std::vector<WordLocation> ret;
  if( locations.size() == 0)
    return ret;
  double mean_length = 0;
  for(unsigned int i = 0; i < locations.size(); ++i)
    mean_length += (locations[i].end - locations[i].start + 1);
  mean_length = mean_length / locations.size();
  double min_length = (1 - length_variance) * mean_length;
  double max_length = (1 + length_variance) * mean_length;
  for(unsigned int i = 0; i < locations.size(); ++i)
  {
    int length = locations[i].end - locations[i].start + 1;
    if( (length >= min_length) && (length <= max_length) )
      ret.push_back(locations[i]);
  }
  return ret;
}

bool LoadDataSegment( const WordLocation &location, 
    const PronunciationParameters &param, utilities::Matrix<double> &data)
{ 
  std::string feature_file = 
      param.datadir + "/" + location.file + "." + param.suffix;
  fileutilities::SpeechFeatures sf;
  sf.ReadHtkFile(feature_file);
  data = sf.frames(0, location.start, location.end);
  data.Transpose();
  return true;
}

std::vector<std::string> ConvertToAlphaPronunciation(
    const std::vector<int> &pronunciation)
{
  std::vector<std::string> ret;
  std::string last("");
  for(unsigned int i = 0; i < pronunciation.size(); ++i)
  {
    if(pronunciation[i] < 0)
    {
      std::cout<<"ERROR"<<std::endl;
      exit(1);
    }
    std::string phone("aa");
    phone[0] += std::floor( static_cast<double>(pronunciation[i]) / 10);
    phone[1] += pronunciation[i] % 10;
    //if( phone != last)
      ret.push_back(phone);
    last = phone;
  }
  return ret;
}

std::vector<std::string> ConvertToHmmNamePronunciation(
    const std::vector<int> &pronunciation, const statistics::HmmSet &htk)
{
  std::vector<std::string> ret;
  std::string last("");
  for(unsigned int i = 0; i < pronunciation.size(); ++i)
  {
    if(pronunciation[i] < 0)
    {
      std::cout<<"ERROR"<<std::endl;
      exit(1);
    }
    std::string phone = htk.Hmm(pronunciation[i]).name();
    //if( phone != last)
      ret.push_back(phone);
    last = phone;
  }
  return ret;
}

std::vector<int> BestPathInSet(
    const std::vector<utilities::Matrix<double> > &score_set,
    statistics::HmmDecoder hdecode)
{
  std::vector<std::vector<int> > path_set;
  std::vector<double> path_score;
  double score;
  for(unsigned int i = 0; i < score_set.size(); ++i)
  {
    std::vector<int> path;
    path = hdecode.ViterbiPath(score_set[i]);
    path_set.push_back(path);
    double total_score = 0;
    // Check the score for every example in pgram_set
    for(unsigned int j = 0; j < score_set.size(); ++j)
    {
      score = hdecode.ForceAlignScore(score_set[j], path);
      total_score += score;
    }
    path_score.push_back(total_score);
    //for(unsigned int x = 0; x < path.size(); ++x)
    //  std::cout<<" "<<path[x];
    //std::cout<<" "<<total_score<<std::endl;
  }

  unsigned int best_index = 0;
  for(unsigned int i = 1; i < path_score.size(); ++i)
    if(path_score[best_index] < path_score[i])
      best_index = i;

  return path_set[best_index];
}

std::vector<std::string> GenerateTriphonePronunciation(
    const std::vector<std::string> &pronunciation)
{
  if(pronunciation.size() == 1)
    return pronunciation;
  std::vector<std::string> ret;
  std::string plus("+");
  std::string minus("-");
  unsigned int last = pronunciation.size() - 1;
  ret.push_back(pronunciation[0] + plus + pronunciation[1]);
  for(unsigned int i = 1; i < last; ++i)
    ret.push_back(pronunciation[i-1] + minus + pronunciation[i] + plus + 
        pronunciation[i+1]);
  ret.push_back(pronunciation[last-1] + minus + pronunciation[last]);
  return ret;
}

std::vector<int> GenerateClusteredPronunciation(
    const std::vector<std::string> &triphone_pronunciation,
    const statistics::HmmSet &htk, const std::vector<int> &cluster_index, 
    const PronunciationParameters &param)
{
  std::vector<int> ret;
  for(unsigned int p = 0; p < triphone_pronunciation.size(); ++p)
  {
    unsigned int index = htk.HmmIndex(triphone_pronunciation[p]);
    ret.push_back(cluster_index[index]);
  }
  return ret;
}

int main(int argc, char* argv[])
{
  if( argc < 11)
  {
    std::cout<<"Usage is <Word Information> <Location Directory>"<<
        " <Full HMM File> <Monophone HMM File>"<<
        " <Cluster Index> <Data Directory>"<<
        " <File Suffix> <Max Examples> <Min Examples> <Length Variance>"<<
        std::endl;
    exit(0);
  }
  PronunciationParameters param;
  param.word_information = std::string(argv[1]);
  param.locationdir = std::string(argv[2]);
  param.hmmfile = std::string(argv[3]);
  param.monohmmfile = std::string(argv[4]);
  param.cluster_file = std::string(argv[5]);
  param.datadir = std::string(argv[6]);
  param.suffix = std::string(argv[7]);
  param.max_examples = utilities::ToNumber<unsigned int>(std::string(argv[8]));
  param.min_examples = utilities::ToNumber<unsigned int>(std::string(argv[9]));
  param.length_variance = utilities::ToNumber<double>(std::string(argv[10]));
  param.pronunciation_type = 1; // 0=clustered, 1=DataDriven !Change to Enum!

  statistics::HmmSet htk, monohtk;
  htk.LoadHtkHmmSet(param.hmmfile);
  monohtk.LoadHtkHmmSet(param.monohmmfile);
  utilities::Matrix<double> lm;
  lm.Initialize(monohtk.HmmCount(), monohtk.HmmCount(), std::log(0.05));
  RemoveSilenceFromLM(lm, monohtk);
  std::vector<statistics::MixtureOfDiagonalGaussians> mog = htk.states();
  // Ugly way of initializing dimension, Find a better way!!!
  param.dimension = mog[0].gaussian(0).dimension();
  std::vector<int> cluster_index = ReadTriphoneLabels(param.cluster_file, htk);
  CorrectIndependentClusters(cluster_index, htk);
  param.total_clusters = cluster_index.size();
  std::ifstream wordlist_fin;
  wordlist_fin.open(param.word_information.c_str());
  
  if( !wordlist_fin.good())
  {
    std::cout<<"File "<<param.word_information<<" could not be opened.\n";
    exit(1);
  }
  // Loop through each word in the word information file
  while( wordlist_fin.good())
  {
    std::string info_line;
    std::getline(wordlist_fin, info_line);
    if( info_line.length() > 0) // Skip blank lines.
    {
      utilities::TrimString(info_line);
      std::vector<std::string> info_tokens, original_pronunciation;
      utilities::TokenizeString(info_line, ' ', info_tokens);

      if(info_tokens.size() < 3) // Need at least 3 tokens.
      {
        std::cout<<"Current line does not contain enough information: "<<
            info_line<<std::endl;
        exit(1);
      }
      std::string word = info_tokens[0];
      std::cout<<word;
      std::string location_file = param.locationdir + "/" + info_tokens[1];
      for(unsigned int i = 2; i < info_tokens.size(); ++i)
        original_pronunciation.push_back(info_tokens[i]);

      std::vector<WordLocation> locations = LoadLocationList(location_file);
      locations = PruneLocationList(locations, param.length_variance);
      std::cout<<" "<<locations.size();

      if( param.pronunciation_type == 1 && 
          locations.size() > param.min_examples)
      {
        statistics::HmmDecoder hdecode;
        hdecode.Initialize(&monohtk, &lm);
        std::vector<int> mono_index_pronunciation;
        std::vector<std::string> mono_pronunciation;
        std::vector<utilities::Matrix<double> > score_set;

        for(unsigned int i = 0; 
            (i < locations.size()) && (i < param.max_examples); ++i)
        {
          utilities::Matrix<double> data, scores;
          LoadDataSegment(locations[i], param, data);
          hdecode.FeaturesStateScores(data, scores);
          score_set.push_back(scores);
        }

        mono_index_pronunciation = BestPathInSet(score_set, hdecode);
        mono_pronunciation = ConvertToHmmNamePronunciation( 
            mono_index_pronunciation, monohtk);
      
        for(unsigned int i = 0; i < mono_pronunciation.size(); ++i)
          std::cout<<" "<<mono_pronunciation[i];
        std::cout<<std::endl;
      }
      else
      {
        std::vector<std::string> final_pronunciation;
        std::vector<int> index_pronunciation;
        std::vector<std::string> triphone_pronunciation = 
          GenerateTriphonePronunciation(original_pronunciation);
        
        index_pronunciation = GenerateClusteredPronunciation(
            triphone_pronunciation, htk, cluster_index, param);
        
        final_pronunciation = ConvertToAlphaPronunciation(index_pronunciation);
        for(unsigned int i = 0; i < final_pronunciation.size(); ++i)
          std::cout<<" "<<final_pronunciation[i];
        std::cout<<std::endl;
      }
    } // end info_line.length() > 0
  } // end wordlist_fin.good()

  wordlist_fin.close();

  return 0;
}
