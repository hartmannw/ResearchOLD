// William Hartmann (hartmannw@gmail.com)
// This is free and unencumbered software released into the public domain.
// See the UNLICENSE file for more information.

// Used for generating pronunciations for a specific word. At the moment it 
// works from feature files. In the future I would like to change it so that it
// works from Posteriorgrams and they are generatedin a separate, earlier step.

#include<cmath>
#include<iostream>
#include<fstream>
#include<vector>
#include<string>
#include<random>

#include "Matrix.h"
#include "SpeechFeatures.h"
#include "PosteriorgramGenerator.h"
#include "StringFunctions.h"
#include "HmmSet.h"
#include "MultiBestPath.h"

// Stores basic information about the location of a segment in a speech file.
typedef struct
{
  std::string file; // File prefix.
  unsigned int start; // Frame where the word begins.
  unsigned int end; // Frame where the word ends.
} WordLocation;

std::vector<int> ReadVectorFromFile(std::string filename)
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

std::vector<std::string> ConvertToAlphaPronunciation(
    const std::vector<int> &pronunciation)
{
  std::vector<std::string> ret;
  std::string last("");
  for(unsigned int i = 0; i < pronunciation.size(); ++i)
  {
    std::string phone("aa");
    phone[0] += std::floor( static_cast<double>(pronunciation[i]) / 10);
    phone[1] += pronunciation[i] % 10;
    if( phone != last)
      ret.push_back(phone);
    last = phone;
  }
  return ret;
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

std::vector<std::vector<double> > GenerateSampleData( 
    std::vector<std::string> &triphone_pronunciation, statistics::HmmSet &htk,
    std::vector<statistics::MixtureOfDiagonalGaussians> &mog,
    std::default_random_engine &generator)
{
  std::vector<std::vector<double> >ret;
  std::uniform_real_distribution<double> distribution(0.0,1.0);
  unsigned int MAX_FRAMES_PER_HMM = 15;
  for(unsigned int p = 0; p < triphone_pronunciation.size(); ++p)
  {
    statistics::HiddenMarkovModel hmm = htk.Hmm(triphone_pronunciation[p]);
    unsigned int h = 0;
    unsigned int count = 0;
    while(h < hmm.NumberOfStates() && count < MAX_FRAMES_PER_HMM)
    {
      count++;
      double trans = distribution(generator);
      if( trans < hmm.transition(h+1, h+1) )
      {
        std::vector<double> frame = mog[hmm.state(h)].sample(generator);
        ret.push_back(frame);
      }
      else
      {
        ++h;
      }
    }
  }
  return ret;
}

int main(int argc, char* argv[])
{
  if( argc < 11)
  {
    std::cout<<"Usage is <Word Information> <Location Directory> <HMM File>"<<
        " <Cluster Index> <Data Directory>"<<
        " <File Suffix> <Max Examples> <Min Examples> <Length Variance>"<<
        " <Self Transition Probability>"<<std::endl;
    exit(0);
  }
  std::string word_information = std::string(argv[1]);
  std::string locationdir = std::string(argv[2]);
  std::string hmmfile = std::string(argv[3]);
  std::string cluster_file = std::string(argv[4]);
  std::string datadir = std::string(argv[5]);
  std::string suffix = std::string(argv[6]);
  unsigned int max_examples = 
      utilities::ToNumber<unsigned int>(std::string(argv[7]));
  unsigned int min_examples = 
      utilities::ToNumber<unsigned int>(std::string(argv[8]));
  double length_variance = utilities::ToNumber<double>(std::string(argv[9]));
  double self_transition = utilities::ToNumber<double>(std::string(argv[10]));
  std::default_random_engine generator; // A single RNG.
  unsigned int min_frames = 3; // Magic number, but this matches the number of 
                               // states in a typical HMM. I don't see this 
                               // changing any time soon.

  statistics::HmmSet htk;
  htk.LoadHtkHmmSet(hmmfile);
  std::vector<statistics::MixtureOfDiagonalGaussians> mog = htk.states();
  // Ugly way of initializing dimension, Find a better way!!!
  unsigned int dimension = mog[0].gaussian(0).dimension();
  statistics::PosteriorgramGenerator pg;
  utilities::Matrix<double> transition;
  std::vector<int> cluster_index = ReadVectorFromFile(cluster_file);
  unsigned int total_clusters = cluster_index.size();
  transition = acousticunitdiscovery::GenerateTransitionMatrix(
      total_clusters, self_transition);
  pg.SetGaussians(mog, cluster_index);
  
  std::ifstream wordlist_fin;
  wordlist_fin.open(word_information.c_str());

  if( !wordlist_fin.good())
  {
    std::cout<<"File "<<word_information<<" could not be opened.\n";
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
      std::string location_file = locationdir + "/" + info_tokens[1];
      for(unsigned int i = 2; i < info_tokens.size(); ++i)
        original_pronunciation.push_back(info_tokens[i]);

      std::ifstream location_fin;
      location_fin.open(location_file.c_str());
      std::vector<WordLocation> locations;
      while( location_fin.good() )
      {
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
      locations = PruneLocationList(locations, length_variance);

      std::vector<std::string> final_pronunciation;
      std::vector<int> index_pronunciation;
      if(locations.size() < min_examples)
      {
        std::vector<std::string> triphone_pronunciation =
            GenerateTriphonePronunciation(original_pronunciation);
        utilities::Matrix<double> data, modified_transition, transition_count;
        modified_transition = transition;
        for(unsigned int r = 0; r < modified_transition.NumRows(); ++r)
          for(unsigned int c = 0; c < modified_transition.NumCols(); ++c)
            modified_transition(r,c) = std::exp(modified_transition(r,c));
        transition_count.Initialize(total_clusters, total_clusters, 1);
        data.Initialize(dimension, triphone_pronunciation.size() * 3);
        std::vector<std::vector<double> > sample_data = 
            GenerateSampleData(triphone_pronunciation, htk, mog, generator);
        for(unsigned int p = 0; p < triphone_pronunciation.size(); ++p)
        {
          statistics::HiddenMarkovModel hmm = 
              htk.Hmm(triphone_pronunciation[p]);
          for(unsigned int h = 0; h < hmm.NumberOfStates(); ++h)
          {
            index_pronunciation.push_back( cluster_index[ hmm.state(h) ]);
            //std::vector<double> mean = mog[hmm.state(h)].WeightedMean();
            std::vector<double> mean = mog[hmm.state(h)].sample(generator);
            for(unsigned int m = 0; m < mean.size(); ++m)
              data(m, (p*3)+h) = mean[m];
          }
          // Handle transition logic
          for(unsigned int h = 0; h < hmm.NumberOfStates() - 1; ++h)
          {
            unsigned int from_index = cluster_index[ hmm.state(h) ];
            unsigned int to_index = cluster_index[ hmm.state(h+1) ];
            modified_transition(from_index, from_index) += 
                hmm.transition(h+1, h+1);
            transition_count(from_index, from_index)++;
            modified_transition(from_index, to_index) += 
                hmm.transition(h+1, h+2);
            transition_count(from_index, to_index)++;
          }
          unsigned int last_hmm = hmm.NumberOfStates() - 1;
          unsigned int from_index = cluster_index[ hmm.state(last_hmm) ];
          modified_transition(from_index, from_index) += 
              hmm.transition(last_hmm+1, last_hmm+1);
          transition_count(from_index, from_index)++;
          for(unsigned int to_index = 0; to_index < total_clusters; ++to_index)
          {
            if( to_index != from_index )
            {
              modified_transition(from_index, to_index) += 
                  hmm.transition(last_hmm+1, last_hmm+2);
              transition_count(from_index, to_index)++;
            }
          }
          // End Transistion Logic
        }
        // Normalize Modified Transition Matrix
        for(unsigned int r = 0; r < modified_transition.NumRows(); ++r)
          for(unsigned int c = 0; c < modified_transition.NumCols(); ++c)
            modified_transition(r,c) = std::log(modified_transition(r,c) / 
                transition_count(r,c));
        // Using data vector generate a posteriorgram.
        //data.Initialize(sample_data);
        //data.Transpose();
        utilities::Matrix<double> pgram = pg.ComputePosteriorgram(data);
        //for(unsigned int p = 0; p < index_pronunciation.size(); ++p)
        //  std::cout<<" "<<index_pronunciation[p];
        //std::cout<<std::endl;
        for(unsigned int r = 0; r < pgram.NumRows(); ++r)
        {
          //std::cout<<r<<": ";
          for(unsigned c = 0; c < pgram.NumCols(); ++c)
          {
            //std::cout<<pgram(r,c)<<" ";
            pgram(r,c) = std::log(pgram(r,c));
          }
          //std::cout<<std::endl;
        }

        double final_score = 0;
        index_pronunciation = acousticunitdiscovery::FindViterbiPath(
            pgram, transition, min_frames, final_score);
      }
      else
      {
        // Add code to do viterbi path based pronunciation generation.
        std::vector<utilities::Matrix<double> > pgram_set;
        for(unsigned int i = 0; (i < locations.size()) &&
            (i < max_examples); ++i)
        {
          std::string feature_file = 
              datadir + "/" + locations[i].file + "." + suffix;
          fileutilities::SpeechFeatures sf;
          sf.ReadHtkFile(feature_file);
          utilities::Matrix<double> data = 
              sf.frames(0, locations[i].start, locations[i].end);
          data.Transpose();
          utilities::Matrix<double> pgram = pg.ComputePosteriorgram(data);
          for(unsigned int r = 0; r < pgram.NumRows(); ++r)
            for(unsigned c = 0; c < pgram.NumCols(); ++c)
              pgram(r,c) = std::log(pgram(r,c));
          pgram_set.push_back(pgram);
        }
        index_pronunciation = acousticunitdiscovery::ApproximateViterbiSet(
            pgram_set, transition, min_frames);
      }
      final_pronunciation = ConvertToAlphaPronunciation(index_pronunciation);
      for(unsigned int i = 0; i < final_pronunciation.size(); ++i)
        std::cout<<" "<<final_pronunciation[i];
      std::cout<<std::endl;
    } // end info_line.length() > 0
  } // end wordlist_fin.good()

  wordlist_fin.close();

  return 0;
}
