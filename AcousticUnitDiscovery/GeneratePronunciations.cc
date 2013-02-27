// William Hartmann (hartmannw@gmail.com)
// This is free and unencumbered software released into the public domain.
// See the UNLICENSE file for more information.

// Used for generating pronunciations for a list of words.

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


// Stores all the basic information that is needed for generating 
// pronunciations. Most of the settings can be set by the command line. By
// keeping them together in a structure, parameter lists are shorter and code is
// cleaner.
typedef struct
{
  std::string word_information;
  std::string locationdir;
  std::string hmmfile;
  std::string cluster_file;
  std::string datadir;
  std::string suffix;
  unsigned int max_examples;
  unsigned int min_examples;
  double length_variance;
  double self_transition;
  unsigned int min_frames;
  unsigned int max_frames_per_hmm;
  bool modify_transition;
  unsigned int pronunciation_type;
  unsigned int attempts_per_example; // Maximum number of attempts to generate
                                     // a valid HMM sample for each example.
  double min_hmm_transition; // Minimum transition value before we assume it
                             // will exist in the state for at least one frame.

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

std::vector<utilities::Matrix<double> > LoadPosteriorgramData(
    const std::vector<WordLocation> &locations,
    const statistics::PosteriorgramGenerator &pg,
    const PronunciationParameters &param)
{
  std::vector<utilities::Matrix<double> > pgram_set;
  for(unsigned int i = 0; (i < locations.size()) &&
      (i < param.max_examples); ++i)
  {
    std::string feature_file = 
        param.datadir + "/" + locations[i].file + "." + param.suffix;
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
  return pgram_set;
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

double ExpectedWordLength(
    std::vector<std::string> &triphone_pronunciation, statistics::HmmSet &htk)
{
  double ret = 0;
  for(unsigned int p = 0; p < triphone_pronunciation.size(); ++p)
    ret += htk.Hmm(triphone_pronunciation[p]).ExpectedDuration();
  return ret;
}

std::vector<std::vector<double> > GenerateSampleData( 
    std::vector<std::string> &triphone_pronunciation, statistics::HmmSet &htk,
    std::vector<statistics::MixtureOfDiagonalGaussians> &mog,
    std::default_random_engine &generator, const PronunciationParameters &param)
{
  std::vector<std::vector<double> >ret;
  std::uniform_real_distribution<double> distribution(0.0,1.0);
  for(unsigned int p = 0; p < triphone_pronunciation.size(); ++p)
  {
    statistics::HiddenMarkovModel hmm = htk.Hmm(triphone_pronunciation[p]);
    unsigned int h = 0;
    unsigned int count = 0;
    while(h < hmm.NumberOfStates() && count < param.max_frames_per_hmm)
    {
      count++;
      double trans = distribution(generator);

      if(hmm.transition(h+1, h+1) > param.min_hmm_transition)
      {
        // Add frame only the first loop iteration because you must spend at
        // least one frame in each state.
        std::vector<double> frame = mog[hmm.state(h)].Sample(generator);
        ret.push_back(frame);
        if( trans > hmm.transition(h+1, h+1) )
          ++h;
      }
      else
        ++h;
    }
  }
  return ret;
}

void AppendSampleData(std::vector<std::string> &triphone_pronunciation, 
    statistics::HmmSet &htk, 
    std::vector<statistics::MixtureOfDiagonalGaussians> &mog,
    const statistics::PosteriorgramGenerator &pg,
    std::default_random_engine &generator, const PronunciationParameters &param,
    double mean_length, std::vector<utilities::Matrix<double> > &pgram_set)
{
  unsigned int attempts = 0;
  double min_length = (1 - param.length_variance) * mean_length;
  double max_length = (1 + param.length_variance) * mean_length;
  while((pgram_set.size() < param.min_examples) && 
        (attempts < (param.min_examples * param.attempts_per_example)) )
  {
    attempts++;
    std::vector<std::vector<double> > data;
    data = GenerateSampleData(triphone_pronunciation, htk, mog, 
        generator, param);
    unsigned int frames = data.size();
    if( frames > min_length && frames < max_length)
    {
      utilities::Matrix<double> sample, pgram;
      sample.Initialize(data);
      sample.Transpose();
      pgram = pg.ComputePosteriorgram(sample);
      for(unsigned int r = 0; r < pgram.NumRows(); ++r)
        for(unsigned c = 0; c < pgram.NumCols(); ++c)
          pgram(r,c) = std::log(pgram(r,c));
      pgram_set.push_back(pgram);
    }
  }
}

std::vector<int> GenerateClusteredPronunciation(
    const std::vector<std::string> &triphone_pronunciation,
    const statistics::HmmSet &htk, const std::vector<int> &cluster_index, 
    const PronunciationParameters &param)
{
  std::vector<int> ret;
  for(unsigned int p = 0; p < triphone_pronunciation.size(); ++p)
  {
    statistics::HiddenMarkovModel hmm = htk.Hmm(triphone_pronunciation[p]);
    //unsigned int best_index = 0;
    std::vector<double> expected_length(param.total_clusters, 0);
    for(unsigned int h = 0; h < hmm.NumberOfStates(); ++h)
    {
      expected_length[ cluster_index[ hmm.state(h) ] ] +=
          hmm.transition(h+1, h+1);
      //if(hmm.transition(h+1, h+1) > hmm.transition(best_index+1, best_index+1))
      //  best_index = h; 
    }
    unsigned best_index = 0;
    for(unsigned int i = 0; i < expected_length.size(); ++i)
      if(expected_length[i] > expected_length[best_index])
        best_index = i;
    if(expected_length[best_index] > param.min_hmm_transition)
      ret.push_back( best_index);
  }
  return ret;
}

std::vector<int> GenerateMogPronunciation(
    const std::vector<std::string> &triphone_pronunciation,
    const statistics::HmmSet &htk, 
    const std::vector<statistics::MixtureOfDiagonalGaussians> &mog,
    std::default_random_engine &generator, 
    const utilities::Matrix<double> &transition,
    const statistics::PosteriorgramGenerator &pg,
    const std::vector<int> &cluster_index,
    const PronunciationParameters &param)
{
  std::vector<int> ret;
  utilities::Matrix<double> data, modified_transition, transition_count;
  
  if( param.modify_transition )
  {
    modified_transition = transition;
    for(unsigned int r = 0; r < modified_transition.NumRows(); ++r)
      for(unsigned int c = 0; c < modified_transition.NumCols(); ++c)
        modified_transition(r,c) = std::exp(modified_transition(r,c));
    transition_count.Initialize(param.total_clusters, 
        param.total_clusters, 1);
  }
        
  data.Initialize(param.dimension, triphone_pronunciation.size() * 3);
  for(unsigned int p = 0; p < triphone_pronunciation.size(); ++p)
  {
    statistics::HiddenMarkovModel hmm = htk.Hmm(triphone_pronunciation[p]);
    for(unsigned int h = 0; h < hmm.NumberOfStates(); ++h)
    {
      std::vector<double> sample = mog[hmm.state(h)].Sample(generator);
      for(unsigned int m = 0; m < sample.size(); ++m)
        data(m, (p*3)+h) = sample[m];
    }
    if( param.modify_transition )
    {
      for(unsigned int h = 0; h < hmm.NumberOfStates() - 1; ++h)
      {
        unsigned int from_index = cluster_index[ hmm.state(h) ];
        unsigned int to_index = cluster_index[ hmm.state(h+1) ];
        modified_transition(from_index, from_index) += hmm.transition(h+1, h+1);
        transition_count(from_index, from_index)++;
        modified_transition(from_index, to_index) += hmm.transition(h+1, h+2);
        transition_count(from_index, to_index)++;
      }
      unsigned int last_hmm = hmm.NumberOfStates() - 1;
      unsigned int from_index = cluster_index[ hmm.state(last_hmm) ];
      modified_transition(from_index, from_index) += 
          hmm.transition(last_hmm+1, last_hmm+1);
      transition_count(from_index, from_index)++;
      for(unsigned int to_index = 0; to_index < param.total_clusters; 
          ++to_index)
      {
        if( to_index != from_index )
        {
          modified_transition(from_index, to_index) += 
              hmm.transition(last_hmm+1, last_hmm+2);
          transition_count(from_index, to_index)++;
        }
      }
    } // End if param.modify_transition
  } // End for triphone_pronunciation.size()

  if( param.modify_transition )
  {
    for(unsigned int r = 0; r < modified_transition.NumRows(); ++r)
      for(unsigned int c = 0; c < modified_transition.NumCols(); ++c)
        modified_transition(r,c) = std::log(modified_transition(r,c) / 
            transition_count(r,c));
  }
  
  utilities::Matrix<double> pgram = pg.ComputePosteriorgram(data);
  for(unsigned int r = 0; r < pgram.NumRows(); ++r)
    for(unsigned c = 0; c < pgram.NumCols(); ++c)
      pgram(r,c) = std::log(pgram(r,c));

  double final_score = 0;
  if( param.modify_transition )
    return acousticunitdiscovery::FindViterbiPath(
        pgram, modified_transition, param.min_frames, final_score);

  return acousticunitdiscovery::FindViterbiPath(
      pgram, transition, param.min_frames, final_score);
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
  PronunciationParameters param;
  param.word_information = std::string(argv[1]);
  param.locationdir = std::string(argv[2]);
  param.hmmfile = std::string(argv[3]);
  param.cluster_file = std::string(argv[4]);
  param.datadir = std::string(argv[5]);
  param.suffix = std::string(argv[6]);
  param.max_examples = utilities::ToNumber<unsigned int>(std::string(argv[7]));
  param.min_examples = utilities::ToNumber<unsigned int>(std::string(argv[8]));
  param.length_variance = utilities::ToNumber<double>(std::string(argv[9]));
  param.self_transition = utilities::ToNumber<double>(std::string(argv[10]));
  param.min_frames = 3;
  param.max_frames_per_hmm = 30;
  param.modify_transition = false;
  param.pronunciation_type = 0; // 0=clustered, 1=MOG, 2=HMM !Change to Enum!
  param.attempts_per_example = 10;
  param.min_hmm_transition = 0.1;
  std::default_random_engine generator; // A single RNG.

  statistics::HmmSet htk;
  htk.LoadHtkHmmSet(param.hmmfile);
  std::vector<statistics::MixtureOfDiagonalGaussians> mog = htk.states();
  // Ugly way of initializing dimension, Find a better way!!!
  param.dimension = mog[0].gaussian(0).dimension();
  statistics::PosteriorgramGenerator pg;
  utilities::Matrix<double> transition;
  std::vector<int> cluster_index = ReadVectorFromFile(param.cluster_file);
  param.total_clusters = cluster_index.size();
  transition = acousticunitdiscovery::GenerateTransitionMatrix(
      param.total_clusters, param.self_transition);
  pg.SetGaussians(mog, cluster_index);
  
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

      std::vector<std::string> final_pronunciation;
      std::vector<int> index_pronunciation;
      std::vector<std::string> triphone_pronunciation = 
          GenerateTriphonePronunciation(original_pronunciation);

      if(locations.size() < param.min_examples && 
          param.pronunciation_type == 0) // clusterd type
      {
        index_pronunciation = GenerateClusteredPronunciation(
            triphone_pronunciation, htk, cluster_index, param);
      }
      else if(locations.size() < param.min_examples &&
          param.pronunciation_type == 1) // MOG
      {
        index_pronunciation = GenerateMogPronunciation( triphone_pronunciation,
            htk, mog, generator, transition, pg, cluster_index, param);
      }
      else if(locations.size() < param.min_examples &&
          param.pronunciation_type == 2) // HMM
      {
        std::vector<utilities::Matrix<double> > pgram_set = 
            LoadPosteriorgramData(locations, pg, param);

        double mean = 0;
        if(pgram_set.size() == 0) // Set mean based on HMM.
        {
          mean = ExpectedWordLength(triphone_pronunciation, htk);
        }
        else // Set mean based on the length of examples
        {
          for(unsigned int i = 0; i < locations.size(); ++i)
            mean += (locations[i].end - locations[i].start + 1);
          mean = mean / locations.size();
        }
        AppendSampleData(triphone_pronunciation, htk, mog, pg,generator, param, 
            mean, pgram_set);
        index_pronunciation = acousticunitdiscovery::BestPathInSet(
            pgram_set, transition, param.min_frames);

      }
      else // Generate data solely from data.
      {
        std::vector<utilities::Matrix<double> > pgram_set = 
            LoadPosteriorgramData(locations, pg, param);
        //index_pronunciation = acousticunitdiscovery::BestPathInSet(
        //    pgram_set, transition, param.min_frames);
        index_pronunciation = acousticunitdiscovery::ApproximateViterbiSet2(
            pgram_set, transition, param.min_frames, 
            original_pronunciation.size());
      }
      std::cout<<" "<<locations.size();
      final_pronunciation = ConvertToAlphaPronunciation(index_pronunciation);
      for(unsigned int i = 0; i < final_pronunciation.size(); ++i)
        std::cout<<" "<<final_pronunciation[i];
      std::cout<<std::endl;
    } // end info_line.length() > 0
  } // end wordlist_fin.good()

  wordlist_fin.close();

  return 0;
}
