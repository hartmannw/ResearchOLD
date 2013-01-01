// William Hartmann (hartmannw@gmail.com)

#include<vector>
#include<string>
#include<iostream>
#include<fstream>
#include<sstream>
#include "DynamicTimeWarp.h"
#include "SpeechFeatures.h"

int stoi(std::string s)
{                                                                               
  int ret;                                                                      
  std::stringstream ss(s);                                                      
  ss >> ret;                                                                    
  return ret;                                                                   
}

double stof(std::string s)
{
  double ret;
  std::stringstream ss(s);
  ss >> ret;
  return ret;
}

bool TokenizeString(std::string line, char delimiter, 
    std::vector<std::string> &tokens);

bool CreateFileNameIndex(std::string filelist, 
    std::vector<std::string> &file_names);

std::vector<bool> ReadSilence(std::string filename);
std::vector<double> ReadFrameScores(std::string filename);

void SaveResultsFile(std::vector<std::vector<double> > scores,
    std::string filename);

bool SegmentAudio(const std::vector<std::vector<double> > &features, 
    const std::vector<double> silence, double threshold, int minimum_length, 
    std::vector<std::vector<std::vector<double> > > &segments, 
    std::vector<int> &offsets);

int main(int argc, char *argv[])
{
  if( argc < 5)
  {
    std::cout<<"Usage: file list, word reference file, data directory,"<<
        " file suffix\n";
    return 1;
  }
  std::vector<std::string> file_names;
  CreateFileNameIndex(std::string(argv[1]), file_names);
  std::ifstream fin(argv[2], std::ios::in);
  std::string line;
  // for testing purposes, skip some lines
  for(int i = 0; i < 358; i++) // previously 358
    getline(fin, line);

  int word = 0;
  while( getline(fin, line) )
  {
    ++word;
    std::vector<std::string> tokens;
    TokenizeString(line, ' ', tokens);
    std::cout<<line<<std::endl;
    for(int i = 1; i < tokens.size(); ++i)
    {
      std::vector<std::string> utterance_count;
      TokenizeString(tokens[i], '.', utterance_count);
      int utterance_one_id = stoi(utterance_count[0]);
      std::string utterance_one, utterance_one_sil;
      utterance_one = std::string(argv[3]) + "/mog_pgram/" + 
          file_names[utterance_one_id] + argv[4];
      utterance_one_sil = std::string(argv[3]) + "/word_locations/" + 
          file_names[utterance_one_id] + ".<silence>";
      fileutilities::SpeechFeatures sf1;
      std::vector<double> silence;
      sf1.ReadHtkFile(utterance_one);
      silence = ReadFrameScores(utterance_one_sil);
      std::vector<std::vector<std::vector<double> > > segments;
      std::vector<int> offsets;
      SegmentAudio(sf1.record(0), silence, 0.99, 25, segments, offsets);
      std::cout<<"Word: "<<word<<" example "<<i<<" of "<<(tokens.size()-1)<<
          " "<<utterance_one<<" "<<sf1.record(0).size()<<" "<<
          segments.size()<<std::endl;
      std::vector<std::vector< double> > best_scores;
      //probably don't really need to evaluate every pair at this point
      for(int j = 1; j < std::min(50, static_cast<int>(tokens.size())); j++)
      {
        if( i != j)
        {
          for(unsigned int s=0; s < segments.size(); ++s)
          {
            std::vector<std::string> utterance_count;
            TokenizeString(tokens[j], '.', utterance_count);
            int utterance_two_id = stoi(utterance_count[0]);
            acousticunitdiscovery::DynamicTimeWarp dtw;
            fileutilities::SpeechFeatures sf2;
            std::string utterance_two;
            utterance_two = std::string(argv[3]) + "/mog_pgram/" + 
               file_names[utterance_two_id] + argv[4];
            sf2.ReadHtkFile(utterance_two);
            //std::cout<<"setting utterance"<<std::endl;
            dtw.set_utterance_one(segments[s]);
            dtw.set_utterance_two(sf2.record(0));
            //std::cout<<"computing matrix"<<std::endl;
            dtw.ComputeSimilarityMatrix();
            //std::cout<<"seg dtw "<<segments[s].size()<<std::endl;
            dtw.ComputeSegmentalDTW(20);
            //std::cout<<"prune"<<std::endl;
            dtw.PrunePathsByLCMA(25, 0.1);
            //std::cout<<"push back"<<std::endl;
            std::vector<double> segment_score = dtw.BestScorePerFrame();
            std::vector<double> score(silence.size(), -1);
            //std::cout<<offsets[s]<<std::endl;
            for(unsigned int  k = 0; k < segment_score.size(); ++k)
              score[k + offsets[s]] = segment_score[k];
            best_scores.push_back(score);
            //std::cout<<"pushed"<<std::endl;
          }
        }
      }
      std::string result_file = std::string(argv[3]) + "/word_locations2/" +
          file_names[utterance_one_id] + "." + tokens[0];
      //std::cout<<result_file<<std::endl;
      SaveResultsFile(best_scores, result_file);
      //std::cout<<"results saved"<<std::endl;
      //fileutilities::WriteBinaryPGM(best_scores, result_file);

    }
  }
  return 0;
}

bool SegmentAudio(const std::vector<std::vector<double> > &features, 
    const std::vector<double> silence, double threshold, int minimum_length, 
    std::vector<std::vector<std::vector<double> > > &segments, 
    std::vector<int> &offsets)
{
  bool segment_started = false;
  int start = 0;

  for(unsigned int i = 0; i < silence.size(); ++i)
  {
    if(silence[i] >= threshold && segment_started)
    {
      segment_started = false;
      if( (i - start + 1) >= minimum_length)
      {
        segments.push_back( std::vector<std::vector<double> >
            (features.begin()+start, features.begin()+i));
        offsets.push_back(start);
      }
    }
    else if(silence[i] < threshold && !segment_started)
    {
      segment_started = true;
      start = i;
    }
  }
  if( segment_started && (silence.size() - start) >= minimum_length)
  {
    segments.push_back( std::vector<std::vector<double> > 
        (features.begin()+start, features.end()-1));
    offsets.push_back(start);
  }
    
  return true;
}
bool TokenizeString(std::string line, char delimiter, 
    std::vector<std::string> &tokens)
{
  std::istringstream iss(line);
  std::string token;
  while(getline(iss, token, delimiter))
    if( token.length() > 0)
      tokens.push_back(token);

  return true;
}

bool CreateFileNameIndex(std::string filelist, 
    std::vector<std::string> &file_names)
{
  std::ifstream fin(filelist.c_str());
  std::string line;

  if(!fin.is_open())
    return false;

  while( fin.good() )
  {
    std::getline(fin, line);
    //std::cout<<line;
    file_names.push_back(line);
  }
  fin.close();
  return true;

}

std::vector<double> ReadFrameScores(std::string filename)
{
  std::vector<double> ret;
  std::ifstream fin(filename.c_str());
  std::string line;

  if(!fin.is_open())
    return ret;
  while( fin.good() )
  {
    std::getline(fin, line);
    if(line.length() > 0)
      ret.push_back( stof(line) );
  }
  return ret;
}

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

void SaveResultsFile(std::vector<std::vector<double> > scores,
    std::string filename)
{
  if( scores.size() < 1) 
    return;

  std::ofstream fout(filename.c_str(), std::ios::out);
  std::vector<double> result;
  result.resize( scores[0].size(), 0);
  for(unsigned int i = 0; i < scores.size(); i++)
    for(unsigned int j = 0; j < scores[i].size(); j++)
      result[j] += (scores[i][j] / scores.size());

  // Find the maximum value.
  double max_value=0;
  for(unsigned int i =0; i < result.size(); i++)
    if(result[i] > max_value)
      max_value = result[i];

  // Normalize and flip the result.
  for(unsigned int i =0; i < result.size(); i++)
  {
    result[i] = 1 - (result[i] / max_value);
    if(result[i] > 1)
      result[i] = 0;
    fout<<result[i]<<std::endl;
  }
  fout.close();
  return;
}

