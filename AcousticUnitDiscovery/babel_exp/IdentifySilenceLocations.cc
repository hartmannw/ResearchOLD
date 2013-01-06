// William Hartmann (hartmannw@gmail.com)

#include<vector>
#include<queue>
#include<string>
#include<iostream>
#include<fstream>
#include<sstream>
#include<algorithm>
#include "DynamicTimeWarp.h"
#include "SpeechFeatures.h"

typedef struct
{
  unsigned int start;
  unsigned int end;
  double score;
  unsigned int utterance_length;
  unsigned int utterance_id;
} SegmentInfo;

class SegmentInfoComparison
{
 private:
  bool reverse;
 public:
  SegmentInfoComparison(const bool& revparam=false) {reverse=revparam;}
  bool operator() (const SegmentInfo& lhs, const SegmentInfo& rhs) const
  {
    if (reverse) return (lhs.score>rhs.score);
    else return (lhs.score<rhs.score);
  }
};

int stoi(std::string s)                                                              
{                                                                               
  int ret;                                                                      
  std::stringstream ss(s);                                                      
  ss >> ret;                                                                    
  return ret;                                                                   
}

bool TokenizeString(std::string line, char delimiter, 
    std::vector<std::string> &tokens);

bool CreateFileNameIndex(std::string filelist, 
    std::vector<std::string> &file_names);

std::vector<bool> ReadSilence(std::string filename);

void SaveResultsFile(std::vector<std::vector<double> > scores,                  
    std::string filename);

std::vector<SegmentInfo> FindSilence(std::vector<std::string> file_names, 
    std::string directory, std::string suffix, int number_of_comparisons,
    int number_of_segments);

std::vector< std::vector<std::vector<double> > > 
    LoadSilence(std::string directory, std::string suffix, 
    std::string silence_file);

int main(int argc, char *argv[])
{
  if( argc < 5)
  {
    std::cout<<"Usage: file list, data directory, file suffix, silence file\n";
    return 1;
  }
  std::vector<std::string> file_names;
  CreateFileNameIndex(std::string(argv[1]), file_names);
  
  std::vector<SegmentInfo> segments;
  std::string silence_file = std::string(argv[4]);
  std::vector< std::vector<std::vector<double> > > silence_exemplars;
  silence_exemplars = LoadSilence(std::string(argv[2]), std::string(argv[3]),
      silence_file);

  for(unsigned int i = 0; i < file_names.size(); ++i)
  {
    std::string utterance = std::string(argv[2]) + "/plp/" + file_names[i] +
        "." + std::string(argv[3]);
    std::cout<<(i+1)<<" "<<utterance<<std::endl;
    fileutilities::SpeechFeatures sf;
    sf.ReadHtkFile(utterance);
    std::vector<std::vector<double> > best_scores;
    for(unsigned j = 0; j < silence_exemplars.size(); ++j)
    {
      acousticunitdiscovery::DynamicTimeWarp dtw;
      dtw.set_utterance_one(sf.record(0));
      dtw.set_utterance_two(silence_exemplars[j]);
      dtw.ComputeSimilarityMatrix();
      dtw.ComputeSegmentalDTW(15);
      dtw.PrunePathsByLCMA(30, 0.0);
      best_scores.push_back(dtw.BestScorePerFrame());
    }
    std::string result_file = std::string(argv[2]) + "/word_locations/" + 
        file_names[i] + ".<silence>";
    SaveResultsFile(best_scores, result_file);
  }

  return 0;
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

std::vector< std::vector<std::vector<double> > > 
    LoadSilence(std::string directory, std::string suffix, 
    std::string silence_file)
{
  std::vector<std::vector<std::vector<double> > > result;
  std::ifstream fin(silence_file.c_str());
  std::string line;

  if(!fin.is_open())
    return result;
  std::getline(fin, line);
  while(line.length() > 5)
  {
    std::vector<std::string> tokens;
    TokenizeString(line, ' ', tokens);
    std::cout<<tokens[1]<<" "<<tokens[2]<<" "<<tokens[5]<<std::endl;

    std::string file_name = directory + "/plp/" + tokens[5] + "." + suffix;
    fileutilities::SpeechFeatures sf;
    sf.ReadHtkFile(file_name);
    std::vector<std::vector<double> > full_vector = sf.record(0);
    result.push_back( std::vector< std::vector<double> > (
        full_vector.begin() + stoi(tokens[1]), 
        full_vector.begin() + stoi(tokens[2])));

    std::getline(fin, line);
  }
  std::cout<<result.size()<<std::endl;
  return result;
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
    fout<<result[i]<<std::endl;                                                 
  }                                                                             
  fout.close();                                                                 
  return;                                                                       
}      
