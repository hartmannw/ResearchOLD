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

void SaveResultsFile(std::vector<SegmentInfo > segments, 
    std::vector<std::string> filelist, std::string filename);

std::vector<SegmentInfo> FindSegments(std::string utterance, 
    std::vector<std::string> file_names, std::string directory, 
    std::string suffix);

int main(int argc, char *argv[])
{
  if( argc < 4)
  {
    std::cout<<"Usage: file list, data directory, file suffix\n";
    return 1;
  }
  std::vector<std::string> file_names;
  CreateFileNameIndex(std::string(argv[1]), file_names);
  
  std::string utterance, result_file;
  std::vector<SegmentInfo> segments;
  utterance = "BABEL_BP_101_11422_20111019_143344_D2_scripted";
  result_file = std::string(argv[2]) + "/segments/BABEL_BP_101_11422_20111019_143344_D2_scripted.list";
  std::cout<<utterance<<std::endl;

  segments = FindSegments(utterance, file_names, std::string(argv[2]), 
      std::string(argv[3]));

  SaveResultsFile(segments, file_names, result_file);

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

std::vector<SegmentInfo> FindSegments(std::string utterance, 
    std::vector<std::string> file_names, std::string directory, 
    std::string suffix)
{
  const unsigned int MAX_PATHS = 1000;
  fileutilities::SpeechFeatures sf1;
  std::priority_queue<SegmentInfo, std::vector<SegmentInfo>, 
      SegmentInfoComparison> best_paths (SegmentInfoComparison(false));
  std::string utterance_one = directory + "/plp/" + utterance + "." + suffix;
  sf1.ReadHtkFile(utterance_one);
  std::string utterance_one_sil = directory + "/silence/" + utterance + ".sil";
  std::vector<bool> silence = ReadSilence(utterance_one_sil);
  for(int i = 0; i < 5000; i++)
  {
    if(utterance != file_names[i])
    {
    std::string utterance_two, silence_file;
    utterance_two = directory + "/plp/" + file_names[i] + "." + suffix;
    std::cout<<i<<" "<<utterance_two<<std::endl;
    acousticunitdiscovery::DynamicTimeWarp dtw;                           
    fileutilities::SpeechFeatures sf2;                                    
    sf2.ReadHtkFile(utterance_two);
    dtw.set_utterance_one(sf1.record(0));
    dtw.set_utterance_two(sf2.record(0));
    dtw.ComputeSimilarityMatrix();
    dtw.ComputeSegmentalDTW(25);
    dtw.IncreaseSilenceCost(silence);
    dtw.PrunePathsByLCMA(50, 0.1);
    std::vector<acousticunitdiscovery::DtwPath> paths = dtw.paths();

    for(unsigned int s = 0; s < paths.size(); s++)
    {
      SegmentInfo segment;
      segment.utterance_length = sf2.record(0).size();
      segment.score = paths[s].total_score;
      segment.start = paths[s].path[0].first;
      segment.end = paths[s].path[ paths[s].path.size() - 1].first;
      segment.utterance_id = i;
      best_paths.push(segment);
      //std::cout<<segment.utterance_id<<" "<<segment.score<<std::endl;
      if(best_paths.size() > MAX_PATHS)
        best_paths.pop();
      //std::cout<<best_paths.top().score<<" "<<best_paths.top().utterance_id<<std::endl;
    }
    }

  }
  std::cout<<best_paths.size()<<std::endl;
  std::vector<SegmentInfo> result;
  unsigned int best_paths_size = best_paths.size();
  for( unsigned int i = 0; i < best_paths_size; i++)
  {
    result.push_back(best_paths.top());
    best_paths.pop();
  }
  std::reverse(result.begin(), result.end());
  return result;
}

void SaveResultsFile(std::vector<SegmentInfo > segments, 
    std::vector<std::string> filelist, std::string filename)
{
  if( segments.size() < 1) 
    return;

  std::ofstream fout(filename.c_str(), std::ios::out);

  for(unsigned int i =0; i < segments.size(); i++)
  {
    
    fout<<segments[i].start<<" "
        <<segments[i].end<<" "
        <<segments[i].utterance_length<<" "
        <<segments[i].utterance_id<<" "
        <<segments[i].score<<" "
        <<filelist[segments[i].utterance_id]<<std::endl;
  }
  fout.close();
  return;
}

