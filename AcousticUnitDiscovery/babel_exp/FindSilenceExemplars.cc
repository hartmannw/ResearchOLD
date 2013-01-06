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

std::vector<SegmentInfo> FindSilence(std::vector<std::string> file_names, 
    std::string directory, std::string suffix, int number_of_comparisons,
    int number_of_segments);

int main(int argc, char *argv[])
{
  if( argc < 6)
  {
    std::cout<<"Usage: file list, data directory, file suffix, number of comparisons, number of examples \n";
    return 1;
  }
  std::vector<std::string> file_names;
  CreateFileNameIndex(std::string(argv[1]), file_names);
  
  std::vector<SegmentInfo> segments;
  std::string result_file = "silence.test";

  segments = FindSilence(file_names, std::string(argv[2]), 
      std::string(argv[3]), stoi(std::string(argv[4])), 
      stoi(std::string(argv[5])));

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

std::vector<SegmentInfo> FindSilence(std::vector<std::string> file_names, 
    std::string directory, std::string suffix, int number_of_comparisons,
    int number_of_segments)
{
  std::priority_queue<SegmentInfo, std::vector<SegmentInfo>, 
      SegmentInfoComparison> best_paths (SegmentInfoComparison(false));

  for(unsigned int i = 0; i < number_of_comparisons; ++i)
  {
    int one = rand() % file_names.size();
    int two = rand() % file_names.size();

    if( one == two)
      continue;

    std::string utterance_one, utterance_two;
    utterance_one = directory + "/plp/" + file_names[one] + "." + suffix;
    utterance_two = directory + "/plp/" + file_names[two] + "." + suffix;
    std::cout<<utterance_one<<std::endl;
    fileutilities::SpeechFeatures sf1, sf2;
    sf1.ReadHtkFile(utterance_one);
    sf2.ReadHtkFile(utterance_two);
    acousticunitdiscovery::DynamicTimeWarp dtw;
    dtw.set_utterance_one(sf1.record(0));
    dtw.set_utterance_two(sf2.record(0));
    dtw.ComputeSimilarityMatrix();
    dtw.ComputeSegmentalDTW(50);
    dtw.PrunePathsByLCMA(100, 0.1);
    std::vector<acousticunitdiscovery::DtwPath> paths = dtw.paths();

    for(unsigned int s = 0; s < paths.size(); s++)
    {
      SegmentInfo segment;
      segment.utterance_length = sf1.record(0).size();
      segment.score = paths[s].total_score;
      segment.start = paths[s].path[0].first;
      segment.end = paths[s].path[ paths[s].path.size() - 1].first;
      segment.utterance_id = one;
      best_paths.push(segment);
      if(best_paths.size() > number_of_segments)
        best_paths.pop();
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
    
    fout<<segments[i].utterance_id<<" "
        <<segments[i].start<<" "
        <<segments[i].end<<" "
        <<segments[i].score<<" "
        <<segments[i].utterance_length<<" "
        <<filelist[segments[i].utterance_id]<<std::endl;
  }
  fout.close();
  return;
}

