#include "SpeechFeatures.h"

namespace fileutilities
{

bool SpeechFeatures::ReadCepFile( std::string filename )
{
  std::ifstream fin;
  fin.open(filename.c_str(), std::ios::in|std::ios::binary);
  if( !fin.is_open() )
    return false;

  CepFileHeader header;
  fin.read (reinterpret_cast<char*>(&header), sizeof(header));
    
  int32_t *number_of_frames;
  number_of_frames = new int32_t [header.number_of_records];
  fin.read(reinterpret_cast<char*>(number_of_frames),
      sizeof(number_of_frames) * header.number_of_records);


  feature_width_ = header.vector_size;
  features_.resize(header.number_of_records);
  for(int r = 0; r < header.number_of_records; ++r)
  {
    features_[r].Initialize(number_of_frames[r], feature_width_);
    float *vector;
    vector = new float [feature_width_];
    for(int f = 0; f < number_of_frames[r]; ++f)
    {
      fin.read(reinterpret_cast<char*> (vector),                                
          sizeof(float) * feature_width_);
      for(int c = 0; c < feature_width_; ++c)
      {
        features_[r](f,c) = static_cast<double>(vector[c]);
      }
    }
    delete[] vector;
  }

  delete[] number_of_frames;

  fin.close();
  return true;
}

bool SpeechFeatures::ReadHtkFile( std::string filename )
{
  std::ifstream fin;
  //int filesize = GetFileSize( filename );
  fin.open(filename.c_str(), std::ios::in|std::ios::binary);
  if( !fin.is_open() )
    return false;

  HtkFileHeader header; 
  fin.read (reinterpret_cast<char*>(&header), sizeof(header));
  EndianSwap(&header.number_of_samples);
  EndianSwap(&header.sample_size);
  EndianSwap(&header.sample_period);
  EndianSwap(&header.parameter_kind);
  feature_width_ = header.sample_size / 4; //we are assuming each individual 
                                           //feature is a 4 byte float
  features_.resize(1); //htk files only store one record at a time.
  features_[0].Initialize(header.number_of_samples, feature_width_);
  float *vector;
  vector = new float [feature_width_];
  for(int f = 0; f < header.number_of_samples; ++f)
  {
    fin.read(reinterpret_cast<char*> (vector), sizeof(float)*feature_width_);
    for(int c = 0; c < feature_width_; ++c)
    {
      EndianSwap(&vector[c]);
      features_[0](f,c) = static_cast<double>(vector[c]);
    }
  }
  delete[] vector;
  return true;
}

bool SpeechFeatures::WritePfileAscii(std::string filename)
{
  std::ofstream fout;
  fout.open(filename.c_str(), std::ios::out);

  for(unsigned int r = 0; r < features_.size(); ++r)
    for(unsigned int f = 0; f < features_[r].NumRows(); ++f)
    {
      fout << r <<" "<< f <<" "<< features_[r](f,0);
      for(unsigned int c = 1; c < features_[r].NumCols(); ++c)
        fout << " " << features_[r](f,c);
      fout << "\n";
    }

  fout.close();
  return true;
}

bool SpeechFeatures::WriteHtkFile(std::string filename)
{
  std::ofstream fout;
  fout.open(filename.c_str(), std::ios::out|std::ios::binary);
  
  HtkFileHeader header;
  header.number_of_samples = features_[0].NumRows();
  header.sample_period = 0;
  header.sample_size = features_[0].NumCols() * 4;
  header.parameter_kind = 9;
  EndianSwap(&header.number_of_samples);
  EndianSwap(&header.sample_size);
  EndianSwap(&header.sample_period);
  EndianSwap(&header.parameter_kind);

  fout.write( reinterpret_cast<char*>(&header), sizeof(header));

  for(unsigned int f = 0; f < features_[0].NumRows(); f++)
    for(unsigned int c = 0; c < features_[0].NumCols(); c++)
    {
      float item = static_cast<float>(features_[0](f,c));
      EndianSwap(&item);
      fout.write( reinterpret_cast<char*>(&item), sizeof(item));
    }

  fout.close();
  return true;
}

void SpeechFeatures::Initialize(std::vector<std::vector<double> > record)
{
  feature_width_ = record[0].size();
  features_.clear();
  features_.resize(1);
  features_[0].Initialize(record);
}

void SpeechFeatures::Initialize( utilities::Matrix<double> record)
{
  feature_width_ = record.NumCols();
  features_.clear();
  features_.push_back(record);
}

int SpeechFeatures::GetFileSize( std::string filename)
{
  struct stat stat_buffer;
  stat(filename.c_str(), &stat_buffer);
  return static_cast<int>(stat_buffer.st_size);
}

template <class T>
void SpeechFeatures::EndianSwap(T *objp)
{
  unsigned char *memp = reinterpret_cast<unsigned char*>(objp);
  std::reverse(memp, memp + sizeof(T));
  return;
}


std::vector<double> SpeechFeatures::CalculateRecordMean(int record)
{
  std::vector<double> mean;
  mean.resize( features_[record].NumCols() );

  for(unsigned int i = 0; i < mean.size(); ++i)
    mean[i] = 0;

  for(unsigned int f = 0; f < features_[record].NumRows(); ++f)
    for( unsigned int c = 0; c < mean.size(); ++c)
      mean[c]+= features_[record](f,c);

  for(unsigned int i = 0; i < mean.size(); ++i)
    mean[i] = mean[i] / features_[record].NumRows();

  return mean;
}

utilities::Matrix<double> SpeechFeatures::frames(int record, int start,
    int end)
{
  utilities::Matrix<double> ret;
  int numRows = (end-start) + 1;
  ret.Initialize(numRows, feature_width_);
  for(int f = 0; f < numRows; ++f)
    for(int c = 0; c < feature_width_; ++c)
      ret(f, c) = features_[record](f+start, c);
  return ret;
}

} // end namespace fileutilities
