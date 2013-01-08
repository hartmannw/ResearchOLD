// William Hartmann (hartmannw@gmail.com)
// This is free and unencumbered software released into the public domain.
// See the UNLICENSE file for more information.

#ifndef FILEUTILITIES_SPEECHFEATURES_H_
#define FILEUTILITIES_SPEECHFEATURES_H_

#include<iostream> //just for testing purposes, REMOVE LATER
#include<vector>
#include<string>
#include<fstream>
#include<algorithm>
#include "stdint.h"
#include "sys/stat.h"

#include "Matrix.h"

namespace fileutilities
{
//Online documentation states vector size and magic number are switched, but
//this is NOT true.
typedef struct 
{
  int32_t number_of_records; //number of records
  int16_t vector_size;   //number of features per vector
  int16_t magic_number;  //0 for regular files,
                         //9997 and 9998 for compressed files
} CepFileHeader;

//Standard 12 byte HTK Header
typedef struct
{
  int32_t number_of_samples; //number of samples
  int32_t sample_period; // in 100ns units
  int16_t sample_size; // number of bytes in the entire feature vector for a 
                       // frame of data.  We assume each individual sample is a
                       // 4 byte float, so this value will be # of features per
                       // frame times 4.  If the file is compressed, each 
                       // sample becomes a 2 byte integer, but we do not handle
                       // compressed htk files.
  int16_t parameter_kind; //can mostly be ignored since we generally generate
                          //data outside of HTK.  In that case the value should
                          //always be 9 for user-defined sample kind (USER).
} HtkFileHeader;


class SpeechFeatures
{
 public:
  SpeechFeatures() { feature_width_ = 0; }
  ~SpeechFeatures() {}; //nothing to clean up
  int num_records() { return features_.size(); }
  int feature_width() { return feature_width_; }
  int num_frames(int record) { return features_[record].NumRows(); }
  int feature(int record, int frame, int feature) {
      return features_[record](frame,feature); }
  std::vector<double> frame(int record, int frame) { 
      return features_[record].GetRow(frame); }
  utilities::Matrix<double> frames(int record, int start, int end);
  utilities::Matrix<double> record(int record){ 
      return features_[record]; }

  void Initialize(std::vector<std::vector<double> > record);
  void Initialize(utilities::Matrix<double> record);

  bool ReadCepFile(std::string filename);
  bool ReadHtkFile(std::string filename);
  bool WritePfileAscii(std::string filename);
  bool WriteHtkFile(std::string filename);
  std::vector<double> CalculateRecordMean(int record);

 private:
  int GetFileSize( std::string );
  template <class T>
  void EndianSwap(T *objp);
  //three layers, record, frame, feature
  std::vector< utilities::Matrix<double> > features_;
  int feature_width_;
};

} //end namespace fileutilities

#endif
