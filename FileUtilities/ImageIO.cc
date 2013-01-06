#include "ImageIO.h"

namespace fileutilities
{

bool WritePlainPGM(const std::vector< std::vector<double> > &image,
    const unsigned int &max, const std::string &filename)
{
  double rangemax= std::numeric_limits<double>::min();
  double rangemin= std::numeric_limits<double>::max(), range=0;
  int value;
  std::vector< std::vector<double> >::const_iterator rit;
  std::vector< double>::const_iterator cit;

  //find maximum value in the entire image for normalization purposes
  for( rit = image.begin(); rit != image.end(); ++rit)
    for( cit = (*rit).begin(); cit != (*rit).end(); ++cit)
    {
      if( *cit > rangemax)
        rangemax = *cit;
      if( *cit < rangemin)
        rangemin = *cit;
    }

  range = rangemax - rangemin;

  std::ofstream fout;
  fout.open(filename.c_str(), std::ios::out);
  fout<<"P2\n";
  fout<<image[0].size()<<" "<<image.size()<<"\n";
  fout<<max<<"\n"; 
  for( rit = image.begin(); rit != image.end(); ++rit)
  {
    value = static_cast<int>(((*(rit->begin()) - rangemin) / rangemax) * max);
    fout<<value;
    for( cit = ((*rit).begin() + 1); cit != (*rit).end(); ++cit)
    {
      value = static_cast<int>(((*cit - rangemin) / rangemax) * max);
      fout<<" "<<value<<std::endl;
    }
    fout<<"\n";
  }
  fout.close();
  
  return true;
}

bool WriteBinaryPGM(const std::vector< std::vector<double> > &image,
    const std::string &filename)
{
  unsigned int max = 255; // Standard binary PGM uses 8 bits per pixel.
  double rangemax= std::numeric_limits<double>::min();
  double rangemin= std::numeric_limits<double>::max(), range=0;
  unsigned char value;
  std::vector< std::vector<double> >::const_iterator rit;
  std::vector< double>::const_iterator cit;

  //find maximum value in the entire image for normalization purposes
  for( rit = image.begin(); rit != image.end(); ++rit)
    for( cit = (*rit).begin(); cit != (*rit).end(); ++cit)
    {
      if( *cit > rangemax)
        rangemax = *cit;
      if( *cit < rangemin)
        rangemin = *cit;
    }

  range = rangemax - rangemin;

  std::ofstream fout;
  fout.open(filename.c_str(), std::ios::out);
  fout<<"P5\n";
  fout<<image[0].size()<<" "<<image.size()<<"\n";
  fout<<max<<"\n";
  fout.close();
  fout.open(filename.c_str(), std::ios::app|std::ios::binary);
  for( rit = image.begin(); rit != image.end(); ++rit)
  {
    for( cit = rit->begin(); cit != rit->end(); ++cit)
    {
      value = static_cast<unsigned char>(((*cit - rangemin) / rangemax) * max);
      fout<<value;
    }
  }
  fout.close();
  
  return true;
}

}// end namespace fileutilities
