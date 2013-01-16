#include<iostream>
#include<string>
#include<vector>
#include "PosteriorgramGenerator.h"
#include "HmmSet.h"
#include "Matrix.h"
#include "ImageIO.h"
#include "SpeechFeatures.h"

int main()
{

  //Testing Hmms from HTK                                                        
  statistics::HmmSet htk;                                                        
  //std::string filename("/people/hartmann/research/SegmentalModel/tidigit_exp/hmm_plp_phone/hmm29/hmmdefs");
  std::string filename("/people/hartmann/research/AASP_CASA/hmm_training/hmm_mfccmvn_cv1/hmm24/hmmdefs");
  htk.LoadHtkHmmSet(filename);                                                   
  
  std::vector<statistics::MixtureOfDiagonalGaussians> mog = htk.states();        
  std::vector<std::vector<std::string> > names = htk.mixture_names();
  utilities::Matrix<double> similarity_matrix, pgram;
  statistics::PosteriorgramGenerator pg;
  pg.SetGaussians(mog);
  similarity_matrix = pg.ComputeSimilarityMatrix();
  fileutilities::WriteBinaryPGM(similarity_matrix.GetVectorOfVectors(), 
      std::string("simmx.pgm") );

  filename = "/people/hartmann/research/AASP_CASA/data/SceneClassification/scenes_stereo/bus01.mfccmvn";
  fileutilities::SpeechFeatures sf;
  sf.ReadHtkFile(filename);
  utilities::Matrix<double> data = sf.record(0);
  data.Transpose();
  pgram = pg.ComputePosteriorgram(data);
  std::vector<int> one_best = pg.BestIndexPerFrame(pgram);
  //for( unsigned int i = 0; i < one_best.size(); ++i)
  //  std::cout<<one_best[i]<<" ";
  //std::cout<<std::endl;
  
  for(unsigned int i = 0; i < names.size(); ++i)
    std::cout<<i<<" "<<names[i][0]<<std::endl;

  std::vector< std::pair<int, double> > posterior_count = pg.BestIndexMass(pgram);
  for(unsigned int i = 0; i < posterior_count.size(); ++i)
    std::cout<<posterior_count[i].first<<" "<<posterior_count[i].second<<
        std::endl;

  fileutilities::WriteBinaryPGM(pgram.GetVectorOfVectors(), 
      std::string("bus01.pgm") );

  // Lets examine this similarity matrix
  for(unsigned int r = 0; r < similarity_matrix.NumRows(); ++r)
  {
    double maxval = 0;
    int maxind = 1;
    for(unsigned int c = 0; c < similarity_matrix.NumCols(); ++c)
    {
      if( r != c && similarity_matrix(r,c) > maxval)
      {
        maxval = similarity_matrix(r,c);
        maxind = c;
      }
    }
    //std::cout<<r<<" "<<maxind<<std::endl;
  }
  return 0;
}
