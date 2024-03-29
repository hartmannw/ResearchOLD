#include "DiagonalGaussian.h"
#include "MixtureOfDiagonalGaussians.h"
#include "HmmSet.h"
#include "SpectralCluster.h"
#include "ImageIO.h"
#include "MatrixFunctions.h"
#include <iostream>
#include <vector>
#include <string>
#include <fstream>

void ComputeSimilarityMatrix(
    std::vector<statistics::MixtureOfDiagonalGaussians> mog, 
    std::vector<std::vector<double> > &matrix)
{
  double max = 0;
  matrix.resize( mog.size() );
  for(unsigned int i = 0; i < mog.size(); ++i)
  {
    matrix[i].resize( mog.size() );
    for(unsigned int j = 0; j < mog.size(); ++j)
    {
      matrix[i][j] = mog[i].CSDivergence(mog[j]);
      if(matrix[i][j] > max)
        max = matrix[i][j];
    }
  }
  for(unsigned int i = 0; i < mog.size(); ++i)
    for(unsigned int j = 0; j < mog.size(); ++j)
      matrix[i][j] = 1 - (matrix[i][j] / max);
}

std::vector<int> ReadVector(std::string filename)
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
  return ret;
}

int main()
{
  statistics::DiagonalGaussian g;
  std::vector<double> mean, variance;

  variance.push_back(1); variance.push_back(1);
  mean.push_back(0); mean.push_back(0);
  for(int i = -20; i <= 20; i++)
  {
    statistics::MixtureOfDiagonalGaussians mog1, mog2;
    //set first mog
    mean[0] = 0; mean[1] = 0;
    g.Initialize(mean, variance); mog1.AddGaussian(g, 0.3);
    mean[0] = 3;
    g.Initialize(mean, variance); mog1.AddGaussian(g, 0.3);
    mean[0] = 8;
    g.Initialize(mean, variance); mog1.AddGaussian(g, 0.4);
    
    //set second mog
    mean[0] = 0; mean[1] = 0+i;
    g.Initialize(mean, variance); mog2.AddGaussian(g, 0.3);
    mean[0] = 3; mean[1] = 0+i;
    g.Initialize(mean, variance); mog2.AddGaussian(g, 0.3);
    mean[0] = 8; mean[1] = 0+i;
    g.Initialize(mean, variance); mog2.AddGaussian(g, 0.4);
    std::cout<<mog1.CSDivergence(mog2)<<std::endl;
    std::cout<<mog1.gaussian(0).Likelihood(mog1.gaussian(0).mean())<<std::endl;
  }

  //Testing Hmms from HTK
  statistics::HmmSet htk;
  //std::string filename("/people/hartmann/research/babel/FlatStartWords/hmm_training/cantonese_plp70/hmm9/hmmdefs");
  std::string filename("/people/hartmann/research/babel/FlatStartGrapheme/hmm_training/cantonese_plp_grapheme00/hmm39/hmmdefs");
  htk.LoadHtkHmmSet(filename);

  std::vector<statistics::MixtureOfDiagonalGaussians> mog = htk.states();
  std::vector<std::vector<std::string> > names = htk.mixture_names();
  std::vector<std::vector<double> > similarity;
  //ComputeSimilarityMatrix(mog, similarity);
  //utilities::SaveMatrixToTextFile(similarity, std::string("sim.txt"));
  //fileutilities::WriteBinaryPGM(similarity, std::string("simmx.pgm"));
  statistics::SpectralCluster sc;
  //sc.SetMatrix(similarity);
  //sc.ConvertToSimpleEigenvalue();
  //fileutilities::WritePlainPGM(sc.similarity_matrix(), 1024, std::string("normal.pgm"));

  unsigned int row, col;
  std::vector<std::vector<double> > mtx, U, S;
  std::vector<double> e;
  utilities::InitializeMatrix(mtx, 12.0, 6, 3);
  utilities::PrintMatrix(mtx, std::string("%f"));
  utilities::CreateIdentityMatrix(mtx, 12, 12);
  utilities::PrintMatrix(mtx, std::string("%f"));
  mtx.resize(4);
  mtx[0].resize(4); mtx[0][0] = 4; mtx[0][1] = -30; mtx[0][2] = 60; mtx[0][3] = -35;
  mtx[1].resize(4); mtx[1][0] = -30; mtx[1][1] = 300; mtx[1][2] = -675; mtx[1][3] = 420;
  mtx[2].resize(4); mtx[2][0] = 60; mtx[2][1] = -675; mtx[2][2] = 1620; mtx[2][3] = -1050;
  mtx[3].resize(4); mtx[3][0] = -35; mtx[3][1] = 420; mtx[3][2] = -1050; mtx[3][3] = 700;
  //mtx.resize(3);
  //mtx[0].resize(3); mtx[0][0] = 1; mtx[0][1] = 0.05; mtx[0][2] = 0.5;
  //mtx[1].resize(3); mtx[1][0] = .05; mtx[1][1] = 1; mtx[1][2] = 0.83;
  //mtx[2].resize(3); mtx[2][0] = 0.5; mtx[2][1] = 0.83; mtx[2][2] = 1;
  std::cout<<utilities::FrobeniusNorm(mtx)<<std::endl;
  std::cout<<utilities::MaxElementInMatrix(mtx, row, col)<<std::endl;
  utilities::PrintMatrix(mtx, std::string("%f"));

 /* 
  sc.SetMatrix(mtx);
  sc.PerformEigenAnalysis(U, e);
  //U = utilities::TransposeMatrix(U);
  //U = utilities::MatrixProduct(mtx, mtx);
  utilities::PrintMatrix(U, std::string("%f"));
  for(unsigned int i = 0; i < e.size(); ++i)
    std::cout<<i<<" "<<e[i]<<std::endl;

  int K = 100;
  S = U;
  for(unsigned int i = 0; i < U.size(); ++i)
  {
    S[i].resize(K);
    for(unsigned int k = 0; k < K; ++k)
      S[i][k] = U[k][i];
  }

  statistics::KMeans km;
  km.Initialize(U, K);
  km.PerformKMeans();*/
  std::vector<int> cluster;// = km.clusters();
  int K = 100;
  cluster = ReadVector(std::string("idx.txt"));
  for(unsigned int i = 0; i < cluster.size(); ++i)
    std::cout<<cluster[i]<<std::endl;
  std::cout<<cluster.size()<<std::endl;
  // We have a problem with the ending, should be fixed later
  cluster.resize(1521);

  for( int k = 0; k < K; ++k)
  {
    std::cout<<k<<":";
    for(unsigned int c = 0; c < cluster.size(); ++c)
    {
      if(cluster[c] == k)
      {
        for(unsigned int i = 0; i < names[c].size(); ++i)
        {
          //std::cout<<" "<<c<<" "<<i<<" "<<names.size()<<std::endl;
          std::cout<<" "<<names[c][i];
        }
      }
    }
    std::cout<<std::endl;
  }

  return 0;
}
