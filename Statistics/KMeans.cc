#include "KMeans.h"

namespace statistics
{

void KMeans::Initialize(std::vector<std::vector<double> > data, unsigned int K)
{
  data_ = data;
  mean_.resize(K);
  cluster_.resize(data.size());

  InitializeMeans();
}

void KMeans::InitializeMeans()
{
  for(unsigned int i = 0; i < mean_.size(); ++i)
  {
    mean_[i].resize(data_[0].size());
    for(unsigned int j = 0; j < mean_[i].size(); ++j)
    {
      unsigned int indexi = std::rand() % data_.size();
      unsigned int indexj = std::rand() % data_[indexi].size();
      mean_[i][j] = data_[indexi][indexj];
    }
  }
}

void KMeans::UpdateMeans()
{
  std::vector<std::vector<double> > newmean = mean_;
  std::vector<int> count;
  count.resize(mean_.size(), 0);
  for(unsigned int i = 0; i < newmean.size(); ++i)
    for(unsigned int j = 0; j < newmean[i].size(); ++j)
      newmean[i][j] = 0;

  for(unsigned int i = 0; i < cluster_.size(); ++i)
  {
    count[cluster_[i]]++;
    for(unsigned int j = 0; j < mean_[0].size(); ++j)
    {
      newmean[cluster_[i]][j] += data_[i][j];
    }
  }

  for(unsigned int i = 0 ; i < mean_.size(); ++i)
  {
    if( count[i] == 0)
    {
      mean_[i] = mean_[0];
    }
    else
    {
      for(unsigned int j = 0; j < mean_[0].size(); ++j)
        mean_[i][j] = newmean[i][j] / count[i];
    }
  }
}

bool KMeans::AssignToClusters()
{
  std::vector<int> newcluster = cluster_;
  std::vector<double> distance;
  distance.resize(cluster_.size());
  bool changed = false;

  for(unsigned int i = 0; i < cluster_.size(); ++i)
  {
    newcluster[i] = 0;
    distance[i] = PointDistance(mean_[0], data_[i]);
    for(unsigned int j = 1; j < mean_.size(); ++j)
    {
      double current_distance = PointDistance(mean_[j], data_[i]);
      if( current_distance < distance[i] )
      {
        distance[i] = current_distance;
        newcluster[i] = j;
      }
    }
    if(cluster_[i] != newcluster[i])
      changed = true;
  }
  cluster_ = newcluster;
  return changed;

}

void KMeans::PerformKMeans()
{
  int iter = 0;
  int maxiter = 1000;

  while(iter < maxiter && AssignToClusters())
  {
    UpdateMeans();
    iter++;
  }
}

std::vector<int> KMeans::clusters()
{
  return cluster_;
}

double KMeans::PointDistance(std::vector<double> &a, std::vector<double> &b)
{
  double ret = 0;
  for(unsigned int i = 0; i < a.size(); ++i)
    ret += pow(a[i] - b[i], 2);
  return std::sqrt(ret);
}

}
