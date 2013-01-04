#ifndef STATISTICS_KMEANS_
#define STATISTICS_KMEANS_

#include<vector>
#include<cstdlib>
#include<cmath>

namespace statistics
{

class KMeans
{
 private:
  std::vector<std::vector<double> > data_;
  std::vector<std::vector<double> > mean_;
  std::vector<int> cluster_;

  double PointDistance(std::vector<double> &a, std::vector<double> &b);

 public:
  KMeans(){}
  ~KMeans(){}
  void Initialize(std::vector<std::vector<double> > data, unsigned int K);
  void InitializeMeans();
  void UpdateMeans();
  bool AssignToClusters();
  void PerformKMeans();

  std::vector<int> clusters();
};

}

#endif
