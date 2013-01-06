// William Hartmann (hartmannw@gmail.com)
// This is free and unencumbered software released into the public domain.
// See the UNLICENSE file for more information.
//
// Definition for the DynamicTimeWarp class.  The general usage for the class 
// is to read in two feature vectors, compute a similarity matrix, and finally
// compute either the DTW path or the set of Segmental DTW paths.  The
// segmental DTW algorithm performs DTW multiple times with different start 
// locations.  Each path is also constrained to a specific window size.  The
// paths can be further refined by only keeping the best subsequence of at 
// least a certain length L within the original path.  For more information
// about the segmental DTW algorithm see: "Unsupervised Pattern Discovery in 
// Speech" by Alex S. Park and James R. Glass, 2008.

#ifndef ACOUSTICUNITDISCOVERY_DYNAMICTIMEWARP_H_
#define ACOUSTICUNITDISCOVERY_DYNAMICTIMEWARP_H_

#include<vector>
#include<cmath>
#include<limits>
#include<algorithm>

#include "Matrix.h"
#include "MatrixFunctions.h"
#include "ImageIO.h" // Functions for writing the similarity matrix and paths
                     // as an image.

namespace acousticunitdiscovery
{

// When computing the best path through the similarity matrix, we need to
// maintain a way to move back along the path.  TrackBackDirection describes 
// whether to move along the first dimension, second dimension, or both.  It 
// also has a value for when we reach the origin (starting location) or an 
// invalid location.
enum TrackBackDirection {INVALID, FIRST_DIMENSION, SECOND_DIMENSION, DIAGONAL, 
    ORIGIN};

// PathPoint stores a particular point in a path.  The score should be
// identical to the value stored in the similarity matrix at index
// [first][second].
typedef struct
{
  unsigned int first; //index into first feature vector
  unsigned int second;//index into second feature vector
  double score;       //distance between the two points based on the similarity
                      //matrix and some distance metric
} PathPoint;

//Represents an entire path through the similarity matrix.
typedef struct
{
  std::vector< PathPoint > path;
  double total_score;
} DtwPath;

// Stores the data and functions required for computing either the standard DTW
// or the segmental DTW.
class DynamicTimeWarp
{
 public:
  DynamicTimeWarp(){}
  ~DynamicTimeWarp(){}

  // Stores the utterances
  void set_utterance_one( const utilities::Matrix<double> &utterance){ 
      utterance_one_ = utterance;}
  void set_utterance_two( const utilities::Matrix<double> &utterance){
      utterance_two_ = utterance;}

  // Access functions
  std::vector<DtwPath> paths(){ return paths_;}

  // Creates a matrix of size length(utterance_one) x length(utterance_two).
  // Each point [i][j] stores the distance between the feature vector at frame
  // i of utterance_one and frame j of utterance_two.  Must be called before 
  // making calls to any of the path finding functions.
  bool ComputeSimilarityMatrix();

  // Computes the best path through the similarity matrix starting at point
  // [0][0] and ending at [length(utterance_one)][length(utterance_two)].  The 
  // computed path is added to the paths_ variable.
  bool ComputeStandardDTW();

  // Computes multiple paths through the similarity matrix.  Each start point 
  // is separated by an interval of (2*constraint)+1 along the boundary of the
  // similarity matrix.  The end point is the point along the diagonal on the 
  // opposite end of the similarity matrix.  All paths must stay within 
  // constraint units of the diagonal.  No two paths can ever overlap.
  bool ComputeSegmentalDTW(unsigned int constraint);

  // Calls the relevant function in ImageIO.h to convert the similarity matrix
  // to a PGM image.  Any paths along the similarity matrix are also shown as 
  // white, the maximum value in the image.
  bool SaveResultAsPGM( std::string filename );

  // Takes any paths in the variable path_ and finds the best subsequence of 
  // length at least minlength.  Park and Glass, 2008, found this to be too
  // conservative and also expanded the edges of the path such that they are 
  // still within a factor of the original average score.  This factor is set 
  // by expansion_factor.  If you do not want to expand the paths, set 
  // expansion_factor less than 0.  Due to the implementation a value of 0 will
  // always expand the path by at least one point.  Also note that if any paths
  // in the original set are less than minlength in size, they will be pruned
  // away completely.
  bool PrunePathsByLCMA(const unsigned int minlength, double expansion_factor);

  // If silence exists in the two utterances, they will recieve a very good 
  // score (low distance).  This behavior is likely unwanted.  Given a boolean
  // vector silence, the length of the first utterance, any path point which 
  // correponds to silence will be given a maximal score.  This maximal score 
  // is determined by looking through the first path and selecting the highest 
  // score.
  bool IncreaseSilenceCost(const std::vector<bool> silence);

  // Returns a vector where each entry corresponds to a frame of utterance_one.
  // The values of the entries are equal to the cost for the best pass that 
  // goes through that frame.  Any frames that do not have a path are given the
  // score for the highest cost path in the set.
  std::vector<double> BestScorePerFrame();
 
 private:
  
  utilities::Matrix<double> utterance_one_; // Utterances are assumed to be
  utilities::Matrix<double> utterance_two_; // organized as frames x features,
                                            // so the first index is into the 
                                            // frame.

  // Stores the distances between the feature vectors for every pair of frames 
  // in utterance_one and utterance_two.
  utilities::Matrix<double> similarity_matrix_;
  std::vector< DtwPath > paths_;  // All computed paths are stored here.

  // Computes the distance between the two given feature vectors.  Currently 
  // the only supported distance metric is Euclidean distance.  Function is 
  // only used to compute the similarity matrix.
  double GetFeatureDistance(const std::vector<double> &one, 
      const std::vector<double> &two);

  // Computes a single DTW path based from startpoint to endpoint.  All points
  // in the path must be within constraint points of the diagonal.  It is 
  // possible to set the endpoint outside of the area covered by the constraint
  // and will result in no path being found.
  bool DTW(const PathPoint &startpoint, 
    const PathPoint &endpoint, const unsigned int &constraint);

  // During the computation of the best path, determines where the best origin
  // for the current_point is.  It considers both the start_point and 
  // constraint in its selection.  direction stores the result.  dp_matrix is 
  // the dynamic programming matrix and stores the cost of the best path to any
  // given point computed so far.
  bool SetBestOrigin(PathPoint start_point, PathPoint current_point, 
      unsigned int constraint, utilities::Matrix<double> &dp_matrix,
      TrackBackDirection &direction,
      const utilities::Matrix<TrackBackDirection> &trackback);

  // Identifies whether a point can be reached given the starting location and
  // constraint.
  bool PointWithinConstraint(const PathPoint &start_point, unsigned int first, 
      unsigned int second, unsigned int constraint);

  // From the given dp_matrix, storing the best path to any given point from
  // the starting point, the best path to the endpoint is found.  The path is 
  // automatically added to paths_.
  bool AddBestPath( const utilities::Matrix<double> &dp_matrix, 
      const utilities::Matrix<TrackBackDirection> &backtrack_matrix,
      const PathPoint &endpoint);

  // Length Constrained Minimum Average (LMCA) susbsequence finds the best 
  // sub-path within a given path.  Information about the path is stored in 
  // section where section.first is the start, section.second is the end, and 
  // section.score stores the average cost of the path.
  bool LCMA(const DtwPath &path, unsigned int minlength, PathPoint &section);

  // Assuming a path has already been pruned, the path is expanded such that 
  // the addition of the points does not increase the average cost by more than
  // expansion_factor * original_average_cost.
  bool ExtendPath(const DtwPath &path, double expansion_factor,
      PathPoint &section);

};

}// end namespace acousticunitdiscovery

#endif
