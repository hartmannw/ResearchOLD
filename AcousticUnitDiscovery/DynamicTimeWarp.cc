// William Hartmann (hartmannw@gmail.com)
//
// Implementation of the DynamicTimeWarp class. For a detailed description of 
// the class, see the corresponding .h file.

#include "DynamicTimeWarp.h"

namespace acousticunitdiscovery
{

bool DynamicTimeWarp::ComputeSimilarityMatrix()
{
  if(utterance_one_.NumRows() < 1 || utterance_two_.NumRows() < 1)
    return false; // We must have two utterances.
  
  similarity_matrix_.Initialize(utterance_one_.NumRows(), 
      utterance_two_.NumRows(), 0);
  for(unsigned int first = 0; first < utterance_one_.NumRows(); ++first)
  {
    for(unsigned int second = 0; second < utterance_two_.NumRows(); ++second)
      similarity_matrix_(first, second) = GetFeatureDistance(
          utterance_one_.GetRow(first), utterance_two_.GetRow(second) );
  }
  return true;
}

bool DynamicTimeWarp::ComputeStandardDTW()
{
  PathPoint start_point, end_point;
  unsigned int constraint;

  if(similarity_matrix_.NumRows() < 1 || similarity_matrix_.NumCols() < 1)
    return false; //similarity_matrix has not been computed.
    
  // Set constraint such that all paths through the similarity matrix are 
  // possible.
  constraint = similarity_matrix_.NumRows() + similarity_matrix_.NumCols();
  start_point.first = 0;  // Assumes we want a path from corner to corner in 
  start_point.second = 0; // the similarity matrix.
  end_point.first = similarity_matrix_.NumRows() - 1;
  end_point.second = similarity_matrix_.NumCols() - 1;

  // Actual logic for computing the best path is in DTW.
  return DTW(start_point, end_point, constraint);
}

bool DynamicTimeWarp::ComputeSegmentalDTW(unsigned int constraint)
{
  unsigned int diagonal;
  // Move the starting point along the first column of the similarity matrix.
  for(unsigned int r = 0; r < similarity_matrix_.NumRows(); 
      r+= (constraint*2)+1)
  {
    PathPoint start_point, end_point;
    start_point.first = r; start_point.second = 0;
    // Note that if the similarity matrix is not square, the first path
    // computed may not be the same as the standard DTW.
    diagonal = similarity_matrix_.NumRows() - r;
    if(similarity_matrix_.NumCols() < diagonal) // Assumes we want the end
      diagonal = similarity_matrix_.NumCols();  // point to be along the
    end_point.first = (r + diagonal - 1);       // from the start point.
    end_point.second = (diagonal - 1);
    DTW(start_point, end_point, constraint);
  }
  // Now move the starting point along the first row of the similarity matrix.
  // Avoids repeating the path with start point [0][0].
  for(unsigned int c = (constraint*2)+1; c < similarity_matrix_.NumCols(); 
      c+= (constraint*2)+1)
  {
    PathPoint start_point, end_point;
    start_point.first = 0; start_point.second = c;
    diagonal = similarity_matrix_.NumCols() - c;
    if(similarity_matrix_.NumRows() < diagonal)
      diagonal = similarity_matrix_.NumRows();
    end_point.first = diagonal -1;
    end_point.second = c + diagonal - 1;
    DTW(start_point, end_point, constraint);
  }
  return true;
}

bool DynamicTimeWarp::SaveResultAsPGM(std::string filename)
{
  double maxvalue = utilities::MaxElementInMatrix(similarity_matrix_);

  // Sets the value for any point in the similarity matrix that corresponds to
  // a path to the maximum value.  This achieves the effect of making the paths
  // white in the resulting image.
  std::vector< std::vector<double> > simmx = 
      similarity_matrix_.GetVectorOfVectors();
  for(unsigned int i = 0; i < paths_.size(); i++)
    for(unsigned int j = 0; j < paths_[i].path.size(); j++)
      simmx[ paths_[i].path[j].first ][ paths_[i].path[j].second] = maxvalue;

  // Writes simmx as an image to filename.  Values are linearly scaled to be 
  // between 0 and 255.
  return fileutilities::WriteBinaryPGM(simmx, filename);
}

bool DynamicTimeWarp::PrunePathsByLCMA(const unsigned int minlength, 
    double expansion_factor)
{
  std::vector<int> to_remove;
  for(unsigned int i = 0; i < paths_.size(); i++)
  {
    PathPoint section;
    if(LCMA(paths_[i], minlength, section))
    {
      ExtendPath(paths_[i], expansion_factor, section);
      // If we have removed points from the end, physically remove them from
      // the path.
      if(section.second+1 < paths_[i].path.size())
        paths_[i].path.erase(paths_[i].path.begin()+section.second+1, 
           paths_[i].path.end());
      // If we have removed points from the start, physically remove them from 
      // the path.
      if(section.first > 0)
        paths_[i].path.erase(paths_[i].path.begin(),
           paths_[i].path.begin()+section.first-1);
      paths_[i].total_score = section.score;
    }
    else // Path was too short; mark it for removal.
    {
      to_remove.push_back(i);
    }
  }
  // Remove paths from back to front as the opposite direction will affect 
  // indexing.
  for(int i = to_remove.size() - 1; i >= 0; --i)
    paths_.erase(paths_.begin() + to_remove[i]);

  return true;
}

bool DynamicTimeWarp::IncreaseSilenceCost(std::vector<bool> silence)
{
  if( paths_.size() < 1 )
    return false;
  // Find the maximal score in the first path.
  double max_score = 0;
  for(unsigned int i=0; i < paths_[0].path.size(); i++)
    if( paths_[0].path[i].score > max_score)
      max_score = paths_[0].path[i].score;

  for(unsigned int i=0; i < paths_.size(); i++)
    for(unsigned int j=0; j < paths_[i].path.size(); j++)
      if( silence[ paths_[i].path[j].first ] )
      {
        paths_[i].path[j].score = max_score;
      }

  return true;
}

std::vector<double> DynamicTimeWarp::BestScorePerFrame()
{
  std::vector<double> result;
  double max_value = 0;
  result.resize(utterance_one_.NumRows(), 0);
  for(unsigned int i=0; i < paths_.size(); i++)
  {
    double total_score = paths_[i].total_score;
    if(total_score > max_value)
      max_value = total_score;
    for(unsigned int j = 0; j < paths_[i].path.size(); j++)
    {
      int frame = paths_[i].path[j].first;
      if( result[frame] == 0 || result[frame] > total_score)
        result[frame] = total_score;
    }
  }
  for(unsigned int i=0; i < result.size(); i++)
    if(result[i] == 0)
      result[i] = max_value;
  return result;
}

/*
double DynamicTimeWarp::GetFeatureDistance(const std::vector<double> &one,
    const std::vector<double> &two)
{
  double distance = 0;
  //for now we will use Euclidean distance
  for(unsigned int i =0; i < one.size(); ++i)
    distance += pow(one[i] - two[i], 2);
  distance = sqrt(distance);
  return distance;
}
*/

double DynamicTimeWarp::GetFeatureDistance(const std::vector<double> &one,
    const std::vector<double> &two)
{
  double distance = 0;
  double one_magnitude = 0;
  double two_magnitude = 0;
  //for now we will use Cosine similarity
  for(unsigned int i =0; i < one.size(); ++i)
  {
    distance += one[i] * two[i];
    one_magnitude += one[i] * one[i];
    two_magnitude += two[i] * two[i];
  }
  distance = 1 - (((distance / ( sqrt(one_magnitude) * 
            sqrt(two_magnitude) ))+1)/2);
  return distance;
}

bool DynamicTimeWarp::DTW(const PathPoint &startpoint,
    const PathPoint &endpoint, const unsigned int &constraint)
{
  utilities::Matrix<double> dp_matrix; // dp as in dynamic programming
  utilities::Matrix<TrackBackDirection> backtrack_matrix;
  PathPoint current_point;

  dp_matrix = similarity_matrix_;

  backtrack_matrix.Initialize(dp_matrix.NumRows(), dp_matrix.NumCols());
  backtrack_matrix.SetCol(0, INVALID);
  
  for(unsigned int r = startpoint.first; r < similarity_matrix_.NumRows(); ++r)
  {
    for(unsigned int c = startpoint.second; c < similarity_matrix_.NumCols();
        ++c) 
    {
      if(PointWithinConstraint(startpoint,r,c,constraint))
      {
        current_point.first = r;
        current_point.second = c;
        SetBestOrigin(startpoint, current_point, constraint, dp_matrix, 
            backtrack_matrix(r,c), backtrack_matrix);
      }
    }
  }
  AddBestPath(dp_matrix, backtrack_matrix, endpoint);
  return true;
}

bool DynamicTimeWarp::AddBestPath(
    const utilities::Matrix<double> &dp_matrix,
    const utilities::Matrix<TrackBackDirection> &backtrack_matrix,
    const PathPoint &endpoint)
{
  unsigned int r = endpoint.first;
  unsigned int c = endpoint.second;
  DtwPath path;
  PathPoint point;

  point.first = r;
  point.second = c;
  point.score = similarity_matrix_(r,c);
  path.path.push_back(point);
  path.total_score = dp_matrix(r,c);

  while(backtrack_matrix(r,c) != ORIGIN)
  {
    if(r < 0 || c < 0) // Best path leads off the similarity matrix.
      return false;

    if(backtrack_matrix(r,c) == FIRST_DIMENSION)
    {
      r = r - 1;
    }
    else if(backtrack_matrix(r,c) == SECOND_DIMENSION)
    {
      c = c -1;
    }
    else if(backtrack_matrix(r,c) == DIAGONAL)
    {
      r = r - 1;
      c = c - 1;
    }
    else // Path leads to an invalid point.  No valid path exists.
    {
      return false;
    }
    PathPoint point;
    point.first = r;
    point.second = c;
    point.score = similarity_matrix_(r,c);
    path.path.push_back(point);
  }
  // Points were added to vector in reverse order.
  reverse(path.path.begin(), path.path.end());
  paths_.push_back(path);
  
  return true;
}

bool DynamicTimeWarp::SetBestOrigin(PathPoint start_point,
    PathPoint current_point, 
    unsigned int constraint, utilities::Matrix<double> &dp_matrix,
    TrackBackDirection &direction, 
    const utilities::Matrix<TrackBackDirection> &trackback)
{
  double best_score = std::numeric_limits<double>::max();
  direction = INVALID;

  // Check if we are at the origin.
  if( current_point.first == start_point.first && 
      current_point.second == start_point.second)
  {
    direction = ORIGIN;
    return true;
  }
  // Check the first dimension.
  if( static_cast<int>(current_point.first) > 0 &&
      trackback(current_point.first-1,current_point.second) != INVALID &&
      dp_matrix(current_point.first-1,current_point.second) < best_score)
  {
    direction = FIRST_DIMENSION;
    best_score = dp_matrix(current_point.first-1,current_point.second); 
  }
  // Check the second dimension.
  if( static_cast<int>(current_point.second) > 0 && 
      trackback(current_point.first,current_point.second-1) != INVALID &&
      dp_matrix(current_point.first,current_point.second-1) < best_score)
  {
    direction = SECOND_DIMENSION;
    best_score = dp_matrix(current_point.first,current_point.second-1);
  }
  // Check along the diagonal.
  if( static_cast<int>(current_point.first) > 0 && 
      static_cast<int>(current_point.second) > 0 && 
      trackback(current_point.first-1,current_point.second-1) != INVALID &&
      dp_matrix(current_point.first-1,current_point.second-1) < best_score)
  {
    direction = DIAGONAL;
    best_score = dp_matrix(current_point.first-1,current_point.second-1);
  }
  if(direction == INVALID) // We could not find a valid point which leads to 
  {                        // the current point, so the current point must also
    return false;          // be invalid.
  }
  dp_matrix(current_point.first,current_point.second)+= best_score;
  return true;  
}

bool DynamicTimeWarp::PointWithinConstraint(const PathPoint &start_point, 
    unsigned int first, unsigned int second, unsigned int constraint)
{
  int distance;
  // Static_cast avoids the undefined behavior of a negative unsigned int.
  distance = static_cast<int>( first - start_point.first );
  distance -= static_cast<int>(second - start_point.second);
  distance = std::abs(distance);
  return (distance <= static_cast<int>(constraint));
}

// Note that asymptotically faster implementations of this algorithm exist.
// However, the potential speed improvements did not seem worth it at this
// time, especially after the performance improvements from constraining the 
// maximum length and keeping a running mean.  For more information on this 
// algorithm see "Efficient algorithms for locating the length-constrained 
// heaviest segments with applications to biomolecular sequence analysis" by 
// Yaw-Ling Lin, Tao Jiang, and Kun-Mao Chao, 2002.
bool DynamicTimeWarp::LCMA(const DtwPath &path, unsigned int minlength,
    PathPoint &section)
{
  section.first=0; section.second=0; 
  section.score=std::numeric_limits<double>::max();
  for(unsigned int s = 0; s < path.path.size(); s++)
  {
    double score = 0;
    // We do not allow any paths longer than 3 x the minimum length.  This 
    // significantly speeds up performance.  In my experience it is unlikely 
    // for these types of paths to exist anyway.  This constraint could be set 
    // lower, such as 2 x minimum length, but that only produces a very small
    // performance improvement.  If the paths are later extended, then the 
    // longer path will still be found.
    for(unsigned int e = s+minlength-1; 
        e < std::min(static_cast<unsigned int>(path.path.size())
        ,s+(3*minlength)); e++)
    {
      // Only calculate the mean for the first iteration of this loop.  Keeping
      // a running calculation of the mean provides a big speed improvement.
      if(score == 0)
      {
        for(unsigned int i = s; i <= e; i++)
        {
          score += path.path[i].score;
        }
        score = score / (e-s+1);
      }
      else
      {
        double length = (e-s);
        score = (length/(length+1)) * score;
        score += (1 / (length+1)) * path.path[e].score;
      }
      if(score < section.score)
      {
        section.first = s; section.second = e;
        section.score = score;
      }
    }
  }
  if(section.first == 0 && section.second == 0) // Only happens if path is 
    return false;                               // already too short.
  return true;
}

bool DynamicTimeWarp::ExtendPath(const DtwPath &path, double expansion_factor,
    PathPoint &section)
{
  int length = section.second - section.first + 1;
  double max_score = (1 + expansion_factor) * section.score;
  bool extend_start = false;  // Determines whether to take a step in the start
                              // or end direction.
  while(section.score < max_score)
  {
    if( section.first == 0 && section.second == (path.path.size() -1) )
    { // The expanded path has reached the original path.
      return true;
    }
    else if( section.first == 0)
    { // No more room to extend at the front.
      extend_start = false;
    }
    else if( section.second == (path.path.size() - 1))
    { // No more room to extend at the end.
      extend_start = true;
    }
    // Make a greedy selection between extending at the front or the back.
    // Implementation favors the point with the lowest cost.
    else if( path.path[section.first-1].score < 
        path.path[section.second+1].score)
    {
      extend_start = true;
    }
    else
    {
      extend_start = false;
    }

    if(extend_start)
    {
      section.first -= 1;
      // Maintains running score instead of recomputing from scratch.
      section.score = ( (length/static_cast<double>(length+1)) * section.score)
          + ( (1/static_cast<double>(length+1)) * 
          path.path[section.first].score);
    }
    else
    {
      section.second += 1;
      section.score = ( (length/static_cast<double>(length+1)) * section.score)
          + ( (1/static_cast<double>(length+1)) * 
          path.path[section.second].score);
    }
    length = section.second - section.first + 1;
  }
  return true;  
}

} //end namespace acousticunitdiscovery
