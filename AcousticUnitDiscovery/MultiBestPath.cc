// William Hartmann (hartmannw@gmail.com)
// This is free and unencumbered software released into the public domain.
// See the UNLICENSE file for more information.

#include "MultiBestPath.h"

namespace acousticunitdiscovery
{

std::vector<int> FindBestPath(const utilities::Matrix<double> &pgram,
    const utilities::Matrix<double> &transition, int min_frames)
{
  std::vector<int> path; // Final path that we return.
  utilities::Matrix<ViterbiInfo> dp_matrix; // Holds the memoization
                                                    // data.
  unsigned int states = pgram.NumRows();
  unsigned int frames = pgram.NumCols();
  double minimum_log = -1000000; // Essentially represents log(0). Used for
                                 // states that should be unreachable.

  dp_matrix.Initialize(states * min_frames, frames );

  // Set the values for the first frame in the dp_matrix. We are expanding the
  // initial number of states by the minimum number of frames a state must span.
  // At the first frame only the first state of each original state can be
  // reached.
  for(unsigned int i = 0; i < states; ++i)
  {
    unsigned int index = i * min_frames;
    dp_matrix(index,0).score = pgram(i,0);
    dp_matrix(index,0).parent = -1;
    for(unsigned int j = 1; j < static_cast<unsigned int>(min_frames - 1); ++j)
    {
      dp_matrix(index+j,0).score = minimum_log;
      dp_matrix(index+j,0).parent = -1;
    }
  }
  
  // Fill in the dp_matrix for each frame.
  for(unsigned int f = 1; f < frames; ++f)
  {
    for(unsigned int s = 0; s < states; ++s)
    {
      unsigned int index = s * min_frames;
      // For the first sub-state it can be reached only from the final sub-state
      // of any possible state.
      double best_score = minimum_log;
      int best_parent = -1;
      for(unsigned int p = 0; p < states; ++p)
      {
        unsigned int parent_index = (p * min_frames) + (min_frames - 1);
        double score = dp_matrix(p,f-1).score + transition(s,p) + pgram(s,f);
        if(score > best_score)
        {
          best_score = score;
          best_parent = parent_index;
        }
      }
      dp_matrix(index,f).score = best_score;
      dp_matrix(index,f).parent = best_parent;

      // All intermediary sub-states can only come from the previous sub-state.
      for(unsigned int i = 1; i < static_cast<unsigned int>(min_frames); ++i)
      {
        dp_matrix(index+i,f).parent = index+i-1;
        dp_matrix(index+i,f).score = dp_matrix(index+i-1,f-1).score + 
            pgram(s,f) + transition(s,s);
      }

      // Final sub-state is actually allowed to do a self loop. Check if self
      // loop provides a better score;
      double self_loop = dp_matrix(index+min_frames-1,f-1).score + pgram(s,f)
          + transition(s,s);
      if(self_loop > dp_matrix(index+min_frames-1,f).score)
      {
        dp_matrix(index+min_frames-1,f).score = self_loop;
        dp_matrix(index+min_frames-1,f).parent = index+min_frames-1;
      }

    }
  }
  return BestPathInDpMatrix(dp_matrix, min_frames);
}

std::vector<int> BestPathInDpMatrix(
    const utilities::Matrix<ViterbiInfo> &dp_matrix, unsigned int min_frames)
{
  std::vector<int> ret;
  unsigned int states = dp_matrix.NumRows();
  unsigned int frames = dp_matrix.NumCols();

  // Find the ending point
  ViterbiInfo endpoint = dp_matrix(min_frames-1, frames-1);
  endpoint.parent = min_frames-1;
  for(unsigned int s = endpoint.parent; s < states; s+=min_frames)
    if( dp_matrix(s, frames-1).score > endpoint.score)
    {
      endpoint.parent = s;
      endpoint.score = dp_matrix(s, frames-1).score;
    }
  //ret.push_back(endpoint.parent);
  //endpoint = dp_matrix(endpoint.parent, frames -1);
  for(unsigned int f = frames-2; f >= 0; --f)
  {
    ret.push_back(endpoint.parent);
    endpoint = dp_matrix(endpoint.parent, f);
  }
  return ret;
}

utilities::Matrix<double> GenerateTransitionMatrix(
    unsigned int states, double self_loop_prob)
{
  utilities::Matrix<double> ret;
  
  double log_self_loop_prob = std::log(self_loop_prob);
  double log_trans_prob = std::log(1 - self_loop_prob);

  ret.Initialize(states, states, log_trans_prob);
  ret.SetDiagonal(log_self_loop_prob);

  return ret;
}

}
