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
  utilities::Matrix<ViterbiInfo> dp_matrix; // Holds the memoization data.
  unsigned int states = pgram.NumRows();
  unsigned int frames = pgram.NumCols();
  double minimum_log = -1000000; // Essentially represents log(0). Used for
                                 // states that should be unreachable.

  ViterbiInfo default_value;
  default_value.parent = -1;
  default_value.score = minimum_log;
  dp_matrix.Initialize(states * min_frames, frames, default_value );

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
        double score = dp_matrix(parent_index,f-1).score + transition(s,p) + 
            pgram(s,f);
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
  utilities::Matrix<double> score;
  score.Initialize(dp_matrix.NumRows(), dp_matrix.NumCols());
  for(unsigned int r = 0; r < score.NumRows(); ++r)
    for(unsigned int c = 0; c < score.NumCols(); ++c)
      score(r,c) = 1 * std::max(dp_matrix(r,c).parent, -1000);
  fileutilities::WriteBinaryPGM(score.GetVectorOfVectors(), 
      std::string("dp1.pgm"));

  std::vector<int> initial_path; 
  return BestPathInDpMatrix(dp_matrix, min_frames, initial_path);
}

std::vector<int> FindViterbiPath(const utilities::Matrix<double> &pgram,
    const utilities::Matrix<double> &transition, int min_frames)
{
  std::vector<int> initial_path;
  initial_path.resize(0); // Make sure the path is of size 0.
  return FindRestrictedViterbiPath(pgram , transition, min_frames, 
      initial_path);
}

// The number of states in the dynamic programming matrix is expanded by the 
// minimum number of states required for each state. The logic to handle keeping
// track of the indices has been handed off to the GetStateScore and 
// GetTransitionScore functions. 
std::vector<int> FindRestrictedViterbiPath(
    const utilities::Matrix<double> &pgram,
    const utilities::Matrix<double> &transition, int min_frames, 
    std::vector<int> initial_path)
{
  std::vector<int> path; // Final path that we return.
  utilities::Matrix<ViterbiInfo> dp_matrix; // Holds the memoization data.
  unsigned int states = (pgram.NumRows() + initial_path.size()) * min_frames;
  unsigned int frames = pgram.NumCols();
  double minimum_log = -1000000; // Essentially represents log(0). Used for
                                 // states that should be unreachable.

  // We initialize the the dp_matrix. Initially, every state is set with the
  // minimum score and considered inaccessible. If there is an initial path 
  // restriction, then only the first overall state is a valid state state.
  // Otherwise, the initial substate for every state is valid.
  ViterbiInfo default_value;
  default_value.parent = -1;
  default_value.score = minimum_log;
  dp_matrix.Initialize(states, frames, default_value);

  if(initial_path.size() > 0) // Only first overall state is a valid start
  {                           // state.
    dp_matrix(0,0).score = std::max(
        GetStateScore(pgram, initial_path, min_frames, 0, 0), minimum_log);
  }
  else // The first substate of every original state is a valid start state.
  {
    for(unsigned int i = 0; i < states; i+= min_frames)
      dp_matrix(i, 0).score = std::max(
          GetStateScore(pgram, initial_path, min_frames, i, 0), minimum_log);
  }
  
  // Fill in the remainder of the dp_matrix
  for(unsigned int f = 1; f < frames; ++f)
  {
    for(unsigned int s = 0; s < states; ++s)
    {
      ViterbiInfo best_point;
      best_point.parent = 0;
      best_point.score = dp_matrix(0, f-1).score + std::max(
          GetTransitionScore(transition, initial_path, min_frames, 0, s),
          minimum_log);
      for(unsigned int p = 1; p < states; ++p)
      {
        double score = dp_matrix(p, f-1).score + std::max(
            GetTransitionScore(transition, initial_path, min_frames, p, s),
            minimum_log);
        if(score > best_point.score)
        {
          best_point.score = score;
          best_point.parent = p;
        }
      } // end for p
      best_point.score += std::max(
          GetStateScore(pgram, initial_path, min_frames, s, f), minimum_log);
      dp_matrix(s, f) = best_point;
    } // end for s
  } // end for f

  utilities::Matrix<double> score;
  score.Initialize(dp_matrix.NumRows(), dp_matrix.NumCols());
  for(unsigned int r = 0; r < score.NumRows(); ++r)
    for(unsigned int c = 0; c < score.NumCols(); ++c)
      score(r,c) = 1 + std::max(dp_matrix(r,c).parent, -1000);
  fileutilities::WriteBinaryPGM(score.GetVectorOfVectors(), 
      std::string("dp2.pgm"));
  
  return BestPathInDpMatrix(dp_matrix, min_frames, initial_path);
}


std::vector<int> BestPathInDpMatrix(
    const utilities::Matrix<ViterbiInfo> &dp_matrix, unsigned int min_frames, 
    const std::vector<int> &initial_path)
{
  std::vector<int> ret;
  unsigned int states = dp_matrix.NumRows();
  unsigned int frames = dp_matrix.NumCols();

  // Find the ending point
  ViterbiInfo endpoint = dp_matrix(initial_path.size() + min_frames-1, 
      frames-1);
  endpoint.parent = initial_path.size() + min_frames-1;
  for(unsigned int s = endpoint.parent; s < states; s+=min_frames)
    if( dp_matrix(s, frames-1).score > endpoint.score)
    {
      endpoint.parent = s;
      endpoint.score = dp_matrix(s, frames-1).score;
    }
  
  int last = -1;
  for(int f = frames-1; f >= 0; --f)
  {
    int index = std::floor(static_cast<double>(endpoint.parent) / 
          static_cast<double>(min_frames));
    if( index < static_cast<int>(initial_path.size()) )
      index = initial_path[index];
    else
      index -= initial_path.size();

    //std::cout<<endpoint.parent<<" "<<index<<std::endl;
    if(index != last)
      ret.push_back(index);
    endpoint = dp_matrix(endpoint.parent, f);
    last = index;
  }
  std::reverse(ret.begin(), ret.end());
  return ret;
}

double GetStateScore(const utilities::Matrix<double> &pgram, 
    const std::vector<int> &initial_path, int min_frames, int state, int frame)
{
  int true_state = std::floor(static_cast<double>(state) / 
      static_cast<double>(min_frames));
  if( true_state < static_cast<int>(initial_path.size()) )
    true_state = initial_path[true_state];
  else
    true_state -= initial_path.size();

  return pgram(true_state, frame);
}

double GetTransitionScore(const utilities::Matrix<double> &transition,
    const std::vector<int> &initial_path, int min_frames, int from_state,
    int to_state)
{
  int true_from_state = std::floor(static_cast<double>(from_state) / 
      static_cast<double>(min_frames));
  int from_substate = from_state % min_frames;
  if( true_from_state < static_cast<int>(initial_path.size()))
    true_from_state = initial_path[true_from_state];
  else
    true_from_state -= initial_path.size();
  
  int true_to_state = std::floor(static_cast<double>(to_state) / 
      static_cast<double>(min_frames));
  int to_substate = to_state % min_frames;
  if( true_to_state < static_cast<int>(initial_path.size()))
    true_to_state = initial_path[true_to_state];
  else
    true_to_state -= initial_path.size();

  // First handle transition logic if we are inside the restricted path
  if( true_to_state < static_cast<int>(initial_path.size()))
  {
    if( (to_state == from_state) || (to_state == (from_state+1)) )
      return transition(true_from_state, true_to_state);
    return std::log(0);
  }
  
  if(to_substate == 0) // First substate
  {
    if(to_state == from_state) // Self loop
      return transition(true_from_state, true_to_state);
    else if(from_substate == (min_frames - 1))           // Final substate from
      return transition(true_from_state, true_to_state); // another state.
  }
  else if ((to_substate > 0) && (to_substate < (min_frames)))
  {                                     // Transition is to an interior substate
    if(to_state == (from_state + 1) ||  // so must come from immediately 
        to_state == from_state)         // preceeding substate or same substate.
      return transition(true_from_state, true_to_state);
  }
  // If we have fallen this far, the transition must not be valid.
  return std::log(0);
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
