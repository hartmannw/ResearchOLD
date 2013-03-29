// William Hartmann (hartmannw@gmail.com)
// This is free and unencumbered software released into the public domain.
// See the UNLICENSE file for more information.

#include "HmmDecoder.h"

namespace statistics
{

void HmmDecoder::Initialize(HmmSet *htk, utilities::Matrix<double> *lm)
{
  htk_ = htk;
  lm_ = lm;
  forward_prob.clear();
  back_prob.clear();
}

bool HmmDecoder::FeaturesStateScores(const utilities::Matrix<double> &data, 
    utilities::Matrix<double> &scores)
{
  unsigned int frames = data.NumCols();
  unsigned int states = htk_->StateCount();
  scores.Initialize(states, frames, zero_log_);
  for(unsigned int f = 0; f < frames; ++f)
  {
    std::vector<double> frame_data = data.GetCol(f);
    for(unsigned int s = 0; s < states; ++s)
      scores(s,f) = std::min(htk_->StateLogLikelihood(s, frame_data), 
          zero_log_);
  }
  return true;
}

ViterbiInfo HmmDecoder::SetMaximumScore(double a, double b, int a_index, 
    int b_index)
{
  ViterbiInfo point;
  if( a > b )
  {
    point.parent = a_index;
    point.score = a;
  }
  else
  {
    point.parent = b_index;
    point.score = b;
  }
  return point;
}

bool HmmDecoder::ForwardPass(const utilities::Matrix<double> &scores)
{
  // Initialize the forward pass structure
  unsigned int frames = scores.NumCols();
  forward_prob.resize( htk_->HmmCount() );
  std::vector<double> exit_prob(htk_->HmmCount(), min_log_);
  std::vector<double> next_prob = exit_prob;
  ViterbiInfo default_vi;
  default_vi.parent = 1000000;
  default_vi.score = min_log_;
  for(unsigned int i = 0; i < forward_prob.size(); ++i)
  {
    forward_prob[i].Initialize(htk_->Hmm(i).NumberOfStates(), frames, 
        default_vi);
    forward_prob[i](0,0).score = scores(htk_->Hmm(i).state(0), 0);
  }

  for(unsigned int f = 1; f < frames; ++f)
  {
    // Handle all movement internal to each HMM
    for(unsigned int h = 0; h < forward_prob.size(); ++h)
    {
      // First state can only self loop
      HiddenMarkovModel hmm = htk_->Hmm(h);
      ViterbiInfo point;
      point.parent = 0;
      point.score = forward_prob[h](0,f-1).score + scores(hmm.state(0), f) + 
        SafeLog(hmm.transition(1,1));
      forward_prob[h](0,f) = point;
      for(unsigned int s = 1; s < hmm.NumberOfStates(); ++s)
        forward_prob[h](s,f) = SetMaximumScore((forward_prob[h](s-1,f-1).score + SafeLog(hmm.transition(s,s+1)) + scores(hmm.state(s), f)),
            (forward_prob[h](s,f-1).score + SafeLog(hmm.transition(s+1,s+1)) + scores(hmm.state(s), f)),
            s-1, s);

      // Check if initial state could have come from anywhere else
      for(unsigned int p = 0; p < exit_prob.size(); ++p)
        forward_prob[h](0,f) = SetMaximumScore(forward_prob[h](0,f).score,
            (exit_prob[p] + scores(hmm.state(0), f) + lm_->At(h,p) ),
            forward_prob[h](0,f).parent, (static_cast<int>(p) + 1) * -1);

      // Update the next_prob
      unsigned int last_state = std::max(0, 
          static_cast<int>(hmm.NumberOfStates()) - 1);
      next_prob[h] = forward_prob[h](last_state,f).score + 
        SafeLog(hmm.transition(last_state+1, last_state+2));
    }
    exit_prob = next_prob;

  }

  return true;
}

std::vector<int> HmmDecoder::BestForwardPath()
{
  std::vector<int> ret;
  if( forward_prob.size() < 1)
    return ret;

  unsigned int frames = forward_prob[0].NumCols();
  ViterbiInfo last_point = forward_prob[0](forward_prob[0].NumRows()-1, 
      forward_prob[0].NumCols()-1);
  last_point.parent = 0;
  for(unsigned int h = 1; h < forward_prob.size(); ++h)
    last_point = SetMaximumScore( last_point.score, 
        forward_prob[h](forward_prob[h].NumRows()-1, 
        forward_prob[h].NumCols()-1).score, last_point.parent, h);

  ret.push_back(last_point.parent);
  double hindex = last_point.parent;
  last_point.parent = forward_prob[hindex]( forward_prob[hindex].NumRows()-1, 
      forward_prob[hindex].NumCols()-1).parent;
  
  for(int f = frames-1; f > 0; --f)
  {
    if(last_point.parent < 0)
    {
      hindex = (last_point.parent + 1) * -1;
      ret.push_back(hindex);
      last_point = forward_prob[hindex](forward_prob[hindex].NumRows()-1, f);
    }
    else
    {
      last_point = forward_prob[hindex](last_point.parent, f);
    }
  }
  std::reverse(ret.begin(), ret.end());
  return ret;
}

std::vector<int> HmmDecoder::ViterbiPath(
    const utilities::Matrix<double> &scores)
{
  ForwardPass(scores);
  return BestForwardPath();
}

double HmmDecoder::ForceAlignScore( const utilities::Matrix<double> &scores,
    const std::vector<int> &sequence )
{
  std::vector<double> ret = ForceAlignExits(scores, sequence);
  return ret[ret.size() - 1];
}

std::vector<double> HmmDecoder::ForceAlignExits( 
    const utilities::Matrix<double> &scores, 
    const std::vector<int> &sequence )
{
  std::vector<double> ret;
  std::vector<utilities::Matrix<ViterbiInfo> > dp_matrix;
  unsigned int frames = scores.NumCols();
  unsigned int labels = sequence.size();

  ret.resize(frames, min_log_);
  dp_matrix.resize(labels);

  // Initialize dp_matrix
  ViterbiInfo default_point;
  default_point.parent = 1000000;
  default_point.score = min_log_;
  for(unsigned int i = 0; i < labels; ++i)
    dp_matrix[i].Initialize(htk_->Hmm(sequence[i]).NumberOfStates(), frames, 
        default_point);
  dp_matrix[0](0,0).score = scores(htk_->Hmm(sequence[0]).state(0), 0);

  for(unsigned int f = 1; f < frames; ++f)
  {
    for(unsigned int h = 0; h < labels; ++h)
    {
      // First state can only self loop
      HiddenMarkovModel hmm = htk_->Hmm(sequence[h]);
      ViterbiInfo point;
      point.parent = 0;
      point.score = dp_matrix[h](0,f-1).score + scores(hmm.state(0), f) + 
        SafeLog(hmm.transition(1,1));
      dp_matrix[h](0,f) = point;

      if(h > 0) // Could have come from previous Hmm
      {
        HiddenMarkovModel lasthmm = htk_->Hmm(sequence[h-1]);
        unsigned int laststate = lasthmm.NumberOfStates()-1;
        dp_matrix[h](0,f) = SetMaximumScore(dp_matrix[h](0,f).score, 
            (dp_matrix[h-1](laststate,f-1).score + SafeLog(lasthmm.transition(laststate+1, laststate+2)) + lm_->At(sequence[h-1], sequence[h])),
            0, -1);
      }

      for(unsigned int s = 1; s < hmm.NumberOfStates(); ++s)
        dp_matrix[h](s,f) = SetMaximumScore((dp_matrix[h](s-1,f-1).score + SafeLog(hmm.transition(s,s+1)) + scores(hmm.state(s), f)),
            (dp_matrix[h](s,f-1).score + SafeLog(hmm.transition(s+1,s+1)) + scores(hmm.state(s), f)),
            s-1, s);
    }
    ret[f] = dp_matrix[labels-1](
        htk_->Hmm(sequence[labels-1]).NumberOfStates()-1, f).score;
  }

  return ret;
}

}
