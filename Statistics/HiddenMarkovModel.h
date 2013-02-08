// William Hartmann (hartmannw@gmail.com)
// This is free and unencumbered software released into the public domain.
// See the UNLICENSE file for more information.

#ifndef STATISTICS_HIDDENMARKOVMODEL_
#define STATISTICS_HIDDENMARKOVMODEL_

#include <string>
#include <vector>
#include "Matrix.h"

namespace statistics
{

// This is simple a container class for a hidden markov model. It is assumed the
// distributions for the states are stored elsewhere. Only an ID for the state
// is maintained here. The transition matrix is fully stored though. The purpose
// is to organize the data, not train or evaluate the HMM.
class HiddenMarkovModel
{
 private:
  std::vector<unsigned int> states_; // IDs for the state that hopefully point
                                     // to the actually distribution stored 
                                     // elsewhere.
  utilities::Matrix<double> transition_matrix_;
  std::string name_;

 public:
  HiddenMarkovModel(){}
  ~HiddenMarkovModel(){}

  // Standard accessor and mutator functions.
  void SetName(std::string name){name_ = name;}
  void AddState(unsigned int state){states_.push_back(state);}
  void SetStates(std::vector<unsigned int> states){states_ = states;}
  void SetTransitions(utilities::Matrix<double> transition){
      transition_matrix_ = transition;}

  std::string name(){return name_;}
  unsigned int state(unsigned int state_index){return states_[state_index];}
  std::vector<unsigned int> states(){return states_;}
  unsigned int NumberOfStates() { return states_.size(); }
  double transition(unsigned int r, unsigned int c) {
      return transition_matrix_(r,c);}
};

}

#endif
