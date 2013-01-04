#ifndef STATISTICS_HIDDENMARKOVMODEL_
#define STATISTICS_HIDDENMARKOVMODEL_

#include <string>
#include <vector>

namespace statistics
{

class HiddenMarkovModel
{
 private:
  std::vector<unsigned int> states_;
  std::vector<std::vector<double> > transition_matrix_;
  std::string name_;

 public:
  HiddenMarkovModel(){}
  ~HiddenMarkovModel(){}

  void SetName(std::string name){name_ = name;}
  void AddState(unsigned int state){states_.push_back(state);}
  void SetStates(std::vector<unsigned int> states){states_ = states;}
  void SetTransitions(std::vector<std::vector<double> > transition){
      transition_matrix_ = transition;}

  std::string name(){return name_;}
  unsigned int state(unsigned int state_index){return states_[state_index];}
  std::vector<unsigned int> states(){return states_;}
};

}

#endif
