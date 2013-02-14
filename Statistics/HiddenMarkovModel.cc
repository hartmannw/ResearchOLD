// William Hartmann (hartmannw@gmail.com)
// This is free and unencumbered software released into the public domain.
// See the UNLICENSE file for more information.

#include "HiddenMarkovModel.h"

namespace statistics
{

double HiddenMarkovModel::ExpectedDuration()
{
  double ret = 0;
  for(unsigned int i = 1; i <  transition_matrix_.NumRows()-1; ++i)
    ret += 1 + (1 / (1 - transition_matrix_(i,i)));
  return ret;
}

}
