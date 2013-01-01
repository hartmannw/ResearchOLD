#ifndef ASTRONOMY_GRAVITATIONALSYSTEM_H_
#define ASTRONOMY_GRAVITATIONALSYSTEM_H_

#include<vector>
#include<sstream>
#include "OrbitalBody.h"

namespace astronomy
{

class GravitationalSystem
{
 public:
  GravitationalSystem(){};
  ~GravitationalSystem(){};

  void AddBody(OrbitalBody body){bodies_.push_back(body);}
  void AdvanceStep();
  std::string GetInformation();
  OrbitalBody body(unsigned int i){return bodies_[i];}
  std::vector<OrbitalBody> bodies(){ return bodies_;}

 private:
  std::vector<OrbitalBody> bodies_;
};

}

#endif
