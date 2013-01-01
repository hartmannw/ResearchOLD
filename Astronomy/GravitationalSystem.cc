#include "GravitationalSystem.h"

namespace astronomy
{

void GravitationalSystem::AdvanceStep()
{
  std::vector<double> ax, ay, az;
  for(unsigned int i = 0; i < bodies_.size(); ++i)
  {
    ax.push_back(0);
    ay.push_back(0);
    az.push_back(0);
    for(unsigned int j=0; j < bodies_.size(); ++j)
    {
      if(i != j)
      {
        double distance = std::sqrt( std::pow(bodies_[i].x() - bodies_[j].x(),2)
            + std::pow(bodies_[i].y() - bodies_[j].y(), 2)
            + std::pow(bodies_[i].z() - bodies_[j].z(), 2));
        ax[i] += - (GRAVITY * (bodies_[j].mass() * 
            (bodies_[i].x() - bodies_[j].x()))) / pow(distance, 3);
        ay[i] += - (GRAVITY * (bodies_[j].mass() * 
            (bodies_[i].y() - bodies_[j].y()))) / pow(distance, 3);
        az[i] += - (GRAVITY * (bodies_[j].mass() * 
            (bodies_[i].z() - bodies_[j].z()))) / pow(distance, 3);
      }
    }
  }
  for(unsigned int i = 0; i < bodies_.size(); i++)
  {
    bodies_[i].UpdatePosition();
    bodies_[i].UpdateVelocity(ax[i], ay[i], az[i]);
  }
  return;
}

std::string GravitationalSystem::GetInformation()
{
  std::ostringstream ss;
  for(unsigned int i = 0; i < bodies_.size(); i++)
  {
    ss << "Body " << (i+1) <<": "<<bodies_[i].GetInformation()<<"\n";
  }
  return ss.str();
}

}
