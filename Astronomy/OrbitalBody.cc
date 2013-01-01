#include "OrbitalBody.h"

namespace astronomy
{

void OrbitalBody::Initialize(std::string name, double mass, 
    double ascension_one, double declination_one, double delta_one, 
    double ascension_two, double declination_two, double delta_two)
{
  name_ = name;
  mass_ = mass;
  z_ =sqrt((delta_one * delta_one) - pow( delta_one * cos(declination_one), 2));
  if(sin(declination_one) < 0)
    z_ = -z_;
  double A = sqrt( pow(delta_one,2) - pow(z_, 2) );
  y_ =sqrt((A * A) - pow( A * cos(ascension_one), 2));
  if(sin(ascension_one) < 0)
    y_ = -y_;
  x_ = sqrt( pow(A,2) - pow(y_, 2));
  if(cos(ascension_one) < 0)
    x_ = -x_;

  // Repeat for the second point.
  double z2 = sqrt((delta_two * delta_two) - 
      pow( delta_two * cos(declination_two), 2));
  if(sin(declination_two) < 0)
    z2 = -z2;
  A = sqrt( pow(delta_two,2) - pow(z2, 2) );
  double y2 =sqrt((A * A) - pow( A * cos(ascension_two), 2));
  double x2 = sqrt( pow(A,2) - pow(y2, 2));
  if(sin(ascension_two) < 0)
    y2 = -y2;
  if(cos(ascension_two) < 0)
    x2 = -x2;
  dx_ = x2 - x_; dx_ = dx_ / 60;
  dy_ = y2 - y_; dy_ = dy_ / 60;
  dz_ = z2 - z_; dz_ = dz_ / 60;
  return;
}

void OrbitalBody::UpdatePosition()
{
  x_ += dx_;
  y_ += dy_;
  z_ += dz_;
}

void OrbitalBody::UpdateVelocity(double ax, double ay, double az)
{
  dx_ += ax;
  dy_ += ay;
  dz_ += az;
}

std::string OrbitalBody::GetInformation()
{
  std::ostringstream ss;
  ss<<name_<<" At ("<<x_<<","<<y_<<","<<z_<<"), "<<
      std::sqrt(pow(x_,2)+pow(y_,2)+pow(z_,2))<<" with velocity ("<<
      dx_<<","<<dy_<<","<<dz_<<"), "<< 
      std::sqrt(pow(dx_,2)+pow(dy_,2)+pow(dz_,2))*KM_PER_AU;
  return ss.str();
}

double OrbitalBody::DistanceFrom(const OrbitalBody &body)
{
  double ret = 0;
  ret += std::pow(x_ - body.x(),2);
  ret += std::pow(y_ - body.y(),2);
  ret += std::pow(z_ - body.z(),2);
  return std::sqrt(ret);
}

}
