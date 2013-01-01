#ifndef ASTRONOMY_ORBITALBODY_
#define ASTRONOMY_ORBITALBODY_
#include <iostream>
#include <cmath>
#include <string>
#include <sstream>

namespace astronomy
{

#define KM_PER_AU 149597870.691
#define GRAVITY 7.17633288458e-17 // defined in terms of AU / 
                                  // (10^24kg * minutes^2)

class OrbitalBody
{
 public:
  OrbitalBody(){}
  ~OrbitalBody(){}
  void Initialize(std::string name, double mass, double ascension_one, 
      double declination_one, double delta_one, double ascension_two, 
      double declination_two, double delta_two);
  double mass() const { return mass_;}
  double x() const { return x_; }
  double y() const { return y_; }
  double z() const { return z_; }
  double dx() const { return dx_; }
  double dy() const { return dy_; }
  double dz() const { return dz_; }
  void UpdatePosition();
  void UpdateVelocity(double ax, double ay, double az);
  std::string GetInformation();

  double DistanceFrom(const OrbitalBody &body);

 private:
  std::string name_;
  double mass_;
  double x_, y_, z_;
  double dx_, dy_, dz_;
};

}

#endif
