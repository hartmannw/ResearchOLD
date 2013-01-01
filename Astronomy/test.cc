#include "OrbitalBody.h"
#include "GravitationalSystem.h"
#define _USE_MATH_DEFINES
#include <cmath>
#include <iostream>

double ConvertAscensionToRadians(double hours, double minutes, double seconds)
{
  double ret;
  ret = (hours / 24) + (minutes / 1440) + (seconds / 86400);
  ret = ret * (M_PI * 2);
  return ret;
}

double ConvertDeclinationToRadians(double degrees, double minutes, 
    double seconds)
{
  double ret;
  ret = (degrees / 360) + (minutes / 21600) + (seconds / 1296000);
  ret = ret * (M_PI * 2);
  return ret;
}

int main()
{
  astronomy::OrbitalBody earth, sun, jupiter, moon;
  astronomy::GravitationalSystem solar_system;
  double ascension_one = ConvertAscensionToRadians(2, 46, 22.74);
  double ascension_two = ConvertAscensionToRadians(2, 46, 32.70);
  double declination_one = ConvertDeclinationToRadians(16, 3, 19.0);
  double declination_two = ConvertDeclinationToRadians(16, 4, 3.7);
  double mass = 5.9736;
  double sun_mass = 1989100;
  double delta_one = 0.99115495427797;
  double delta_two = 0.99114452058751;
  earth.Initialize("Earth", mass, ascension_one, declination_one, delta_one, 
      ascension_two, declination_two, delta_two);
  sun.Initialize("Sun", sun_mass, 0, 0, 0, 0, 0, 0);
  ascension_one = ConvertAscensionToRadians(2, 16, 3.75);
  ascension_two = ConvertAscensionToRadians(2, 16, 4.63);
  declination_one = ConvertDeclinationToRadians(12, 23, 8.5);
  declination_two = ConvertDeclinationToRadians(12, 23, 13.1);
  double jupiter_mass = 1898.13;
  delta_one = 4.96442858953864;
  delta_two = 4.96443425534284;
  jupiter.Initialize("Jupiter", jupiter_mass, ascension_one, declination_one, delta_one, 
      ascension_two, declination_two, delta_two);
  ascension_one = ConvertAscensionToRadians(2, 45, 56.98);
  ascension_two = ConvertAscensionToRadians(2, 46, 7.16);
  declination_one = ConvertDeclinationToRadians(16, 2, 14.1);
  declination_two = ConvertDeclinationToRadians(16, 2, 59.8);
  double moon_mass = .07349;
  delta_one = 0.99316716681236;
  delta_two = 0.99317120335859;
  moon.Initialize("Moon", moon_mass, ascension_one, declination_one, delta_one, 
      ascension_two, declination_two, delta_two);

  solar_system.AddBody(sun);
  solar_system.AddBody(earth);
  solar_system.AddBody(jupiter);
  solar_system.AddBody(moon);
  std::cout<< solar_system.GetInformation();
  bool greater = true;
  double initialx = solar_system.body(1).x();
  double initialy = solar_system.body(1).y();
  for(int i = 0; i < 525949 * 10; ++i)
  {
    solar_system.AdvanceStep();
    double curx = solar_system.body(1).x();
    double cury = solar_system.body(1).y();
    if(greater == false)
    {
      if( (cury / curx) > (initialy / initialx) )
      {
        greater = true;
        std::cout<<(static_cast<double>(i) / 525949)<<std::endl;
        std::cout<< solar_system.GetInformation()<<std::endl;
      }
    }
    else
    {
      if( (cury / curx) < (initialy / initialx) )
        greater = false;
    }
  }
  std::cout<<solar_system.GetInformation()<<std::endl;
  /*
  for(int y=0; y < 10; y++)
  {
    for(int i =0; i < 525949; i++)
    {
      solar_system.AdvanceStep();
    }
    std::cout<< solar_system.GetInformation();
  }
  */
  return 0;
}
