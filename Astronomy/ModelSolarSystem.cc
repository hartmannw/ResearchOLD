#include "OrbitalBody.h"
#include "GravitationalSystem.h"
#include "StringFunctions.h"
#define _USE_MATH_DEFINES
#include <cmath>
#include <iostream>
#include <vector>
#include <string>
#include <fstream>

double ConvertAscensionToRadians(double hours, double minutes, double seconds);
double ConvertDeclinationToRadians(double degrees, double minutes,
    double seconds);

void TestEarthOrbit( astronomy::GravitationalSystem s, double max_years)
{
  astronomy::OrbitalBody original_earth;
  original_earth = s.body(3);
  for(int i = 0; i < (525949 * max_years); ++i)
  {
    if( (i % 525949) == 0)
    {
      std::cout<<"Year "<<(i / 525949.0)<<std::endl;
      //std::cout<<s.GetInformation()<<std::endl;
      std::cout<<original_earth.DistanceFrom(s.body(3))<<std::endl;
    }
    s.AdvanceStep();
  }
  std::cout<<original_earth.GetInformation()<<std::endl;
  std::cout<<s.body(3).GetInformation()<<std::endl;
}

double AverageSystemDistance(astronomy::GravitationalSystem gs_start,
    astronomy::GravitationalSystem gs_end)
{
  double ret=0;
  int total = 0;
  for(int i = 1; i < gs_start.bodies().size(); ++i)
  {
    double distance = (gs_start.body(i).DistanceFrom(gs_end.body(i))); 
        //gs_end.body(0).DistanceFrom(gs_end.body(i));
    ret += distance;
    std::cout<<distance<<" "<<
        distance / gs_end.body(0).DistanceFrom(gs_end.body(i))<< std::endl;
    total++;
  }
  return ret / total;
}

void TestSystem(astronomy::GravitationalSystem gs_start,
    astronomy::GravitationalSystem gs_end, double days)
{
  for(int i = 0; i < (1440 * days); ++i)
    gs_start.AdvanceStep();

  std::cout<<AverageSystemDistance(gs_start, gs_end)<<std::endl;
  std::cout<<gs_start.GetInformation()<<std::endl;
}

astronomy::GravitationalSystem LoadSolarSystem(std::string filename)
{
  astronomy::GravitationalSystem ret;
  std::ifstream fin;
  fin.open(filename.c_str(), std::ios::in);
  std::string line;
  while( getline(fin, line) )
  {
    std::vector<std::string> token;
    astronomy::OrbitalBody o;
    std::string name;
    double weight;
    double ascension_one, ascension_two;
    double declination_one, declination_two;
    double delta_one, delta_two;
    line = utilities::TrimString(line);
    utilities::TokenizeString(line, ' ', token);
    name = token[0];
    weight = utilities::ToNumber<double>(token[1]);
    getline(fin, line);
    token.clear();
    line = utilities::TrimString(line);
    utilities::TokenizeString(line, ' ', token);
    ascension_one = ConvertAscensionToRadians(
        utilities::ToNumber<double>(token[0]), 
        utilities::ToNumber<double>(token[1]), 
        utilities::ToNumber<double>(token[2]));
    declination_one = ConvertDeclinationToRadians(
        utilities::ToNumber<double>(token[3]), 
        utilities::ToNumber<double>(token[4]), 
        utilities::ToNumber<double>(token[5]));
    delta_one = utilities::ToNumber<double>(token[6]);
    getline(fin, line);
    token.clear();
    line = utilities::TrimString(line);
    utilities::TokenizeString(line, ' ', token);
    ascension_two = ConvertAscensionToRadians(
        utilities::ToNumber<double>(token[0]), 
        utilities::ToNumber<double>(token[1]), 
        utilities::ToNumber<double>(token[2]));
    declination_two = ConvertDeclinationToRadians(
        utilities::ToNumber<double>(token[3]), 
        utilities::ToNumber<double>(token[4]), 
        utilities::ToNumber<double>(token[5]));
    delta_two = utilities::ToNumber<double>(token[6]);
    o.Initialize(name, weight, ascension_one, declination_one, delta_one, 
        ascension_two, declination_two, delta_two);
    ret.AddBody(o);
  }
  fin.close();
  return ret;
}

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
    if(degrees < 0)
      ret = (degrees / 360) - (minutes / 21600) - (seconds / 1296000);
    else
      ret = (degrees / 360) + (minutes / 21600) + (seconds / 1296000);
    ret = ret * (M_PI * 2);
    return ret;
}

int main()
{
    astronomy::GravitationalSystem solar_system, ss10, ss2061;
    solar_system = LoadSolarSystem(std::string("planetdata02091986.txt"));
    ss10 = LoadSolarSystem(std::string("planetdata02091996.txt"));
    ss2061 = LoadSolarSystem(std::string("planetdata07282061.txt"));
    std::cout<<solar_system.GetInformation()<<std::endl;
    std::cout<<ss10.GetInformation()<<std::endl;
    std::cout<<ss2061.GetInformation()<<std::endl;

    TestSystem(solar_system, ss10, 27563);
    //TestEarthOrbit(solar_system, 10.0);
    //std::cout<<solar_system.body(3).DistanceFrom(solar_system.body(9))<<std::endl;
    return 0;
}
