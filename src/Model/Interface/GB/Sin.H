#ifndef MODEL_INTERFACE_GB_SIN_H
#define MODEL_INTERFACE_GB_SIN_H

#include <iostream>
#include <fstream>

#include "AMReX.H"
#include "GB.H"
#include "Set/Set.H"
#include "Util/Util.H"

#define PI 3.14159265 

namespace Model
{
namespace Interface
{
namespace GB
{
class Sin : public GB
{
public:
	Sin()
	{
	};
	Sin(amrex::Real _theta0, amrex::Real _sigma0, amrex::Real _sigma1, amrex::Real _frequency)
	{
		Define(_theta0,_sigma0,_sigma1,_frequency);
	};
	void Define(amrex::Real _theta0, amrex::Real _sigma0, amrex::Real _sigma1, amrex::Real _frequency)
	{
		theta0 = _theta0;
		sigma0 = _sigma0;
		sigma1 = _sigma1;
		frequency = _frequency;
	};
	void Randomize()
	{
		theta0 = Util::Random()*Set::Constant::Pi;
		sigma0 = Util::Random();
		sigma1 = Util::Random();
		frequency = 4;
	};
	amrex::Real W(amrex::Real theta)
	{
		//if(std::isnan(sigma0 + 0.5*sigma1*(1.0 - cos(frequency*(theta-theta0))))) Util::Abort(INFO, "sigma0: ", sigma0, "\n sigma1: ", sigma1, "\n frequency: ", frequency, "\n theta0: ", theta0, "\n theta: ", theta);
		return sigma0 + 0.5*sigma1*(1.0 - cos(frequency*(theta-theta0)));
	};
	amrex::Real DW(amrex::Real theta)
	{
		return (frequency / 2) *sigma1*sin(frequency*(theta-theta0));
	};
	amrex::Real DDW(amrex::Real theta)
	{
		return (frequency*frequency / 2)*sigma1*cos(frequency*(theta-theta0));
	};
  
private:
	amrex::Real theta0 = NAN, sigma0 = NAN, sigma1 = NAN, frequency = NAN;
};
}
}
}
#endif
