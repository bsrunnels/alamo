#ifndef MODEL_INTERFACE_GB_SIN_H
#define MODEL_INTERFACE_GB_SIN_H

#include <iostream>
#include <fstream>

#include "AMReX.H"
#include "GB.H"
#include "Set/Set.H"
#include "Util/Util.H"

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
		return sigma0 + 0.5*sigma1*(1.0 - cos(frequency*(theta-theta0)));
	};
	amrex::Real DW(amrex::Real theta)
	{
		return (0.5 * frequency) *sigma1*sin(frequency*(theta-theta0));
	};
	amrex::Real DDW(amrex::Real theta)
	{
		return (0.5 * frequency*frequency)*sigma1*cos(frequency*(theta-theta0));
	};
  
private:
	amrex::Real theta0 = NAN, sigma0 = NAN, sigma1 = NAN, frequency = NAN;
};
}
}
}
#endif
