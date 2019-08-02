#ifndef MODEL_INTERFACE_GB_ABSSIN_H
#define MODEL_INTERFACE_GB_ABSSIN_H

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
class AbsSin : public GB
{
public:
  AbsSin()
  {
  };
  AbsSin(amrex::Real _theta0, amrex::Real _sigma0, amrex::Real _sigma1, amrex::Real _frequency)
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
    //sigma(theta)=sigma0+sigma1*fabs(sin(n*(theta-theta0)))
    //n=2: 

    if (std::isnan(sigma0+sigma1*fabs(sin((frequency / 2)*(theta-theta0))))) Util::Abort(INFO, "theta: ", theta);
    return sigma0+sigma1*fabs(sin((frequency / 2)*(theta-theta0)));  
    
  };
  amrex::Real DW(amrex::Real theta)
  {
    //sigma'(theta)=(n*sigma1*sin(n*(theta-theta0))*cos(n*(theta-theta0)))/fabs(sin(n*(theta-theta0)))
    //n=2:
    if (std::isnan(((frequency / 2)*sigma1*sin((frequency / 2)*(theta-theta0))*cos((frequency / 2)*(theta-theta0)))/fabs(sin((frequency / 2)*(theta-theta0)))) Util::Abort(INFO, "theta: ", theta);
    return ((frequency / 2)*sigma1*sin((frequency / 2)*(theta-theta0))*cos((frequency / 2)*(theta-theta0)))/fabs(sin((frequency / 2)*(theta-theta0)));
  };
  amrex::Real DDW(amrex::Real theta)
  {
    //sigma''(theta)=-n^2*sigma1*fabs(sin(n*(theta-theta0)))
    //n=2:
    if (std::isnan(-1*(frequency*frequency / 4)*sigma1*fabs(sin((frequency / 2)*(theta-theta0)))) Util::Abort(INFO, "theta: ", theta);
    return -1*(frequency*frequency / 4)*sigma1*fabs(sin((frequency / 2)*(theta-theta0)));
  };
 
private:
  amrex::Real theta0 = NAN, sigma0 = NAN, sigma1 = NAN, frequency = NAN;
};
}
}
}

#endif
