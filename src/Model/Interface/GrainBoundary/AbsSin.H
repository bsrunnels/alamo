#ifndef MODEL_INTERFACE_GRAINBOUNDARY_ABSSIN_H
#define MODEL_INTERFACE_GRAINBOUNDARY_ABSSIN_H

#include <iostream>
#include <fstream>

#include "AMReX.H"
#include "GrainBoundary.H"

#define PI 3.14159265 

namespace Model
{
namespace Interface
{
namespace GrainBoundary
{
class AbsSin : public GrainBoundary
{
public:
  AbsSin(amrex::Real _theta0, amrex::Real _sigma0, amrex::Real _sigma1, amrex::Real _frequency)
  {       
    theta0 = _theta0;
    sigma0 = _sigma0;
    sigma1 = _sigma1;
    frequency = _frequency;
  };
  amrex::Real W(amrex::Real theta)
  {
    //sigma(theta)=sigma0+sigma1*fabs(sin(n*(theta-theta0)))
    //n=2:  
    return sigma0+sigma1*fabs(sin( (frequency / 2) *(theta-theta0)));  
    
  };
  amrex::Real DW(amrex::Real theta)
  {
    //sigma'(theta)=(n*sigma1*sin(n*(theta-theta0))*cos(n*(theta-theta0)))/fabs(sin(n*(theta-theta0)))
    //n=2:
    return ( (frequency / 2) *sigma1*sin( (frequency / 2) *(theta-theta0))*cos( (frequency / 2) *(theta-theta0)))/fabs(sin( (frequency / 2) *(theta-theta0)));
  };
  amrex::Real DDW(amrex::Real theta)
  {
    //sigma''(theta)=-n^2*sigma1*fabs(sin(n*(theta-theta0)))
    //n=2:
    return -1 * (frequency * frequency / 4) *sigma1*fabs(sin( (frequency / 2) *(theta-theta0)));
  };
 
private:
  amrex::Real theta0, sigma0, sigma1, frequency;
};
}
}
}

#endif
