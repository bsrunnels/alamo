#ifndef MODEL_INTERFACE_GB_READ_H
#define MODEL_INTERFACE_GB_READ_H

#include <iostream>
#include <fstream>
#include <vector>
#include <typeinfo> 

#include "AMReX.H"
#include "GB.H"

#define PI 3.14159265 

namespace Model
{
namespace Interface
{
namespace GB
{
class Read : public GB
{
public:
  Read()
  {
    
  };
  Read(std::string filename)
  {
    std::ifstream ifs{filename};
    amrex::Real theta, W, res, thetamin, thetamax;

    if (!ifs){ Util::Abort(INFO, "Unable to open the file."); };

    while(ifs >> theta >> W){                   // reads data.
      thetas.push_back(fmod(theta, 360) * PI / 180);
      Ws.push_back(W);
    };

    size = thetas.size();
    res = thetas[1] - thetas[0];
    thetamin = thetas[0];
    thetamax = thetas[0];

    for (i = 1; i < size; i++){                 // checks for full range of thetas.
      if (thetas[i] < thetamin){thetamin = thetas[i]; };
      if (thetas[i] > thetamax){thetamax = thetas[i]; };
    };
    if ((thetamax - thetamin + res) < 2.0*PI){
      Util::Abort(INFO, "Input file doesn't cover full range of thetas (0 to 2PI).");
    };

    DWs.push_back((Ws[1] - Ws[0]) / res);       // takes derivatives.
    for (i = 1; i < size; i++){
      DWs.push_back((Ws[(i+1) % size] - Ws[i]) / res);
      DDWs.push_back((DWs[i] - DWs[i-1]) / res);
    };
    DDWs.push_back((DWs[0] - DWs[size-1]) / res);
  };

  amrex::Real W(amrex::Real theta)
  {
    amrex::Real Wumbo{interpolate(theta, Ws)};
    return Wumbo;
  };
  amrex::Real DW(amrex::Real theta)
  {
    amrex::Real DWumbo{interpolate(theta, DWs)};
    return DWumbo;
  };
  amrex::Real DDW(amrex::Real theta)
  {
    amrex::Real DDWumbo{interpolate(theta, DDWs)};
    return DDWumbo;
  };
 
private:
  int i, size;
  std::vector<amrex::Real> thetas, Ws, DWs, DDWs;

  double interpolate(amrex::Real theta, std::vector<amrex::Real> data)  //interpolates to get data between the points.
  {
    int revs, lo_index{0}, hi_index{size - 1};

    revs = fabs(theta) / (2*PI);
    if (theta > 2*PI){ theta -= (revs*2*PI); };
    if (theta < 0){ theta += ((revs+1)*2*PI); };

    for (i = 0; i < size; i++){
      if ((thetas[i] > thetas[lo_index]) && (thetas[i] < theta)) { lo_index = i; };
      if ((thetas[i] < thetas[hi_index]) && (thetas[i] > theta)) { hi_index = i; };
    };

    if (hi_index == lo_index){
      return data[lo_index];
    }
    return ((theta - thetas[lo_index]) / (thetas[hi_index]- thetas[lo_index])) * (data[hi_index] - data[lo_index]) + data[lo_index];
  };
};
}
}
}

#endif