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
    amrex::Real theta, W, res;

    if (!ifs){ Util::Abort(INFO, "Unable to open the file."); };

    while(ifs >> theta >> W){                   // reads data.
      thetas.push_back(fmod(theta, 360) * PI / 180);
      Ws.push_back(W);
    };

    size = thetas.size();
    res = thetas[1] - thetas[0];
    thetamin = thetas[0];
    thetamax = thetas[0];
    
    for (i = 0; i < size; i++){
      if (thetas[i] > thetamax) { thetamax = thetas[i]; };
    }

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
    return interpolate(theta, Ws);
  };
  amrex::Real DW(amrex::Real theta)
  {
    return interpolate(theta, DWs);
  };
  amrex::Real DDW(amrex::Real theta)
  {
    return interpolate(theta, DDWs);
  };
 
private:
  int i, size;
  std::vector<amrex::Real> thetas, Ws, DWs, DDWs;
  amrex:: Real thetamin, thetamax;

  double interpolate(amrex::Real theta, std::vector<amrex::Real> data)  //interpolates to get data between the points.
  {
    int revs, hi_index;
    double lo_index;

    revs = fabs(theta) / (2*PI);
    if (theta > 2*PI){ theta -= (revs*2*PI); };
    if (theta < 0){ theta += ((revs+1)*2*PI); };

    lo_index = theta / thetamax * size;
    if (fmod(lo_index, 1) == 0) { return data[lo_index]; };
    
    lo_index = int(lo_index);
    hi_index = lo_index + 1;

    if (thetas[hi_index] == thetas[lo_index]) { return data[lo_index]; };
    return ((theta - thetas[lo_index]) / (thetas[hi_index]- thetas[lo_index])) * (data[hi_index] - data[lo_index]) + data[lo_index];
  };
};
}
}
}
#endif