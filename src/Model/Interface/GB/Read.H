#ifndef MODEL_INTERFACE_GB_READ_H
#define MODEL_INTERFACE_GB_READ_H

#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>

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
  Read(std::string filename, bool translate, amrex::Real _sigma0, amrex::Real _sigma1, amrex::Real _FWHM)
  {
    std::ifstream ifs{filename};
    amrex::Real theta, W;
    int j;

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

    if (translate){                           // smooths, scales, and shifts data
      bool smooth{true};
      Define(_sigma0, _sigma1, _FWHM);

      amrex::Real sigma = FWHM / sqrt(8.0 * log(2.0));

      if(sigma < res){Util::Abort(INFO, "std. dev. of kernel must be greater than resolution, try increasing FWHM"); };     // 

      std::vector<amrex::Real> kernal_integrals(size, 0);
      std::vector<amrex::Real> smoothed_Ws(size, 0);
      
      for (i = 0; i < size; i++){
        std::vector<amrex::Real> theta_kernels(size, 0);

        for (j = 0; j < size; j++){           // smooths data
          theta_kernels[j] = exp(-0.5 * pow((((thetas[j] - thetas[i]) * 180 / PI) / sigma), 2));
          kernal_integrals[i] += theta_kernels[j];
        };

        for (j = 0; j < size; j++){
          theta_kernels[j] /= kernal_integrals[i];
          smoothed_Ws[i] += theta_kernels[j] * Ws[j];
        };
      };
      
      amrex::Real sWmin = smoothed_Ws[0];    // scales and shifts smoothed data
      amrex::Real sWmax = smoothed_Ws[0];
      
      for (i = 0; i < size; i++){
        if (smoothed_Ws[i] < sWmin){sWmin = smoothed_Ws[i]; };
        if (smoothed_Ws[i] > sWmax){sWmax = smoothed_Ws[i]; };
      };
      
      amrex::Real scale = sigma1 / (sWmax - sWmin);
      amrex::Real shift = sigma0 - (sWmin * scale);
      
      for (i = 0; i < size; i++){
        Ws[i] = smoothed_Ws[i] * scale + shift;
      };

      modfile(filename, thetas, Ws);
    };

    DWs.push_back((Ws[1] - Ws[0]) / res);   // takes derivatives.
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

  void Define(amrex::Real _sigma0, amrex::Real _sigma1, amrex::Real _FWHM)
	{
		sigma0 = _sigma0;
		sigma1 = _sigma1;
		FWHM = _FWHM;
    if (std::isnan(sigma0)){sigma0 = 0.15; };       // default values which generally work.
    if (std::isnan(sigma1)){sigma1 = 0.175; };
    if (std::isnan(FWHM)){FWHM = 5.0; };
	};
 
private:
  int i, size;
  double res;
  std::vector<amrex::Real> thetas, Ws, DWs, DDWs;
  amrex::Real thetamin, thetamax, sigma0, sigma1, FWHM;

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

  int modfile(std::string filename, std::vector<amrex::Real> data1, std::vector<amrex::Real> data2)   // writes modified data to a file in the same directory as the input file
  {
    for (i = 0; i < 4; i++){ filename.pop_back(); };
    filename.append("-translated.txt");

    std::ofstream ofs{filename};
    
    for (i = 0; i < size; i++){
      ofs << data1[i] * 180 / PI << " " << data2[i] << "\n";
    };

    ofs.close();
    return 0;
  };
  
};
}
}
}
#endif
