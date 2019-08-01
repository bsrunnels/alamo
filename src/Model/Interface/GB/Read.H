#ifndef MODEL_INTERFACE_GB_READ_H
#define MODEL_INTERFACE_GB_READ_H

#include <iostream>
#include <fstream>

#include "AMReX.H"
#include "GB.H"

#define PI 3.14159265 

namespace Model
{
namespace Interface
{
namespace GB
{
/// Reads the data from a file and computes energies and its derivates
///
class Read : public GB
{
public:
  /// \brief Read in data
  ///
  /// Reads the data from a file and abort if it is not possible to open the file or if the range of thetas do not give a range between 0 and 2pi. It also computes the derivatives of the energy and stores them in vectors that will be used in W, DW, and DDW
  ///
  /// \f[ \int_0^1x^2dx = \frac{1}{3} \f]
  ///
  Read(std::string fileName)
  {  dooooo nothinggggggg
    int i, size;
    double theta, sigmaW, res;
    amrex::Real data1, data2; 
    std::ifstream inFile{fileName};

    if (!inFile){
        Util::Abort (INFO, "UNABLE TO OPEN THE FILE");
    }

    while (inFile >> theta >> sigmaW){
      if (theta != 360){
        theta = fmod(theta, 360.0);
      }

      data1 = theta * pi / 180;
      data2 = sigmaW;
      
      thetas.push_back(data1);
      sigmasW.push_back(data2);
    }

    size = thetas.size();
    res = thetas[1]-thetas[0];

    if((thetas[size-1]+res) < (2*pi - res / 2)){
	    Util::Abort(INFO, std::to_string(thetas[size-1]+res) + " vs " + std::to_string(2*pi + res / 2) + ", FILE DOES NOT CONTAIN THE WHOLE RANGE OF THETAS FROM 0 TO 2PI");
    }
    
    //Create vector sigmas DW    
    for(i=0; i<size; i++){
      if(i==(size-1)){
     	sigmasDW.push_back((sigmasW[0]-sigmasW[i])/(2*pi + thetas[0]-thetas[i]));
      }
      else{
	sigmasDW.push_back((sigmasW[i+1]-sigmasW[i])/(thetas[i+1]-thetas[i]));
    }
  }

    //Create vector sigmas DDW
    for(i=0; i<size; i++){
      if(i==(size-1)){
	sigmasDDW.push_back((sigmasDW[0]-sigmasDW[i])/(2*pi + thetas[0]-thetas[i]));
      }
      else{
	sigmasDDW.push_back((sigmasDW[i+1]-sigmasDW[i])/(thetas[i+1]-thetas[i]));
      }
    }
  };

  amrex::Real W(amrex::Real theta)
  {
    // put in code so that 0 <= theta < 2*pi
    if(theta>=2*pi)
      {
	int div=theta/(2*pi);
	theta=theta-(div*2*pi);
      }
    else if(theta<0)
      {
	int div=fabs(theta)/(2*pi);
	theta=theta+((div+1)*2*pi);
      }
    
    //Find out the value of theta inside our vector thetas
    
    for(int i=0; i<size; i++){
      if(i==(size-1)){
	if(theta >= thetas[i]){
	  //Interpolate to get the value of sigma
	  return ((theta-thetas[i])*(sigmasW[0]-sigmasW[i])/(2*pi+thetas[0]-thetas[i]))+sigmasW[i];
	}
      }
      else{
	if(theta >= thetas[i] && theta < thetas[i+1]){
	  //Interpolate to get the value of sigma	  
	  return ((theta-thetas[i])*(sigmasW[i+1]-sigmasW[i])/(thetas[i+1]-thetas[i]))+sigmasW[i];
	}
      } 
    }
    return 0;
  };
  amrex::Real DW(amrex::Real theta)
  {
   // put in code so that 0 <= theta < 2*pi
    if(theta>=2*pi)
      {
	int div=theta/(2*pi);
	theta=theta-(div*2*pi);
      }
    else if(theta<0)
      {
	int div=fabs(theta)/(2*pi);
	theta=theta+((div+1)*2*pi);
      }
    
    //Find out the value of theta inside our vector thetas
    
    for(int i=0; i<size; i++){
      if(i==(size-1)){
	if(theta >= thetas[i]){
	  //Interpolate to get the value of sigma
	  return ((theta-thetas[i])*(sigmasDW[0]-sigmasDW[i])/(2*pi+thetas[0]-thetas[i]))+sigmasDW[i];
	}
      }
      else{
	if(theta >= thetas[i] && theta < thetas[i+1]){
	  //Interpolate to get the value of sigma	  
	  return ((theta-thetas[i])*(sigmasDW[i+1]-sigmasDW[i])/(thetas[i+1]-thetas[i]))+sigmasDW[i];
	}
      } 
    }
    return 0;
  };
  amrex::Real DDW(amrex::Real theta)
  {
    // put in code so that 0 <= theta < 2*pi
    if(theta>=2*pi)
      {
	int div=theta/(2*pi);
	theta=theta-(div*2*pi);
      }
    else if(theta<0)
      {
	int div=fabs(theta)/(2*pi);
	theta=theta+((div+1)*2*pi);
      }
    
    //Find out the value of theta inside our vector thetas
    
    for(int i=0; i<size; i++){
      if(i==(size-1)){
	if(theta >= thetas[i]){
	  //Interpolate to get the value of sigma
	  return ((theta-thetas[i])*(sigmasDDW[0]-sigmasDDW[i])/(2*pi+thetas[0]-thetas[i]))+sigmasDDW[i];
	}
      }
      else{
	if(theta >= thetas[i] && theta < thetas[i+1]){
	  //Interpolate to get the value of sigma	  
	  return ((theta-thetas[i])*(sigmasDDW[i+1]-sigmasDDW[i])/(thetas[i+1]-thetas[i]))+sigmasDDW[i];
	}
      } 
    }
    return 0;
  }; 

private:
  std::vector<amrex::Real> thetas, sigmasW, sigmasDW, sigmasDDW;
  int size;
};
}
}
}
#endif
