#ifndef _compound
#define _compound

#include <string>
#include <iostream>
#include "legendre.h"
#include "channel.h"

using namespace std;
/**
 *\brief compound-elastic contribution to elastic scattering
 *
 *Deals with the decay of the compound nucleus and the determiniation
 *of the compound elastic contribution
 */

class compound
{
 public:
  virtual double  TransCoef(int,double)=0; //!< transmission coeff
  compound();
  compound(string*);
  ~compound();
  double statistical(double,double);
  double DifferentialXsectionCE(double);
  void corrections();
  double decayElastic(level);

  int okay;
  double sigmaAbsorption; //total absorption cross section 
  double sigmaCompoundElastic;  //total
  int Ntll;
  double *xsec; //array of cross section for lwave of the exit channel
  double *xsec0; //array of cross section for lwave of the exit channel
  double *xsec1; //array of cross section for lwave of the exit channel

 private:
  static int const  dim;
  channel *elastic;  // channel of elastic scatter
  channel *other;

  // if proton scattering then elastic==proton, other==neutron
  // if neutron scattering then elastic==neutron, other==proton


  double T_tot; // sum of all transmission coefficients



  double *norm0;
  double *norm1;
  double *T; //array of transmission coefficients
  double *mult; //multiplicity of channel - for continuum level density
  double *nu; //array of degrees of freedom
  int Nelastic; // stores the location in the T array for aligned elastic decay
  int Lelastic; //stores the exit channel L wave for aligned decay
  int Aligned; 
  int Nch; // number of open channels

  double decay(level,channel*);




  string elasticState;
  string elasticTrans;
  string otherState;
  string otherTrans;
  double Qother;


};

#endif
