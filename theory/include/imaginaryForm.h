#ifndef imaginaryForm_
#define imaginaryForm_
#include "expInt.h"
#include <cmath>

#include <algorithm>
#include "asy.h"
using namespace std;



double const EiTerm[] = {0.,1.,4.,18.,96.,600.,4320.,3520.,322560.,3265920.,
		     36288000.,439084800.,5748019200.,80951270400.,
		     1220496076800.,19615115520000.,334764638208000.,
		     6046686277632000.,115242726703104000.,
		     2311256907767808000.};

/**
 *\brief energy-dependence of imaginary potentials and disperive corr.
 *
 *This class deals with a parametrized energy dependence of the imaginary
 *potentials and the dispersive corrections to the real potential
 * \f$ \Delta V(E) = \frac{1}{\pi}P\int W(E') \left( \frac{1}{E`-E}-\frac{1}{E'-E_{Fermi}}\right) dE'\f$
 */

class imaginaryForm
{
  public:
  static double const pi;
  static double const EulerGamma;


  imaginaryForm();
  imaginaryForm(double,double,double,double,double,int,int,double,double);
  void init(double,double,double,double,double,int,int,double,double); 
  double ImaginaryPotential(double);

  double DerImaginaryPotential(double);
  double DeltaRealPotential(double);
  double DeltaRealHole(double);
  double DeltaRealParticle(double);

  double DerDeltaHole(double);
  double DerDeltaParticle(double);
  void Print();

  double Delta(double);

  double A;
  double B;
  double C;
  double Ep;
  double Efermi;
  int m;
  double DerDelta; //derivative of DeltaRealPotential w.r.t. energy

  int asymmetric; // logical to indicate if asymmetrical form is used
  double Ea;
  double El; //=Ea+Efermi
  double alpha;
  double Em;  //Efermi - Ea
  double Wm;  


  private:

  double DeltaAsymmetricHole(double);
  double DeltaAsymmetricParticle(double);
  double DerDeltaAsymmetricHole(double);
  double DerDeltaAsymmetricParticle(double);
  double AsymmetricImaginaryPotential(double);
  double DeltaZero4(double);
  double Delta2(double);
  double Delta4(double);
  double SinIntegral(double);
  double CosIntegral(double);
  double ExpIntegralEi(double);

  double DeltaZero2Hole(double);
  double DeltaZero2Particle(double);
  double Delta2Hole(double);
  double Delta2Particle(double);
  double DeltaZero2(double);

  double DerDeltaZero2Hole(double);
  double DerDeltaZero2Particle(double);
  double DerDelta2Hole(double);
  double DerDelta2Particle(double);

  double DeltaZero4Hole(double);
  double DeltaZero4Particle(double);
  double Delta4Hole(double);
  double Delta4Particle(double);

  double DerDeltaZero4Hole(double);
  double DerDeltaZero4Particle(double);
  double DerDelta4Hole(double);
  double DerDelta4Particle(double);

  double Delta4HoleExpand(double);
  double Delta4ParticleExpand(double);
  double DerDelta4HoleExpand(double);
  double DerDelta4ParticleExpand(double);


  double Cos; //cos(B*C)
  double Sin; //sin(B*C)
  double Si; //SinIntegral(B*C) - pi/2.
  double Ci; //CosIntegral(B*C)
  double EpB2; //pow(Ep,2)+pow(B,2)
  double EpB4; //pow(Ep,4)+pow(B,4)
  double Integral[4];

  // used for energy-asymmetric correction if C!=0
  double Cfact0;
  double Cfact1;
  double Cfact2;

  expInt ExpIntegral;
  asy Asy;
  };
#endif 
