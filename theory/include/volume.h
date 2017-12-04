


#ifndef volume_
#define volume_
#include "potPara.h"
#include "imaginaryForm.h"
#include "surVolume.h"

/**
 *\brief imaginay volume potential and dispersive potential
 *
 *This class deals with all aspects of the volume imaginary potential, its
 *dispersive correction and contributions to effetive mass and occupation
 *probabilities
 *
 *there is a main volume component (with constant geometry, i.e. r 
 *dependence), plus a surface correction with acts to approximate a change of
 *radius with energy. At Ecm = Ef the radius is Rzero + deltaR,
 * at Ecm = (+-) infinity the radius is Rzero. 
 */


class volume
{
 public:
  void load(double,double,double,double,double,double,double,double,
   int,int,double,double);
  void load(double strength, double R0, double a0);
  void SetEnergy(double);
  double ImaginaryPot(double);
  double DispersiveCorrection(double);
  double DerivativeDispersive(double);
  double DerivativeParticleHole(double);

  //main volume component
  potPara ImaginaryVolume;  //imaginary potential 
  potPara DispersiveVolume; //dispersive coorection to real potential
  potPara DerDispersiveVolume; //derivative of dispersive needed
                                // for effective mass
  potPara ParticleHoleVolume; // derivative of particle of hole 
                              //contribution to dispersive correction 
                              //needed for occupation probabilities

  //surface correction
  potPara ImaginarySurface;  //imaginary potential 
  potPara DispersiveSurface; //dispersive coorection to real potential
  potPara DerDispersiveSurface; //derivative of dispersive needed
                                // for effective mass
  potPara ParticleHoleSurface; // derivative of particle of hole 
                              //contribution to dispersive correction 
                              //needed for occupation probabilities

  imaginaryForm volume; // main volume contribution with fixed radius
  //imaginaryForm surface; // surface correction to give change of radius
  surVolume surface; // surface correction to give change of radius
   double deltaR; // energy dependence of radius
   double Rzero; // radius at E= Efermi
   double expR; //exponential decrease of radius parameter
   double a;  //diffuseness
   double A;
   double B;
   double Ep;
   double Efermi;  //Fermi energy
   int Asy;  //switch for energy-asymmetry contribution
   int m;
   double alpha;
   double Ea;
   double Ecm; //Energy 

};
#endif
