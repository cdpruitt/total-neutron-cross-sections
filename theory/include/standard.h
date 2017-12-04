#include "disperse.h"

/**
 *\brief imaginary surface potential and dispersive correction
 *
 *this class deals with the energy-dependence of the 
 * magnitude of the imaginary potential, which is parametrized as
 *\f$ W(E) = \frac{A}{2\sqrt{2\pi}\sigma} \frac{\|X\|^{m}}{\|X\|^{m}+B^m} \exp\left(-C \|X\|\right) \f$
 * where \f$ X = E-E_{Fermi} \f$ 
 * It also calcuated the dispersive correction
 */


class standard : public disperse
{
 protected:
  double A;
  double B;
  double C;
  double E0;
  double Ef; //!< Fermi energy in MeV
  double m;

 public:
  standard(){};
  standard(double A, double B, double C0,  double E0, double m, double Ef); 
  void init(double A, double B, double C0, double E0, double m, double Ef); 

  double funct(double E);
  double derFunct(double E);

  double functX(double x);
  double derFunctX(double x);

  double deltaV(double E);
  double derDeltaV(double E);

  double derDeltaVHole(double E);
  double derDeltaVParticle(double E);

};
