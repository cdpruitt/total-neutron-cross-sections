#ifndef disperse_
#define disperse_

#include <cmath>
#include <cstdlib>

using namespace std;


/**
 *!\brief base class for numerical dispersive corrections
 *
 * This base class is responsible for numerically deriving the
 * dispersive correction, the energy derivative of this correction
 * and the factor used used to calculate the occupanacy in the DOM
 */
class disperse
{
 public:
  disperse(){};
  virtual ~disperse(){};
  virtual double functX(double)=0;
  virtual double derFunctX(double)=0;
  double deltaVX(double);
  double derDeltaVX(double);
  double derDeltaVX();
  double derDeltaVHoleX(double);
  double derDeltaVParticleX(double);
  double derDeltaVParticleHoleX();  
  double deltaVXAboveBelow(double);
  double deltaVAbove();
  double deltaVBelow();
  double derDeltaVAbove();
  double derDeltaVBelow();
  double findHalfIntegral();
  double halfIntegral;

 private:
  double sumDerivative; //!< energy derivative of dispersive correction.
  double sumPH; //!<particle hole factor
  static double const pi; //!< 3.14159
  double functXzeroBelow(double x);
  double functXzeroAbove(double x);
  double derFunctXzeroBelow(double x);
  double derFunctXzeroAbove(double x);
  double sumDerivativeAbove;
  double sumDerivativeBelow;
  double sumDispersiveAbove;
  double sumDispersiveBelow;


};

#endif
