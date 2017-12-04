#include "../include/imaginaryForm.h"

#include <iostream>
using namespace std;
//this class deals with the imaginary potential in the dispersive
//optical model. It does not deal with the position dependence, but only the 
//energy dependence. Also give the correction to the real potential from
//the imaginary due to the dispersive relation.
//It also gives quantities needed to calculate the occupancy of levels
//in the dispersive optical model

double const imaginaryForm::pi = acos(-1.);
double const imaginaryForm::EulerGamma= 0.577215;

/**
 * default constructor
 */
imaginaryForm::imaginaryForm()
{
}
//****************************************************
/**
 * Constructor: energy dependence of imgainary potential has form
 * \f$W(E) = A0\, \Theta(X) \frac{X^{m0}}{X^{m0}+B0^{m0}} exp(-C0\, X) \f$
 * where \f$X=|E-E_{fermi}| - E_{p}\f$ and \f$\Theta(X)\f$ is Heaviside's
 * step function
\param A0 defines the magnitude
\param B0 defines the rise away from the Fermi energy
\param C0 defines the exponential decay
\param Ep gap width, \f$W=0\f$ for \f$|E-E_{fermi}| < E_{p} \f$
\param Efermi0 is the Fermi Energy in MeV
\param mo is the power, either 2 or 4 are allowed
\param asymmetric0 indicates the asymmetric correction is applied
\param alpha0 gives the magnitude of the correction for \f$E>Efermi+Ea0\f$
\param Ea0 is the gap for the correction
*/


imaginaryForm::imaginaryForm(double A0, double B0, double C0,
                             double Ep0, double Efermi0, int m0,
    int asymmetric0,double alpha0, double Ea0)
{
  init(A0,B0,C0,Ep0,Efermi0,m0,asymmetric0,alpha0,Ea0);
}
//************************************************************************
/**
 * initialize the  energy dependence of imgainary potential which has form
 * \f$W(E) = A0\, \Theta(X) \frac{X^{m0}}{X^{m0}+B0^{m0}} exp(-C0\, X) \f$
 * where \f$X=|E-E_{fermi}| - E_{p}\f$ and \f$\Theta(X)\f$ is Heaviside's
 * step function
\param A0 defines the magnitude
\param B0 defines the rise away from the Fermi energy
\param C0 defines the exponential decay
\param Ep gap width, \f$W=0\f$ for \f$|E-E_{fermi}| < E_{p} \f$
\param Efermi0 is the Fermi Energy in MeV
\param mo is the power, either 2 or 4 are allowed
\param asymmetric0 indicates the asymmetric correction is applied
\param alpha0 gives the magnitude of the correction for \f$E>Efermi+Ea0\f$
\param Ea0 is the gap for the correction
*/
void imaginaryForm::init(double A0, double B0, double C0,
                             double Ep0, double Efermi0, int m0, 
                       int asymmetric0, double alpha0, double Ea0)
{
  A = A0;
  B = B0;
  C = C0;
  Ep = Ep0;
  Efermi = Efermi0;
  m = m0;
  if (C != 0.)
    {
      double CB = C*B;
      if (m==2)
	{
          Sin = sin(CB);
          Cos = cos(CB);
          ExpIntegral.SiCi(CB);
          Si = ExpIntegral.Si - pi/2.;
          Ci = ExpIntegral.Ci;
	}
      else if (m == 4)
	{

         if (CB > .5)
           {

             double const coef[4][8] = {
 {.999993,-1.10947,.768440,-.445884,.223283,-.0847743,.0200403,-.00212104},
 {1.00002,-.00335547,-.738522,.796650,-.527578,.230252,-.0589380,.00657945},
 {1.99997,.00642687,-.0901795,-.498564,.590940,-.321721,.0909345,-.0106643},
 {6.00002,-.00480085,.0723877,-.536746,.243927,-.0194388,-.0129787,.00281113}};

             if (CB > 2. || CB < 0.) 
	       cout << "bad CB value " << CB << " C= " << C << " B= " << B 
                    <<  endl;
             for (int i=0;i<4;i++)
	      {
	        Integral[i] =0.;
	        for (int j=0;j<8;j++) Integral[i] += coef[i][j]*pow(CB,j);
	      }
	   }
	 else  // expansion about CB = 0
	   {
	     Integral[0] = 1. - pi/2./sqrt(2.)*CB + pi/4.*pow(CB,2) - 
	       pi/4./sqrt(2.)*pow(CB,3) 
               + pow(CB,4)/36.*(11.-6.*EulerGamma-6.*log(CB)) + 
	       pi/48./sqrt(2.)*pow(CB,5);

	     Integral[1] = 1. - pi/4.*pow(CB,2) + pi/2./sqrt(2.)*pow(CB,3) +
	       pow(CB,4)/4.*(-3.+2.*EulerGamma+2.*log(CB)) -
	       pi/12./sqrt(2.)*pow(CB,5);

             Integral[2] = 2. - pi/2./sqrt(2.)*pow(CB,3) +
	       (1.-EulerGamma-log(CB))*pow(CB,4) + pi/4./sqrt(2.)*pow(CB,5);

	     Integral[3] = 6. + (EulerGamma + log(CB))*pow(CB,4) - 
	       pi/2./sqrt(2.)*pow(CB,5);
	   } 
	}

    }
  if (m == 2)EpB2 = pow(Ep,2) + pow(B,2);
  else if (m == 4) EpB4 = pow(Ep,4) + pow(B,4);
  else cout << " unknown m value " << endl;

  asymmetric = asymmetric0;
  if (asymmetric)
    {
      Ea = Ea0;
      alpha = alpha0;
      El = Efermi + Ea;
      Em = Efermi - Ea;
      //      asymmetric = 0;
      //Wm = ImaginaryPotential(Em);
      //asymmetric = 1;
      if (C == 0.) Wm = -A;
      else 
	{
         Wm = -A*exp(-C*fabs(Ea-Ep));
         ExpIntegral.SiCi(C*Ea);
         Cfact0 = cos(Ea*C)*(ExpIntegral.Si - pi/2.)- sin(Ea*C)*ExpIntegral.Ci;
         Cfact1 = cos(Ea*C)*ExpIntegral.Ci + sin(Ea*C)*(ExpIntegral.Si-pi/2.);
         Cfact2 = exp(C*Ea)*ExpIntegralEi(-C*Ea);
	}
      Asy.init(alpha,C,Ea,Efermi);
    }
}
//****************************************************************************
  /**
   * returns the imaginary potential at the specified energy
   \param E is the energy in MeV
  */
double imaginaryForm::ImaginaryPotential(double E)
{
  if (abs(E-Efermi) < Ep) return 0.;
  double Estar = abs(E-Efermi) - Ep;
  double pot= A*pow(Estar,m)/(pow(Estar,m)+pow(B,m))*exp(-C*Estar);
  // cout << "A= " << A << " B = " << B << " c= " << C << "Estar= " << Estar
  //     << " m = " << m << " Ep = " << Ep << endl;

  if (asymmetric) pot += AsymmetricImaginaryPotential(E);

  return pot;
}
//****************************************************************************
  /**
   * returns the derivative of the imaginary potential \f$ dW/dE \f$
    \param E is the energy in MeV
  */

double imaginaryForm::DerImaginaryPotential(double E)
{
  if (abs(E-Efermi) < Ep) return 0.;
  double Estar = abs(E-Efermi) - Ep;
 
  double fact =  A*exp(-C*Estar)*pow(Estar,m-1)*
  ((float)m*pow(B,m)-C*(pow(B,m)+pow(Estar,m))*Estar)
  /pow(pow(B,m)+pow(Estar,m),2);
  if (E < Efermi) fact = -fact;

  return fact;
}

//***********************************************************************
  /**
   *returns the dispersive correction for m0=2 and C0=0
   * used the analytic expression of VanderKam, J. Phys. G. 26 (2000) 1787
   \param E is the energy in MeV
  */
double imaginaryForm::DeltaZero2(double E)
{
  double Ex = E - Efermi;
  if (Ex == 0.) return 0.;
  else if (Ep ==0) return A*B*Ex/(pow(Ex,2)+pow(B,2));
  else if (abs(Ex) == Ep)return A/pi/(pow(Ex+Ep,2)+pow(B,2))
     /(pow(Ex-Ep,2)+pow(B,2))*
    (pi*B*Ex*(Ex*Ex-Ep*Ep+B*B) + 2.*pow(B*Ep,2)*log(4.*pow(Ep/B,2)));
  else return A/pi/(pow(Ex+Ep,2)+pow(B,2))/(pow(Ex-Ep,2)+pow(B,2))*
    (pi*B*Ex*(Ex*Ex-Ep*Ep+B*B) +(pow(Ex*Ex-Ep*Ep,2)+B*B*(Ex*Ex+Ep*Ep))
     *log(abs((Ex+Ep)/(Ex-Ep))) + 2.*B*B*Ex*Ep*log(abs((Ex*Ex-Ep*Ep)/B/B)));
}
//***********************************************************************
  /**
   *returns the dispersive correction for m0=4 and C0=0
   * used the analytic expression of VanderKam, J. Phys. G. 26 (2000) 1787
   \param E is the energy in MeV
  */
double imaginaryForm::DeltaZero4(double E)
{
  double Ex = E - Efermi;
  if (Ex == 0.) return 0.;

  double fact =pi*B*B*Ex*Ep*(pow(Ex*Ex-Ep*Ep,2)-pow(B,4))    +
       pi*B*Ex/sqrt(2.)*(pow(Ex*Ex-Ep*Ep,3)+B*B*(pow(B,4)+3.*pow(B*Ep,2)
       + pow(B*Ex,2)+pow(Ex,4)+2.*pow(Ex*Ep,2)-3.*pow(Ep,4))) ;

  if (Ex == Ep) fact += 8.*pow(B*Ep,4)*log(4.*pow(Ep/B,2));
  else if (Ex == -Ep) fact -= 8.*pow(B*Ep,4)*log(4.*pow(Ep/B,2));
  else if (Ep!=0.)
     fact += 4.*pow(B,4)*Ex*Ep*(Ex*Ex+Ep*Ep)*log(abs((Ex*Ex-Ep*Ep)/B/B)) +
       (pow(Ex*Ex-Ep*Ep,4)+pow(B,4)*(pow(Ep,4)+6.*pow(Ep*Ex,2)+pow(Ex,4)))*
    log(abs((Ex+Ep)/(Ex-Ep)));

  return A/pi/(pow(Ex+Ep,4)+pow(B,4))/(pow(Ex-Ep,4)+pow(B,4))*fact;

  /*
  return A/pi/(pow(Ex+Ep,4)+pow(B,4))/(pow(Ex-Ep,4)+pow(B,4))*
    (pi*B*B*Ex*Ep*(pow(Ex*Ex-Ep*Ep,2)-pow(B,4))    +
       pi*B*Ex/sqrt(2.)*(pow(Ex*Ex-Ep*Ep,3)+B*B*(pow(B,4)+3.*pow(B*Ep,2)
       + pow(B*Ex,2)+pow(Ex,4)+2.*pow(Ex*Ep,2)-3.*pow(Ep,4)))   +
       4.*pow(B,4)*Ex*Ep*(Ex*Ex+Ep*Ep)*log(abs((Ex*Ex-Ep*Ep)/B/B)) +
       (pow(Ex*Ex-Ep*Ep,4)+pow(B,4)*(pow(Ep,4)+6.*pow(Ep*Ex,2)+pow(Ex,4)))*
       log(abs((Ex+Ep)/(Ex-Ep))));
  */
}
//*********************************************************************
  /**
   *returns the dispersive correction for m0=2 and C0!=0
   * used the analytic expression of VanderKam, J. Phys. G. 26 (2000) 1787
   \param E is the energy in MeV
  */

double imaginaryForm::Delta2(double E)
{
 
  double Ex = E - Efermi;
  if (Ex == 0.) return 0.;
  double Eplus = Ex + Ep;
  double Eminus = Ex - Ep;

  double Cplus = Eplus*C;
  double Cminus = Eminus*C;

  //exponential integrals
  double EiPlus = ExpIntegralEi(-Cplus);
  double EiMinus = ExpIntegralEi(Cminus);

  double fact;
  double denom = (Eminus*Eminus+B*B)*(Eplus*Eplus+B*B);
  fact = -(Cos*Ci+Sin*Si)*4.*B*B*Ex*Ep/denom +
    (Sin*Ci-Cos*Si)*2.*Ex*B*(Ex*Ex+B*B-Ep*Ep)/denom;

  if (abs(Ex) == Ep)
    {
      double stuff = 4.*pow(Ep,2)*exp(2*C*Ep)*ExpIntegralEi(-2.*C*Ep)/
	(4.*pow(Ep,2)+pow(B,2));
      if (Ex == Ep) fact += stuff;
      else fact -= stuff;
    }
  else fact +=  (exp(Cplus)*EiPlus-exp(-Cminus)*EiMinus)*
    (1.-B*B*(Ex*Ex+B*B+Ep*Ep)/denom) +
    (exp(Cplus)*EiPlus+exp(-Cminus)*EiMinus)*2.*B*B*Ep*Ex/denom;
  return fact*A/pi; 
}
//*********************************************************************
  /**
   *returns the dispersive correction for m0=2 and C0!=0
   * used the analytic expression of VanderKam, J. Phys. G. 26 (2000) 1787
   \param E is the energy in MeV
  */
double imaginaryForm::Delta4(double E)
{
  double Ex = E - Efermi;
  if (Ex == 0.) return 0.;
  double Eplus = Ex + Ep;
  double Eminus = Ex - Ep;
  double k = C*B;
  double F1 = pow(Eminus,4)/(pow(Eminus,4)+pow(B,4));
  double F2 = pow(Eplus,4)/(pow(Eplus,4)+pow(B,4));
  double S1 = -(F1-F2)/2./Ex;
  double S2 = 1. - (F1+F2)/2. - Ep*S1;
  double S3 = -2.*Ep + Ep*(F1+F2) - (Ep*Ep+Ex*Ex)/2./Ex*(F1-F2);
  double S4 = 3*Ep*Ep + Ex*Ex - (3.*Ep*Ep+Ex*Ex)/2.*(F1+F2) +
    Ep*(3.*Ex*Ex+Ep*Ep)*(F1-F2)/2./Ex;

  double R1 = (pow(k,4)/12.-2.)*(EulerGamma+log(k)) + pi*(
     -pow(k,2)/4. - 25.*pow(k,4)/144.+pow(k,6)/1440.+ 
     k/sqrt(2.)*(1.+pow(k,2)/6.-pow(k,4)/120.));

  double R2 = (2.*k-pow(k,5)/60.)*(EulerGamma + log(k)) + pi/sqrt(2.)*
    (1.-k*k/2.-pow(k,4)/24.+pow(k,6)/720.)
    -2.*k + pi*pow(k,3)/12. + 137.*pow(k,5)/3600.;

  double R3 = (pow(k,6)/360.-k*k)*(EulerGamma + log(k)) + pi/2. + 3.*k*k/2. 
    -pi*pow(k,4)/48. - 49.*pow(k,6)/7200. - pi*k/sqrt(2.) 
    + pi*pow(k,3)/6./sqrt(2.) + pi*pow(k,5)/120./sqrt(2.);

  double R4 = pow(k,3)/3.*(EulerGamma + log(k)) + pi/sqrt(2.)*(1.+k*k/2.
     -pow(k,4)/24.-pow(k,6)/720.) - pi*k/2.-11.*pow(k,3)/18.+pi*pow(k,5)/240.;


  //exponential integrals
  double EiMinus = 0.;
  double Cminus = Eminus*C;
  if (Eminus != 0.)
    {
     EiMinus = EulerGamma + log(abs(Cminus));
     for (int i=1;i<=19;i++)EiMinus += pow(Cminus,i)/EiTerm[i];
    }
  double EiPlus = 0.;
  double Cplus = Eplus*C;
  if (Eplus != 0.)
    {
     EiPlus = EulerGamma + log(abs(-Cplus));
     for (int i=1;i<=19;i++)EiPlus += pow(-Cplus,i)/EiTerm[i];
    }

  return A*Ex/pi*(R1*S1 + R2*S2/B + R3*S3/B/B + R4*S4/B/B/B) 
    + A/pi*(exp(Cplus)*EiPlus*F2-exp(-Cminus)*EiMinus*F1);
}
//*********************************************************************
  /**
   *returns the dispersive correction to the real potential
   * \f$ \Delta V(E) = \frac{1}{\pi}P\int_{-\infty}^{\infty} W(E') \left( \frac{1}{E`-E}-\frac{1}{E'-E_{Fermi}}\right) dE'\f$
   *see Eq. B.1a in Mahaux+Sator Nucl. Phys. A503(1989)page 525
   \param E is the energy in MeV
   */
double imaginaryForm::DeltaRealPotential(double E)
{

  double out=0.;
  if (C == 0.)
    {
     if (m == 2) out = DeltaZero2(E);
     else if (m == 4) out = DeltaZero4(E);
    }
  else
    {
     if (m==2) out = Delta2(E);
     else if (m==4) out = Delta4Hole(E) + Delta4Particle(E);
    }

  if (asymmetric) out += DeltaAsymmetricHole(E) + DeltaAsymmetricParticle(E);
   return out;

}
//*********************************************************************
  /**
   * returns the value of the integral Eq. B.1b 
   *in Mahaux+Sator Nucl. Phys. A503 (1989) page 525
   * \f$ \frac{1}{\pi}P\int_{E_{Fermi}}^{\infty} W(E') \left( \frac{1}{E`-E}-\frac{1}{E'-E_{Fermi}}\right) dE'\f$
\param E is the energy in MeV
  */
double imaginaryForm::DeltaRealHole(double E)
{
  double out=0.;
  if (C == 0.)
    {
     if (m == 2) out =  DeltaZero2Hole(E);
     else if (m == 4) out =  DeltaZero4Hole(E);
    }
  else
    {
     if (m==2) out =  Delta2Hole(E);
     else if (m==4) out =  Delta4Hole(E);
    }

  if (asymmetric) out += DeltaAsymmetricHole(E);
  return out;
}
//*********************************************************************
  /**
   * returns the value of the integral Eq. B.1c 
   *in Mahaux+Sator Nucl. Phys. A503 (1989) page 525
   * \f$  \frac{1}{\pi}P\int_{-\infty}^{E_{Fermi}} W(E') \left( \frac{1}{E`-E}-\frac{1}{E'-E_{Fermi}}\right) dE'\f$
\param E is the energy in MeV
  */
double imaginaryForm::DeltaRealParticle(double E)
{
  double out=0.;
  if (C == 0.)
    {
     if (m == 2) out =  DeltaZero2Particle(E);
     else if (m == 4) out = DeltaZero4Particle(E);
    }
  else
    {
     if (m==2) out = Delta2Particle(E);
     else if (m==4) out = Delta4Particle(E);
    }
  if (asymmetric) out += DeltaAsymmetricParticle(E);
  return out;
}
//*********************************************************
  /**
   * returns the energy derivative of the function DeltaRealHole 
   * needed to get the occupation probability for hole states.
    \param E is the energy in MeV
   */
double imaginaryForm::DerDeltaHole(double E)
{

  double out=0.;
  if (C == 0.)
    {
     if (m == 2) out = DerDeltaZero2Hole(E);
     else if (m == 4) out = DerDeltaZero4Hole(E);
    }
  else 
    {
     if (m==2) out = DerDelta2Hole(E);
     else if (m==4) out = DerDelta4Hole(E);
    }

  if (asymmetric) out += DerDeltaAsymmetricHole(E);
  return out;
}
//***************************************************
  /**
   * returns the energy derivative of the function DeltaRealParticle 
   * needed to get the occupation probability for particle states.
    \param E is the energy in MeV
   */
double imaginaryForm::DerDeltaParticle(double E)
{
  double out=0.;
  if (C == 0.)
    {
     if (m == 2) out = DerDeltaZero2Particle(E);
     else if (m == 4) out = DerDeltaZero4Particle(E);
    }
  else
    {
     if (m==2) out = DerDelta2Particle(E);
     else if (m==4) out = DerDelta4Particle(E);
    }
  if (asymmetric) out += DerDeltaAsymmetricParticle(E);
  return out;
}
//***************************************************
  /**
   * returns the integral
   * \f$ \frac{1}{\pi}P\int_{E_{Fermi}}^{\infty} W(E') \left( \frac{1}{E`-E}-\frac{1}{E'-E_{Fermi}}\right) dE'\f$ 
   * when m0=2 and C0=0
   */
double imaginaryForm::DeltaZero2Hole(double E)
{
  double Ex = E - Efermi;
  if (Ex == 0.) return 0.;
  if (Ep == 0.) 
     return A*Ex/(pow(Ex,2)+pow(B,2))
    *(B/2. - Ex*log(abs(Ex))/pi);

  else if (Ex == Ep) return A*Ep/EpB2
    *(B/2.+Ep/pi*log(Ep/B));
  else 
    {
      double Eminus = Ex-Ep;
  
      return A/pi/EpB2*(-pow(Ep,2)*log(abs(Eminus/Ep))
      -(Ex-2.*Ep)*Ex*pow(B,2)/(pow(Eminus,2)+pow(B,2))*log(abs(Eminus/B))
      +Ex*B*(Eminus*Ep+pow(B,2))*pi/2./(pow(Eminus,2)+pow(B,2)));
    }
}
//********************************************************************
  /**
   * returns the integral
   * \f$ \frac{1}{\pi}P\int_{-\infty}^{E_{Fermi}} W(E') \left( \frac{1}{E`-E}-\frac{1}{E'-E_{Fermi}}\right) dE'\f$ 
   * when m0=2 and C0=0
   */
double imaginaryForm::DeltaZero2Particle(double E)
{
  double Ex = E - Efermi;
  if (Ex == 0.) return 0.;
  if (Ep == 0) return A*Ex/(pow(Ex,2)+pow(B,2))
    *(B/2. + Ex*log(abs(Ex))/pi);
  else if (Ex == -Ep) return -A*Ep/EpB2*(B/2.+Ep/pi*log(Ep/B));
  else 
    {
      double Eplus = Ex + Ep;
      return A/pi/EpB2*(pow(Ep,2)*log(abs(Eplus/Ep))
      +(Ex+2.*Ep)*Ex*pow(B,2)/(pow(Eplus,2)+pow(B,2))*log(abs(Eplus/B))
      -Ex*B*(Eplus*Ep-pow(B,2))*pi/2./(pow(Eplus,2)+pow(B,2)));
    }
}
//*********************************************************************
  /**
   * returns the energy derivatibe of the function DeltaZero2Hole()
   \param E is the energy in MeV
   */
double imaginaryForm::DerDeltaZero2Hole(double E)
{
  double Ex = E - Efermi;
  double Eminus = Ex - Ep;
  double EB = pow(Eminus,2)+pow(B,2);

  if (Ex == Ep) return A/2./B;
  else if (Ep == 0.) return A*(pow(B,2)-pow(Ex,2))/pow(EB,2)
      *(B/2.-Ex/pi*log(abs(Ex))) -
      A*Ex/pi/EB*(1.+log(abs(Ex)));

  return A/2./pi/pow(EB,2)*(-2.*Eminus*EB + B*(pow(B,2)-pow(Eminus,2))*pi 
			   -4.*pow(B,2)*Eminus*log(abs(Eminus/B)));
}
//***********************************************************************
  /**
   *returns the energy derivative of the function DeltaZero2Particle()
   \parameter E is the energy in MeV
  */
double imaginaryForm::DerDeltaZero2Particle(double E)
{
  double Ex = E - Efermi;
  double Eplus = Ex + Ep;
  double EB = pow(Eplus,2)+pow(B,2);

  if (Ex == -Ep) return A/2./B;

  else if (Ep == 0.) return A*(pow(B,2)-pow(Ex,2))/pow(EB,2)*
    (B/2.+Ex/pi*log(abs(Ex))) + A*Ex/pi/EB*(1.+log(abs(Ex)));

  return A/2./pi/pow(EB,2)*(2.*Eplus*EB + B*(pow(B,2)-pow(Eplus,2))*pi 
			   +4.*pow(B,2)*Eplus*log(abs(Eplus/B)));
}
//**********************************************************************
  /**
   *returns the  Sine Integral function \f$ Si(x)=\int_{0}^{x} \frac{sin(t)}{t}dt \f$
   * this is just an expansion for low x
   \param is the arguement
  */
double imaginaryForm::SinIntegral(double x)
{
  //sine integral of CB
  return x - pow(x,3)/18. + pow(x,5)/600. - pow(x,7)/35280. 
    + pow(x,9)/3265920.;
}
//************************************************************************
  /**
   *returns the  Cosine Integral function \f$ Ci(x)=-\int_{x}^{\infty} \frac{cos(t)}{t}dt \f$
   * this is just an expansion for low x
   \param is the arguement
  */
double imaginaryForm::CosIntegral(double x)
{
  return  EulerGamma + log(x) - pow(x,2)/4. + pow(x,4)/96. 
    - pow(x,6)/4320. + pow(x,8)/322560.;
}
//*************************************************************************
/**
  *returns the exponential Integral \f$Ei(x) = \int_{-\infty}^{x} \frac{exp(t)}{t} dt \f$
  * this is an expansion for small x
  */
double imaginaryForm::ExpIntegralEi(double x)
{

  if (x < 0.) return -ExpIntegral.E1(-x);
  else return ExpIntegral.Ei(x);

  /*
  double y = EulerGamma + log(abs(x));

  for (int i=1;i<=19;i++)
      y += pow(x,i)/EiTerm[i];

  return y;
  */
}
//***********************************************************************
  /**
   * returns the integral
   * \f$ \frac{1}{\pi}P\int_{E_{Fermi}}^{\infty} W(E') \left( \frac{1}{E`-E}-\frac{1}{E'-E_{Fermi}}\right) dE'\f$ 
   * when m0=2 and C0!=0
   */
double imaginaryForm::Delta2Hole(double E)
{
  double Ex = E - Efermi;
  double Eminus = Ex - Ep;

  double SE2 = pow(Eminus,2)+pow(B,2);

  double fact= -A*B*Ex*(Eminus*Ep+pow(B,2))/EpB2/SE2/pi*(Cos*Si-Sin*Ci) +
    A*Ex*pow(B,2)*(Ex - 2.*Ep)/EpB2/SE2/pi*(Cos*Ci+Sin*Si);

  if (Eminus != 0.) 
    fact -= A*pow(Eminus,2)/pi/SE2*exp(-C*Eminus)*ExpIntegralEi(C*Eminus);
  if (Ep != 0.)
    fact +=  A*pow(Ep,2)/pi/EpB2*exp(C * Ep)*ExpIntegralEi(-C*Ep);
  return fact;
}
//*************************************************************************
  /**
   * returns the integral
   * \f$ \frac{1}{\pi}P\int_{-\infty}^{E_{Fermi}} W(E') \left( \frac{1}{E`-E}-\frac{1}{E'-E_{Fermi}}\right) dE'\f$ 
   * when m0=2 and C0!=0
   */
double imaginaryForm::Delta2Particle(double E)
{
  double Ex = E - Efermi;
  double Eplus = Ex + Ep;

  double SE2 = pow(Eplus,2)+pow(B,2);

  double fact= A*B*Ex*(Eplus*Ep-pow(B,2))/EpB2/SE2/pi*(Cos*Si-Sin*Ci) -
    A*Ex*pow(B,2)*(Ex + 2.*Ep)/EpB2/SE2/pi*(Cos*Ci+Sin*Si);
  if (Eplus != 0) fact +=
    A*pow(Eplus,2)/pi/SE2*exp(C*Eplus)*ExpIntegralEi(-C*Eplus);
  if (Ep != 0.) fact -=
      A*pow(Ep,2)/pi/EpB2*exp(C * Ep)*ExpIntegralEi(-C*Ep);
  return fact;
}
//****************************************************************
  /**
   * returns the energy derivative of Delta2Hole()
   \param E is the energy in MeV
  */

double imaginaryForm::DerDelta2Hole(double E)
{
  double Ex = E - Efermi;
  double Eminus = Ex - Ep;
  double SE1 = pow(B,2) - pow(Eminus,2);
  double SE2 = pow(B,2) + pow(Eminus,2);
  double cEminus = C* Eminus;

  if (Ex == Ep) return A/pi/B*(Sin*Ci - Cos*Si);


  return A/pi/pow(SE2,2)*(
   Eminus*(-SE2 - exp(-cEminus)*ExpIntegralEi(cEminus)
          *(-C*pow(Eminus,3)+pow(B,2)*(2.-cEminus))) +
	  B*Sin*(SE1*Ci+2.*B*Eminus*Si) +B*Cos*(2.*B*Eminus*Ci-SE1*Si)  );
  
}
//***************************************************************
  /**
   * returns the energy derivative of the function Delta2Particle()
\param E is the energy in MeV
  */
double imaginaryForm::DerDelta2Particle(double E)
{
  double Ex = E - Efermi;
  double Eplus = Ex + Ep;
  double SE1 = pow(B,2) - pow(Eplus,2);
  double SE2 = pow(B,2) + pow(Eplus,2);
  double cEplus = C* Eplus;

  if (Ex == -Ep) return A/pi/B*(Sin*Ci - Cos*Si);

  else return A/pi/pow(SE2,2)*(
   Eplus*(SE2 + exp(cEplus)*ExpIntegralEi(-cEplus)
          *(C*pow(Eplus,3)+pow(B,2)*(2.+cEplus))) +
	  B*Sin*(SE1*Ci-2.*B*Eplus*Si) -B*Cos*(2.*B*Eplus*Ci+SE1*Si)  );

}
//********************************************************************
  /**
   * returns the integral
   * \f$ \frac{1}{\pi}P\int_{E_{Fermi}}^{\infty} W(E') \left( \frac{1}{E`-E}-\frac{1}{E'-E_{Fermi}}\right) dE'\f$ 
   * when m0=4 and C0=0
   */

double imaginaryForm::DeltaZero4Hole(double E)
{
  double Ex = E - Efermi;
  double Eminus = Ex - Ep;
  double SE1 = pow(Ep,4) + pow(B,4);
  double SE2 = pow(Eminus,4) + pow(B,4);
  
  double fact = (pow(B,7) + B*pow(Ep*Eminus,3))*Ex*sqrt(2.)/4. + 
    (pow(Eminus,2)-pow(Ep,2))*(pow(B,6)-pow(B*Ep*Eminus,2))/4. +
    (pow(Eminus,3)+pow(Ep,3))*sqrt(2.)*(pow(B,5)+pow(B,3)*Ep*Eminus)/4.;
  if (Eminus != 0. && Ep!= 0.) fact -= pow(Ep*Eminus,4)/pi*log(abs(Eminus/Ep));
  if (Eminus !=0) fact -= pow(B*Eminus,4)/pi*log(abs(Eminus/B)) ;
  if (Ep != 0.) fact += pow(Ep*B,4)/pi*log(Ep/B);
  return A*fact/SE1/SE2; 
}
//***********************************************************************
  /**
   * returns the integral
   * \f$ \frac{1}{\pi}P\int_{-\infty}^{E_{Fermi}} W(E') \left( \frac{1}{E`-E}-\frac{1}{E'-E_{Fermi}}\right) dE'\f$ 
   * when m0=4 and C0=0
   */
double imaginaryForm::DeltaZero4Particle(double E)
{
  double Ex = E - Efermi;
  double Eplus = Ex + Ep;
  double SE1 = pow(Ep,4) + pow(B,4);
  double SE2 = pow(Eplus,4) + pow(B,4);

  double fact = sqrt(2.)*Ex/4.*(B*pow(Ep*Eplus,3)-pow(B,7)) +
     (pow(Eplus,2)-pow(Ep,2))/4.*(pow(B,6)-pow(B*Ep*Eplus,2)) +
      (pow(Eplus,3)-pow(Ep,3))*sqrt(2.)/4.*(pow(B,3)*Ep*Eplus-pow(B,5));

  if (Eplus != 0.) fact -= pow(B*Eplus,4)/pi*log(abs(Eplus/B));
  if (Eplus != 0. && Ep != 0.) fact -= pow(Ep*Eplus,4)/pi*log(abs(Eplus/Ep));
  if (Ep != 0.) fact += pow(Ep*B,4)/pi*log(Ep/B);

  return -A*fact/SE1/SE2;
}
//**********************************************************************
  /**
   * returns the integral
   * \f$ \frac{1}{\pi}P\int_{E_{Fermi}}^{\infty} W(E') \left( \frac{1}{E`-E}-\frac{1}{E'-E_{Fermi}}\right) dE'\f$ 
   * when m0=4 and C0!=0
   */
double imaginaryForm::Delta4Hole(double E)
{
 double Ex = E- Efermi;
  double Eminus = Ex - Ep;
  double SE2 = pow(Eminus,4)+pow(B,4);

  double fact[4];

  fact[0] = (-pow(Ep*Eminus,3) 
              - pow(B,4)*(pow(Ep,2)-Ep*Eminus+pow(Eminus,2)))/C;
  fact[1] = (Ex-2.*Ep)*(pow(Ep*Eminus,2)-pow(B,4))/pow(C,2);
  fact[2] = (-pow(B,4) - Ep*Eminus*(pow(Ep,2)-Ep*Eminus+pow(Eminus,2)))
             /pow(C,3);
  fact[3] = (-pow(Ep,3) + pow(Ep,2)*Eminus - Ep*pow(Eminus,2) + pow(Eminus,3))
             /pow(C,4);




  double factMinus = (6.+2.*C*Eminus+pow(C*Eminus,2)+pow(C*Eminus,3))/pow(C,4);
 if (Eminus != 0.) factMinus +=  -pow(Eminus,4)*exp(-C*Eminus)*
		    ExpIntegralEi(C*Eminus);

  double factp =  (6.-2.*C*Ep+pow(C*Ep,2)-pow(C*Ep,3))/pow(C,4);
  if (Ep > 0.) factp += -pow(Ep,4)*exp(C*Ep)*ExpIntegralEi(-C*Ep);
	        

  //  cout << "Hole " << fact[0] << " " << fact[1] << " " << fact[2] << 
  //    " " << fact[3] << " " << factMinus << " " << factp << " " << Eminus <<endl;

  return A*Ex/pi*(Integral[0]*fact[0] +Integral[1]*fact[1]
		+Integral[2]*fact[2]+Integral[3]*fact[3])/EpB4/SE2 + 
                 A/pi/SE2*factMinus - A/pi/EpB4*factp;

}
//***********************************************************************
  /**
   * returns the integral
   * \f$ \frac{1}{\pi}P\int_{-\infty}^{E_{Fermi}} W(E') \left( \frac{1}{E`-E}-\frac{1}{E'-E_{Fermi}}\right) dE'\f$ 
   * when m0=4 and C0!=0
   */
double imaginaryForm::Delta4Particle(double E)
{
  double Ex = E - Efermi;
  double Eplus = Ex + Ep;
  double SE2 = pow(Eplus,4) + pow(B,4);

  double fact[4];
  
  fact[0] = (pow(Ep*Eplus,3) - pow(B,4)*(pow(Ep,2)+Ep*Eplus+pow(Eplus,2)))/C;
  fact[1] = (Ex + 2.*Ep)*(pow(B,4)-pow(Eplus*Ep,2))/pow(C,2);
  fact[2] = (-pow(B,4) + Ep*Eplus*(pow(Ep,2)+Ep*Eplus+pow(Eplus,2)))/pow(C,3);
  fact[3] = (-pow(Ep,3)-pow(Ep,2)*Eplus-Ep*pow(Eplus,2)-pow(Eplus,3))/pow(C,4);


  
  double factPlus = (6.-2.*C*Eplus + pow(C*Eplus,2) - pow(C*Eplus,3))/pow(C,4);

  if (Eplus != 0.) factPlus += -pow(Eplus,4)*exp(C*Eplus)
		     *ExpIntegralEi(-C*Eplus);
   
  double factp = (6.-2.*C*Ep+pow(C*Ep,2) - pow(C*Ep,3))/pow(C,4);
  if (Ep > 0.) factp+= -pow(Ep,4)*exp(C*Ep)*ExpIntegralEi(-C*Ep);
 
  // cout << "Particle " << fact[0] << " " << fact[1] << " " << fact[2] << 
  //   " " << fact[3] << " " << factPlus << " " << factp << " " << Eplus << endl;

  return A*Ex/pi*(Integral[0]*fact[0] +Integral[1]*fact[1]
		  +Integral[2]*fact[2]+Integral[3]*fact[3])/EpB4/SE2 - 
                  A/pi/SE2*factPlus + A/pi/EpB4*factp;

}
//**********************************************************************
  /**
   *returns the energy derivative of DeltaZero4Hole()
   \param E is the energy in MeV
  */
double imaginaryForm::DerDeltaZero4Hole(double E)
{
  double Ex = E - Efermi;
  double Eminus = Ex - Ep;
  double SE2 = pow(Eminus,4) + pow(B,4);
  double SE1 = pow(Ep,4) + pow(B,4);
  if (Ex == Ep) return A/2./sqrt(2.)/B;

  else if (Ep == 0.) return -A*(4.*pow(B,4)*pow(Ex,3) + 4.*pow(Ex,7) - 
      sqrt(2.)*pow(B,7)*pi - 2.*pow(B,6)*Ex*pi - 
      3.*sqrt(2.)*pow(B,5)*pow(Ex,2)*pi + 
      3.*sqrt(2.)*pow(B,3)*pow(Ex,4)*pi + 2.*pow(B,2)*pow(Ex,5)*pi + 
      sqrt(2.)*B*pow(Ex,6)*pi + 16.*pow(B,4)*pow(Ex,3)*log(abs(Ex/B)))/
      (4.*pow(pow(B,4) + pow(Ex,4),2)*pi);
  
  else return 
   A*(B*SE1*(sqrt(2.)*(pow(B,4) + 4.*pow(B,2)*pow(Eminus,2) + pow(Eminus,4)) + 
    2.*B*(pow(B,2) + pow(Eminus,2))*Eminus)*
    (pow(B,2) - pow(Eminus,2))*pi - 4.*pow(Eminus,3)*(SE1*SE2 + 
    4.*pow(B,4)*(pow(B,4)*log(abs(Eminus/B)) + 
    pow(Ep,4)*(log(Ep/B) + log(abs(Eminus/Ep))))))/
   (4.*(pow(B,4) + pow(Ep,4))*pow(SE2,2)*pi);

}
//***********************************************************************
  /**
   * returns the energy derivative of the function DeltaZero4Paqrticle()
\param E is the energy in MeV
  */
double imaginaryForm::DerDeltaZero4Particle(double E)
{
  double Ex = E - Efermi;
  double Eplus = Ex + Ep;
  double SE1 = pow(Ep,4) + pow(B,4);
  double SE2 = pow(Eplus,4) + pow(B,4);

  if (Ex == -Ep) return A/2./sqrt(2.)/B;

  else if (Ep == 0.) return A*(B*(B - Ex)*(B + Ex)*
     (-2.*B*Ex*(pow(B,2)+pow(Ex,2)) + 
     sqrt(2.)*(pow(B,4) + 4.*pow(B,2)*pow(Ex,2) + pow(Ex,4)))*pi + 
     4.*pow(Ex,3)*(SE2 + 4.*pow(B,4)*log(abs(Ex/B))))/
    (4.*pow(pow(B,4) + pow(Ex,4),2)*pi);

  else return A*(SE1*(4.*pow(Eplus,3)*SE2 + B*(B - Eplus)*(B + Eplus)*
       (-2.*B*Eplus*(pow(B,2) + pow(Eplus,2)) + sqrt(2.)*(pow(B,4) + 
        4.*pow(B,2)*pow(Eplus,2) + pow(Eplus,4)))*pi) + 
       16.*pow(B,4)*pow(Eplus,3)*(pow(B,4)*log(abs(Eplus)/B) + 
        pow(Ep,4)*(log(Ep/B) + log(abs(Eplus)/Ep))))/(4.*SE1*pow(SE2,2)*pi);

}
//**********************************************************************
  /**
   *returns the energy derivative of the function Delta4Hole(E)
   \param E is the energy in MeV
  */
double imaginaryForm::DerDelta4Hole(double E)
{
  double Ex = E - Efermi;
  double Eminus = Ex - Ep;

  double dfact[4];

  dfact[0] = (A*(-3.*pow(B,4) + pow(Eminus,4))*pow(Eminus,2))/
             (C*pow(pow(B,4) + pow(Eminus,4),2)*pi);

  dfact[1] = (2.*A*(-pow(B,4) + pow(Eminus,4))*Eminus)/
             (pow(C,2)*pow(pow(B,4) + pow(Eminus,4),2)*pi);

  dfact[2] = -((A*(pow(B,4) - 3.*pow(Eminus,4)))/
	     (pow(C,3)*pow(pow(B,4) + pow(Eminus,4),2)*pi));

  dfact[3] = (4.*A*pow(Eminus,3))/
             (pow(C,4)*pow(pow(B,4) + pow(Eminus,4),2)*pi);

  double dfactMinus = (pow(B,4)*(2.*C - pow(C,4)*pow(Eminus,3) + 
             pow(C,2)*Eminus*(2. + 3.*C*Eminus))- pow(Eminus,3)*
             (24. + 6.*C*Eminus + pow(C,4)*pow(Eminus,4) + 
	      pow(C,2)*pow(Eminus,2)*(2. + C*Eminus)))/pow(C,4);

	     if (Eminus != 0.) dfactMinus  -= 
	     exp(-C*Eminus)*pow(Eminus,3)*(-C*pow(Eminus,5) + 
	     pow(B,4)*(4. - C*Eminus))*ExpIntegralEi(C*Eminus);

             dfactMinus *= A/pi/pow(pow(Eminus,4)+pow(B,4),2);


 return (Integral[0]*dfact[0] +Integral[1]*dfact[1]
	 +Integral[2]*dfact[2]+Integral[3]*dfact[3])+dfactMinus;

}
//**********************************************************************
  /**
   *returns the energy derivative of the function Delta4Particle(E)
   \param E is the energy in MeV
  */
double imaginaryForm::DerDelta4Particle(double E)
{
  double Ex = E - Efermi;
  double Eplus = Ex + Ep;

  double dfact[4];
  dfact[0] = (A*pow(Eplus,2)*(-3.*pow(B,4) + pow(Eplus,4)))/
              (C*pow(pow(B,4) + pow(Eplus,4),2)*pi);

  dfact[1] = (-2.*A*Eplus*(-pow(B,4) + pow(Eplus,4)))/
             (pow(C,2)*pow(pow(B,4) + pow(Eplus,4),2)*pi);

  dfact[2] = (A*(-pow(B,4) + 3*pow(Eplus,4)))/
             (pow(C,3)*pow(pow(B,4) + pow(Eplus,4),2)*pi);

  dfact[3] = (-4.*A*pow(Eplus,3))/
            (pow(C,4)*pow(pow(B,4) + pow(Eplus,4),2)*pi);


  double dfactPlus = -((pow(Eplus,3)*(24. - 6.*C*Eplus + 
               pow(C,4)*pow(Eplus,4) - pow(C,2)*pow(Eplus,2)*
               (-2. + C*(Eplus))) + pow(B,4)*(2.*C + pow(C,4)*pow(Eplus,3) 
               + pow(C,2)*Eplus*(-2. + 3.*C*Eplus)))/pow(C,4));

    if (Eplus != 0.) dfactPlus -= exp(C*Eplus)*pow(Eplus,3)*(C*pow(Eplus,5) + 
	       pow(B,4)*(4. + C*Eplus))*ExpIntegralEi(-(C*Eplus));


  dfactPlus *= A/pi/pow(pow(B,4) + pow(Eplus,4),2);

  return (Integral[0]*dfact[0] +Integral[1]*dfact[1]
	 +Integral[2]*dfact[2]+Integral[3]*dfact[3]) - dfactPlus;

}
//*****************************************
  /**
   * calculates the asymmetric term in the imaginary potential
\para E is the energy in MeV
   */
double imaginaryForm::AsymmetricImaginaryPotential(double E)
{
  if (E > El) 
    {
     if (C == 0.) return alpha*(sqrt(E) + pow(El,3./2.)/2./E - 3./2.*sqrt(El));
     return Asy.W(E);
    }
  if (E < Em)
    { 
      return Wm*pow(E-Em,2)*exp(-C*fabs(E-Em))
                            /(pow(E-Em,2)+pow(Ea,2));
    }
  else return 0.;
}
//******************************************
  /**
   * returns the value of the integral 
   * \f$ \frac{1}{\pi}P\int_{E_{Fermi}}^{\infty} W_{asy}(E') \left( \frac{1}{E`-E}-\frac{1}{E'-E_{Fermi}}\right) dE'\f$
   * where \f$W_{asy}\f$ is the asymmetric correction
\param E is the energy in MeV
  */
double imaginaryForm::DeltaAsymmetricHole(double E)
{
  if(C!=0.) return Asy.dispersive(E);
  double term;
  if (E == El) 
    {
     term = sqrt(El)/2.*(log(16.)+3.*log(abs(El/Ea)));
     term -= pow(El,3./2.)/2./Efermi*log(abs(El/Ea));
     term += 2.*sqrt(abs(Efermi))*(pi/2.-atan(sqrt(abs(El/Efermi))));
     return term*alpha/pi;
    }  

  if (E < 0.) term = -2.*sqrt(abs(E))*(pi/2.-atan(sqrt(abs(El/E))));
  else term = sqrt(E)*log(abs((sqrt(El)+sqrt(E))/(sqrt(El)-sqrt(E))));

  term += 2.*sqrt(abs(Efermi))*(pi/2.-atan(sqrt(El/abs(Efermi))));

  if (E == 0.) term += pow(El,3./2.)/2.*(1./El - log(abs(El/Ea))/Efermi);
  else term += pow(El,3./2.)/2./E/Efermi*(Efermi*log(abs(El/(El-E))) 
					  - E*log(abs(El/Ea)));

  term += 3./2.*sqrt(El)*log(abs((El-E)/Ea));

  return term*alpha/pi;

}
//*******************************************
  /**
   * prints out the values of the paramters defining the imaginary potenetial
   */
void imaginaryForm::Print()
{
  cout << "A= " << A << endl;
  cout << "B= " << B << endl;
  cout << "C= " << C << endl;
  cout << "m= " << m << endl;
  cout << "asymmetric= " << asymmetric << endl;
  if (asymmetric)
    {
     cout << "alpha = " << alpha << endl;
     cout << "Ea = " << Ea << endl;
    }
}
//*********************************************************
  /**
   * returns the energy derivative of the function DeltaAsymmetricHole
   \param E is the energy in MeV
  */
double imaginaryForm::DerDeltaAsymmetricHole(double E)
{
  if (C != 0) return Asy.DerDispersive(E);
  if (E == El) return alpha*(1.+log(4.))/2./pi/sqrt(El);
  else if (abs(E) < 0.01 )  return 3.*alpha/4./sqrt(El)/pi
   - 5./2.*E/pow(El,3./2.)/pi*alpha;

  else if (E > 0.)   return alpha/2./pi/pow(E,2)*
      (pow(E,3./2.)*log((sqrt(El)+sqrt(E))/abs(sqrt(El)-sqrt(E))) + 
       sqrt(El)*(E-El*log(abs(El/(El-E)))));

  else   return alpha/2./pi*(2.*sqrt(abs(El))/(abs(El)-E) + 
   (pi-2.*atan(sqrt(abs(El/E))))/sqrt(abs(E)) + 
   sqrt(El)*((El-3.*E)*E + El*(E-El)*log(El/(El-E)))/pow(E,2)/(El-E));

}
//******************************************
  /**
   * returns the value of the integral 
   * \f$  \frac{1}{\pi}P\int_{-\infty}^{E_{Fermi}} W_{asy}(E') \left( \frac{1}{E`-E}-\frac{1}{E'-E_{Fermi}}\right) dE'\f$ 
   *where \f$ W_{asy} \f$ is the asymmetric correction
\param E is the energy in MeV
  */
double imaginaryForm::DeltaAsymmetricParticle(double E)
{

  double Ex = E - Efermi;
  if (Ex == 0.) return 0.;
  if (C == 0) 
    {
     if (Ea == 0) return Wm*Ex/(pow(Ex,2)+pow(Ea,2))
       *(Ea/2. + Ex*log(abs(Ex))/pi);
     else if (Ex == -Ea) return -Wm/4.;
     else 
       {

         //return 0.;
	 double Eplus = Ex + Ea;
         return Wm/pi/2.*(log(abs(Eplus/Ea))
         +(Ex+2.*Ea)*Ex/(pow(Eplus,2)+pow(Ea,2))*log(abs(Eplus/Ea))
         -pow(Ex,2)*pi/2./(pow(Eplus,2)+pow(Ea,2))); 
       }
    }
  else
    {
      double Eplus = Ex + Ea;
      
      return Wm/pi/(pow(Eplus,2)+pow(Ea,2))*(pow(Ex,2)/2.*Cfact0 - 
      Ex*(Ex+2.*Ea)/2.*Cfact1 +
       pow(Eplus,2)*exp(C*Eplus)*ExpIntegralEi(-C*Eplus)) 
     -Wm/pi*Cfact2/2.;

    }
}
//********************************************
  /**
   * returns the energy derivative of the function DeltaAsymmetricParticle(E)
\param E is the energy in MeV
  */
double imaginaryForm::DerDeltaAsymmetricParticle(double E)
{
  double Ex = E - Efermi;
  double Eplus = Ex + Ea;
  double EB = pow(Eplus,2)+pow(Ea,2);
  if (C == 0)
    {
     if (Ex == -Ea) return Wm/2./Ea;

      return Wm/2./pi/pow(EB,2)*(2.*Eplus*EB + Ea*(pow(Ea,2)-pow(Eplus,2))*pi 
			   +4.*pow(Ea,2)*Eplus*log(abs(Eplus/Ea)));
    }
  else
    {
 return -2.*Eplus*Wm/pi/pow(EB,2)*(pow(Ex,2)/2.*Cfact0 - 
      Ex*(Ex+2.*Ea)/2.*Cfact1 +
       pow(Eplus,2)*exp(C*Eplus)*ExpIntegralEi(-C*Eplus)) 

       +Wm/pi/EB*(Ex*Cfact0 - Eplus*Cfact1 +
       2.*Eplus*exp(C*Eplus)*ExpIntegralEi(-C*Eplus) + 
       C*pow(Eplus,2)*exp(C*Eplus)*ExpIntegralEi(-C*Eplus) +
       Eplus) 
     ;
    }
}
//******************************************************************
  /**
   * returns the function Delta4Hole as an expansion
  \param E is the energy in MeV
   */
double imaginaryForm::Delta4HoleExpand(double E)
{
  double Ex = E- Efermi;
  double Eminus = Ex - Ep;
  double SE2 = pow(Eminus,4)+pow(B,4);

  double term0 =(B*Ex*((pow(B,2) + Ep*Eminus)*
   (-(B*(pow(B,2) - Ep*Eminus)*(2.*Ep - Ex)) + sqrt(2.)*(pow(B,4) + 
   pow(Ep,2)*pow(Eminus,2) + pow(B,2)*pow(-2.*Ep + Ex,2)))*pi - 
   4.*pow(B,3)*(2.*Ep - Ex)*(2.*pow(Ep,2) - 2.*Ep*Ex + pow(Ex,2))*log(B) - 
   4.*pow(B,3)*(2.*Ep - Ex)*(2.*pow(Ep,2) - 2.*Ep*Ex + pow(Ex,2))*
   (EulerGamma + log(C))))/(4.*(pow(B,4) + pow(Ep,4))*
   (pow(B,4) + pow(Eminus,4)));

 if (Ex != Ep) term0 += -pow(Eminus,4)*exp(-C*Eminus)*
		    ExpIntegralEi(C*Eminus) /SE2;

 if (Ep > 0.)  term0 += pow(Ep,4)*exp(C*Ep)*ExpIntegralEi(-C*Ep)/EpB4;

  double term1 = (pow(B,2)*Ex*(sqrt(2.)*B*
    (-pow(B,4) + pow(Ep,2)*pow(Eminus,2))*(-2.*Ep + Ex)*pi + 
    sqrt(2.)*pow(B,3)*(2.*Ep - Ex)*(2.*pow(Ep,2) - 2.*Ep*Ex + pow(Ex,2))*pi + 
    (-pow(Ep,3)*pow(Eminus,3) - pow(B,4)*(3.*pow(Ep,2) - 3.*Ep*Ex + 
    pow(Ex,2)))*pi + 4.*pow(B,2)*(pow(B,4) + Ep*Eminus*
    (3.*pow(Ep,2) - 3.*Ep*Ex + pow(Ex,2)))*(-1.+EulerGamma+log(B) + log(C))))/
    (4.*(pow(B,4) + pow(Ep,4))*(pow(B,4) + pow(Eminus,4)));

 return A/pi*(term0 + C*term1);

}
//******************************************************************
  /**
   * gives the function Delta4Hole as an expansion
   */

double imaginaryForm::Delta4ParticleExpand(double E)
{
  double Ex = E- Efermi;
  double Eplus = Ex + Ep;
  double SE2 = pow(Eplus,4)+pow(B,4);

  double term0 = (B*Ex*((pow(B,2) - Ep*(Ep + Ex))*
        (-(B*(2.*Ep + Ex)*(pow(B,2) + Ep*(Ep + Ex))) + 
          sqrt(2.)*(pow(B,4) + 
             pow(Ep,2)*pow(Ep + Ex,2) + 
             pow(B,2)*pow(2.*Ep + Ex,2)))*pi - 
       4.*pow(B,3)*(2.*Ep + Ex)*
        (2.*pow(Ep,2) + 2.*Ep*Ex + pow(Ex,2))*log(B) - 
       4.*pow(B,3)*(2.*Ep + Ex)*
        (2.*pow(Ep,2) + 2.*Ep*Ex + pow(Ex,2))*
        (EulerGamma + log(C))))/
   (4.*(pow(B,4) + pow(Ep,4))*
    (pow(B,4) + pow(Ep + Ex,4)));

 if (Ex != -Ep) term0 +=  pow(Eplus,4)*exp(C*Eplus)
		     *ExpIntegralEi(-C*Eplus)/SE2;

 if (Ep > 0.)  term0 +=  -pow(Ep,4)*exp(C*Ep)*ExpIntegralEi(-C*Ep)/EpB4;


double term1 = (pow(B,2)*Ex*(sqrt(2.)*pow(B,3)*(2.*Ep + Ex)*
        (2.*pow(Ep,2) + 2.*Ep*Ex + pow(Ex,2))*pi + 
       sqrt(2.)*B*(2.*Ep + Ex)*
        (pow(B,4) - pow(Ep,2)*pow(Ep + Ex,2))*pi + 
       (pow(Ep,3)*pow(Ep + Ex,3) - 
          pow(B,4)*(3.*pow(Ep,2) + 3.*Ep*Ex + 
             pow(Ex,2)))*pi + 
       4.*pow(B,2)*(pow(B,4) - 
          Ep*(Ep + Ex)*
           (3.*pow(Ep,2) + 3.*Ep*Ex + pow(Ex,2)))*
        (-1. + EulerGamma + log(B) + log(C))))/
   (4.*(pow(B,4) + pow(Ep,4))*
    (pow(B,4) + pow(Ep + Ex,4)));


 return A/pi*(term0 + C*term1);
}

