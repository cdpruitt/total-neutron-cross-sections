#include "../include/scat.h"
#include "../include/waves.h"




//***********************************************************************
double const scat::kconstant =  .048192;
double const scat::e2 = 1.44; // Coulomb constant
double const scat::m0 = 931.5; // nucleon mass
double const scat::pi = acos(-1.);
double const scat::EminRel = 1.;
//*************************************************************************
//default constructor
scat::scat(){}
//**************************************************************************
  /**
   * loads the potential so subsequent operations
\param HartreeFock points to the Hartree-Fock potential
\param Volume0 point to the volume imaginary potential
\param Surface0 points to the  Surface imaginary pot
\param SpinOrbit0 points to the real and imaginary spin-orbit potential
\param Rc0 is the radius in fm for the Coulomb potential
\param energCM0 is the energy in the center-of-mass frame in units of MeV
\param energyLab0 is the energy in the labortory frame
   */
void scat::loadPotential(hartreeFock *HartreeFock0, 
	   volume *Volume0, surfaceTF *SurfaceLow0, surfaceSTD *SurfaceHigh0, 
           spinOrbit * SpinOrbit0, 
            double Rc0,double energyCM0, double energyLab0)
{
  HartreeFock = HartreeFock0;
  Volume = Volume0;
  SurfaceLow = SurfaceLow0;
  SurfaceHigh = SurfaceHigh0;
  SpinOrbit = SpinOrbit0;
  Rc = Rc0;

  energyCM = energyCM0;
  energyLab = energyLab0;
  Kwave2 = kconstant*mu*energyCM;



  //relativistic correction
  if (energyCM > EminRel) Kwave2 = kconstant*pow(A/(1.+ A + energyCM/m0),2)*
		       energyLab*(energyLab/2./m0 + 1.);  

  Kwave = sqrt(abs(Kwave2)); //wave number in fm**-1

  if (energyCM > EminRel) gammaRel = 2.*(energyCM/m0 +1.)/(energyCM/m0 + 2.);
  else gammaRel = 1.;

  if (energyCM > 0.) muhbar = Kwave2/energyCM;
  else muhbar = kconstant*mu;

  gamma = fabs(gammaRel*Z*Zp*e2*Kwave/energyCM/2.); //Sommerfeld parameter
  //gamma = mu*Z*Zp/Kwave/28.820;  // Sommerfeld parameter
  konst = pi/pow(Kwave,2); //const for cross sections in units of Fermi squared
  konst *= 10.; //umits of mb

}
//***************************************************************************
  /**
   * Constructor
\param Zp0 is atomic number of projectile(0=neutron,1=proton)
\param Z0 is atomic number of target
\param A0 is the mass number of the target
\param flag0 indicates integrating of wavefunctions takes place
  */
scat::scat(double Zp0, double Z0, double A0, bool flag0, string *title0)
  :compound(title0)
{
  init(Zp0,Z0,A0,flag0);
}
//************************************************************************
  /**
   * initializes for a new calculation
\param Zp0 is atomic number of projectile(0=neutron,1=proton)
\param Z0 is atomic number of target
\param A0 is the mass number of the target
\param flag0 indicates integrating of wave functions takes place
  */
void  scat::init(double Zp0,double Z0, double A0, bool flag0)

{
  Zp = Zp0;
  flag = flag0;
  A = A0;
  Z = Z0;
  mu = A/(1.+A);
  if (Zp!=0.)proton = 1;
  else proton = 0;
  if (flag == 0) return;
  elastic = 1;
  initIntegration();
}
//*************************************************************************
  /**
   *returns the real total potential (centripetal not included) in MeV
   *Includes the dispersive correction
   *It is scaled by gamma \f$\sqrt{1-\left(\frac{v}{c}\right)^2}\f$
\param r is the radial distance in fm
   */
double scat::RealPotential(double r)
{
  // first nuclear potential

  double one = HartreeFock->RealPotential(r);
  double two = Volume->DispersiveCorrection(r);
  double three = SurfaceLow->DispersiveCorrection(r);
  double four = SurfaceHigh->DispersiveCorrection(r); 
  double fact = one + two + three + four;



  if (proton == 1) // finally the Coulomb potential
    {
      if (r > Rc) fact += e2*Z*Zp/r;
      else fact += e2*Z*Zp/2./Rc*(3.-pow(r/Rc,2));
    }

  if (r > 0.) 
    {
      // next the spin-orbit potential
     fact += SpinOrbit->RealPotential(r,LdotSigma);
    }

  return fact*gammaRel; 
}
//**************************************************************************
  /**
   * return the imaginary potental in MeV
   *It is scaled by gamma \f$\sqrt{1-\left(\frac{v}{c}\right)^2}\f$
   \param r is the radial distance in fm
   */
double scat::ImaginaryPotential(double r)
{
  double fact =  SurfaceLow->ImaginaryPot(r) 
     + SurfaceHigh->ImaginaryPot(r) 
     + Volume->ImaginaryPot(r);
  if (r > 0.) 
    {
      // next the spin-orbit potential
      fact += SpinOrbit->ImaginaryPotential(r,LdotSigma);
    }
  return fact*gammaRel; 
}
//**************************************************************************
valarray<double> scat::diff(double r , valarray<double> w)
{
 //       REAL (kind=8) :: pot_real,pot_imag   
 // this subroutine is used by the mersion subroutine 
 // w[0] = REAL(u(r)) real part of  radial wave function
 //w[1] = REAL(du/dr) derrivative of real part of radial wave function
 //w[2] = IMAG(u(r)) Imaginary part of radial wave function
 //w[3] = IMAG(du/dr) derrivative of imaginary part
 //F(0:3) are the derivaties of w(0:3)




  int n = w.size();
  valarray<double> f(n);
  f[0] = w[1]; //these are equal by definition

  double potReal = RealPotential(r)*muhbar;
  if ( l > 0) potReal += (double)(l*(l+1))/pow(r,2);

  f[1] = -(Kwave2 - potReal)*w[0]; 


  if (n == 4)
    {
      f[2] = w[3]; // equal by definition
      double potImag = ImaginaryPotential(r)*muhbar;
      f[1] -= potImag*w[2]; 
      f[3] = -(Kwave2 - potReal)*w[2] + potImag*w[0];

    }

  return f;
}

//*****************************************************************
  /**
   * For positive energies, this Integrates the wave functions out from 
   * the origin to the matching radius, where the wavefunction is match to
   * Coulomb wavefunctions or spherical Bessel functions.
   * It determines the scattering phase shift
   */

int scat::integrateWave()
{


  // initialize the coupled differential equation solver
  initMerson(.001,.00001,.1,0);

  // find initial and matching wave functions
  //initialise wave function
  double rhoStop = rStop*Kwave; 


  waves outWave(rhoStop,gamma,lMax);
  for (int i=0;i<lMax;i++) Sigma[i] = outWave.Sigma[i];

  valarray<double> WaveFunct(4);

  l = 0;
  lStop = -1;
  for (;;) // loop over orbital angular momentum
    {
      for (int plusMinus = -1;plusMinus<=1;plusMinus+=2)// spin up or spin down
	{
          j = (float)l + (float)plusMinus*0.5;
          if (j < 0.) continue;
          LdotSigma = j*(j+1.) - double(l*(l+1)) - 0.5*1.5 ;


          // potential at start
          double Vstart = RealPotential(rStart);
          double Wstart = ImaginaryPotential(rStart);
          //derivative of potential at start
          double dVstart = (RealPotential(rStart+.01) - Vstart)/0.01;
          double dWstart = (ImaginaryPotential(rStart+.01) - Wstart)/0.01;

	  // initialize wavefunctions
          double fact = pow(rStart,l+3)/2./(double)(2*l+3);
          WaveFunct[0] = pow(rStart,l+1) 
                                - muhbar*(energyCM-Vstart)*fact; // real part
          WaveFunct[2] =  Wstart*muhbar*fact;              //imaginary part

          // derivative of wavefunction
          fact = (double)(l+3)*pow(rStart,l+2)/2./double(2*l+3);
          WaveFunct[1] = (double)(l+1)*pow(rStart,l) 
                           - muhbar*(energyCM-Vstart)*fact; // real
          WaveFunct[3] = muhbar*Wstart*fact;          // imaginary

          fact = muhbar*pow(rStart,l+3)/2./(double)(2*l+3);
          WaveFunct[1] += dVstart*fact;
          WaveFunct[3] += dWstart*fact;


	  //integrate wavefunctions out to matching radius
          solveMerson(&WaveFunct,rStart,rStop);

          if (ok == 0) 
	    {
              cout << "j= " << j << " l= " << l << " Ecm= " <<
		     energyCM << endl;
              return 0;
	    }

          //outWave gives derivates with respect to rho = Kwave*r
	  //but at this point WaveFunctOut have derivative with respect to r
          // so make them with respect to rho=kwave*r

	  WaveFunct[1] /= Kwave;
          WaveFunct[3] /= Kwave;


	  // match wave functions 
	  //real WaveFunct = AA*F + BB*G
	  double  BB = outWave.dF[l]*WaveFunct[0] 
                     - outWave.F[l]*WaveFunct[1];
	  double AA = -outWave.dG[l]*WaveFunct[0] 
                      + outWave.G[l]*WaveFunct[1];

	  // imaginary part => Wavefunct  = CC*F + DD*G
	  double DD = outWave.dF[l]*WaveFunct[2] 
                    - outWave.F[l]*WaveFunct[3];
	  double CC = -outWave.dG[l]*WaveFunct[2] 
                    + outWave.G[l]*WaveFunct[3];


	  double denominator = pow(AA+DD,2) + pow(CC-BB,2);
	  double etaReal = (pow(AA,2)-pow(DD,2) + pow(CC,2) 
                    - pow(BB,2))/denominator;
	  double etaImag = 2.*(AA*BB+CC*DD)/denominator;

          int up = (plusMinus+1)/2;  //=0 form down and =1 for up
	  phaseShift[l][up] = atan2(etaImag,etaReal);
          if (phaseShift[l][up] < 0.) phaseShift[l][up] += 2.*pi;
          phaseShift[l][up] /=2.;
         
	  eta2[l][up] = pow(etaReal,2) + pow(etaImag,2);
          eta[l][up] = complex<double>(etaReal,etaImag);


	}
      //check to see if we have included enough l-waves
      if (1.-eta2[l][1] < 0.0005 && 1.-eta2[l][0] < 0.0005 )
	{
          lStop = l;
	  //if (l > 19) break;
	   break;
	}
      l++;
      if (l > lMax) break;
    }
      if (lStop == -1) cout << "increase lMax" << endl;
      return 1;
}
//********************************************************
  /**
   * returns the elastic center-of-mass differntial cross section.
   * the function integrateWave() must be executed first
   \param theta center-of-mass angle in radians 
  */
double scat::DifferentialXsection(double theta)
{
  complex<double> A(0.,0.);
  complex<double> B(0.,0.);
  complex<double>tempA;
  complex<double>tempB;



  int ll = lMax;
  legendre Poly(ll);


  // start with l=0
  // eta is the S matrix

  A = eta[0][1];
  if (elastic) A -= 1.;

  if (proton) A *=  exp(complex<double>(0.,2.*Sigma[0])); 
  A *= Poly.LegendreP0(0,theta);


  //now the other l's

  for (int i=1;i<=lStop;i++)
    {
      //spin nonflip
      l = i;
      tempA = (double)(l+1)*eta[l][1] + (double)l*eta[l][0] ;
      if (elastic) tempA -=  complex<double>(2*l+1);
      tempA *= Poly.LegendreP0(l,theta);
      

      //spin flip
      tempB = (eta[l][1]-eta[l][0])*Poly.LegendreP1(l,theta);

      if (proton) // include Coulomb phase shift
	{
	  tempA *= exp(complex<double>(0.,2.*Sigma[l]));
          tempB *= exp(complex<double>(0.,2.*Sigma[l]));
	}

      A += tempA;
      B += tempB;
    }


  if (proton == 0 && theta < .01) cout << A << " " << B << endl;

  A /= complex<double>(0.,2.*Kwave);
  B /= complex<double>(0.,2.*Kwave);  



  if (proton)
    {
      complex<double>Acoulomb = -gamma/2./Kwave/pow(sin(theta/2.),2)*
	exp(complex<double>(0.,2.*Sigma[0] 
        - 2.*gamma*log(sin(theta/2.))));
      A += Acoulomb;
    }

  //analysing power
  complex<double>temp = 2.*A*conj(B);
  AnalyzePower = -imag(temp)/(pow(abs(A),2)+pow(abs(B),2));

  SpinRotation = real(temp)/(pow(abs(A),2)+pow(abs(B),2));



  double xsec = pow(abs(A),2) + pow(abs(B),2); //xsec in fm**2/sr
  //return xsec/100.; // units barns/sr

  return xsec*10.;    // units of mb/sr
}
//*****************************************************************
  /**
   * returns the absorbtion cross section in mb
   * also calculates the array scat.SigmaAb - i.e. the cross section as 
   * a function of compound nucleus spin
   * the function integrateWave() must be executed first
   */
double scat::AbsorptionXsection()
{
  //l==0 contribution
  SigmaAb[0] = (1. - eta2[0][1]) *konst;
  SigmaAb_pos[0] = SigmaAb[0];
  SigmaAb_neg[0] = 0.; 
  double tot = SigmaAb[0];
  //contribution from higher l waves
 
  for (int i=1;i<=lStop;i++)
    {
      l = i;
      int parity = 1-2*(i%2);
      double up = (double)(l+1)*(1.-eta2[l][1])*konst;
      double down = (double)l*(1.-eta2[l][0])*konst;
      SigmaAb[i-1] += down;
      SigmaAb[i] = up;
      if (parity == 1)
	{
	  SigmaAb_pos[i-1] += down;
	  SigmaAb_pos[i] = up;
	  SigmaAb_neg[i] = 0.;
	}
      else
	{
	  SigmaAb_neg[i-1] += down;
	  SigmaAb_neg[i] = up;
	  SigmaAb_pos[i] = 0.;
	}
      tot += up + down;
    }
  return tot;
  
}
//*************************************************************
  /**
   * return the total Elastic cross section, 
   * infinite for protons (don't try), finite for neutrons in mb
   * The function integrateWave() must be called first
   */

double scat::ElasticXsection()
{
  complex<double> one(1.,0.);
  //l==0 contribution
  double tot = pow(abs(one-eta[0][1]),2);
  //contribution from higher l waves
 
  for (int i=1;i<=lStop;i++)
    {
      tot += (double)(i+1)*pow(abs(one-eta[i][1]),2) + 
	(double)i*pow(abs(one-eta[i][0]),2);
    }
  return tot*konst;
  
} 


//*************************************************************
  /**
   * returns the total cross section, infinite for protons (don't try),
   * finite for neutrons in mb
   * The function integrateWave() must be called first
   */
double scat::TotXsection()
{
  //l==0 contribution
  double tot = 1-real(eta[0][1]);
  //contribution from higher l waves
 
  for (int i=1;i<=lStop;i++)
    {
      tot += (double)(i+1)*(1.-real(eta[i][1])) +
      (double)i*(1.-real(eta[i][0]));
    }
  return tot*konst*2.;
  
} 

//**************************************************************
  /**
   * returns the Rutherford differential cross section in mb
\param theta is the center-of-mass scattering angle in radians
  */
double scat::Rutherford(double theta)
{
  return 10.* pow(gamma,2)/4./Kwave2/pow(sin(theta/2.),4);
}
//*********************************************************************
  /**
   * After integration the wavefunction to rStop, this function returns 
   *Wavefunction magnitude and its derivative
   \param l0 is the orbital angular momentum of the nucleon
   \param j0 is the total angular momentum of the nucleon
   */
valarray<double> scat:: IntegrateBound( double j0, int l0)
{
 
  l = l0;
  j = j0;

  // initialize the coupled differential equation solver
  initMerson(.0001,.00001,.1,0);

  LdotSigma = j*(j+1.) - double(l*(l+1)) - 0.5*1.5 ;
  valarray<double> WaveFunct(2);

  // find initial and matching wave functions
  //initialise wave function
  // potential at start
  double Vstart = RealPotential(rStart);

  //derivative of potential at start
  double dVstart = (RealPotential(rStart+.01) - Vstart)/0.01;

  // initialize wavefunctions
  double fact = pow(rStart,l+3)/2./(double)(2*l+3);
  WaveFunct[0] = pow(rStart,l+1) - muhbar*(energyCM-Vstart)*fact; 
  // derivative of wavefunction
  fact = (double)(l+3)*pow(rStart,l+2)/2./double(2*l+3);
  WaveFunct[1] = (double)(l+1)*pow(rStart,l) 
                      - muhbar*(energyCM-Vstart)*fact;
  fact = pow(rStart,l+3)*muhbar/2./(double)(2*l+3);
  WaveFunct[1] += dVstart*fact;

  //integrate wavefunctions out to matching radius
  solveMerson(&WaveFunct,rStart,rStop);

  if (ok == 0) 
     {
       cout << "j= " << j << " l= " << l << " Ecm= " << energyCM << endl;
       valarray<double> finish(0);
       return finish;
     }

  return WaveFunct;
}
//*****************************************************************
  /**
   *calculates the Wavefunction and stores it in an array WaveArray
   */

void scat::GetWaveFunctionArray(double j0, int l0,
     double& derivative)
{

  j = j0;
  l = l0;
  LdotSigma = j*(j+1.) - double(l*(l+1)) - 0.5*1.5 ;


  WaveArray[0] = 0.;

 // initialize the coupled differential equation solver
  initMerson(.001,.00001,.1,0);
  valarray<double> WaveFunct(2);



  //initialise wave function
 // potential at start
  double Vstart = RealPotential(rStart);
     
  //derivative of potential at start
  double dVstart = (RealPotential(rStart+.01) - Vstart)/0.01;

  // initialize wavefunctions
  double fact = pow(rStart,l+3)/2./(double)(2*l+3);
  WaveFunct[0] = pow(rStart,l+1) - muhbar*(energyCM-Vstart)*fact; 
  // derivative of wavefunction
  fact = (double)(l+3)*pow(rStart,l+2)/2./double(2*l+3);
  WaveFunct[1] = (double)(l+1)*pow(rStart,l) 
                      - muhbar*(energyCM-Vstart)*fact;
  fact = pow(rStart,l+3)*muhbar/2./(double)(2*l+3);
  WaveFunct[1] += dVstart*fact;


  //integrate wavefunctions out to matching radius, but storing value
  // at intervales of deltaR

  double r1 = rStart;
  double r2 = r1+deltaR;
  WaveArray[0] = WaveFunct[0];


  for (int i=1;i<mWave;i++)
    {
     solveMerson(&WaveFunct,r1,r2);
     if (WaveFunct.size() == 0)
       {
	 cout << "problem in GetWaveFunctionArray " << endl;
	 abort();
       }
     WaveArray[i] = WaveFunct[0];
     WaveFunct = WaveFunct;



     r1 = r2;
     r2 += deltaR;
    }

  derivative = WaveFunct[1];
}
//*************************************************************
  /**
   * Normalizes the wavefunction
   */

void scat::normalizeWaveFunction(double E, double Efermi)
{

      double SumBar = 0.;
      double r = rStart;
      for(int kk=0;kk<nWave;kk++)
	{
	  WaveBar[kk] = WaveArray[kk]*sqrt(eMassHF(r));
          double delta = deltaR;
          if (kk == 0 || kk == nWave-1) delta /=2.;
	  SumBar += pow(WaveBar[kk],2)*delta;
          r += deltaR;
	}


      ANC = pow(ANC,2)/SumBar; 
      // normalize wave function and calcalates wave function averages

      AverageMassEff = 0.;
      AverageInvMassBar = 0.;
      AverageW = 0.;
      Rrms = 0.; // RMS radius
      Occupation = 0.; // occupation probability

      r = rStart;
      for (int kk=0;kk<nWave;kk++)
	{
          double delta = deltaR;
	  if (kk == 0 || kk == nWave-1) delta /= 2.;
	  WaveBar[kk] /= sqrt(SumBar);
	  WaveArray[kk] /= sqrt(SumBar);

          double MassEff = eMass(r);


          double MassHF = eMassHF(r);
          double MassBar = MassEff/MassHF;
          double WW = Volume->ImaginaryPot(r) 
             + SurfaceLow->ImaginaryPot(r) 
             + SurfaceHigh->ImaginaryPot(r) 
             + SpinOrbit->ImaginaryPotential(r,LdotSigma);
          AverageMassEff += pow(WaveBar[kk],2)*delta*MassEff;
          AverageInvMassBar +=  pow(WaveBar[kk],2)*delta/MassBar;
          AverageW +=   pow(WaveBar[kk],2)*delta*WW;
          Rrms += pow(WaveBar[kk],2)*delta*pow(r,2);


          if (E > Efermi) // Particles
	    {
              Occupation += -pow(WaveBar[kk],2)/MassHF*delta*DerDelta(r);
	    }
	  else //Holes
	    {
	      Occupation += pow(WaveBar[kk],2)*delta*
                            (1.+1./MassHF*DerDelta(r));

	    }
          r += deltaR;

          // leave as sqaure of wavefunction
          //WaveBar[kk] = pow(WaveBar[kk],2);
	}
      Width = -2.*AverageW/AverageMassEff;
      SpectFactor = AverageInvMassBar;
      ANC *= SpectFactor;
      Rrms = sqrt(Rrms);
}

//********************************************************************
  /**
   *returns the effective mass relative to nucleon mass
\param r is the radial distance in fm
   */
double scat::eMass(double r)
{
  return 1. - HartreeFock->DerivativePotential(r) 
   - Volume->DerivativeDispersive(r) 
    - SurfaceLow->DerivativeDispersive(r)
    - SurfaceHigh->DerivativeDispersive(r);

}
//*******************************************************************
  /**
   * returns the Hartree-Fock effective mass relative to the mucleon mass
   * This is the momentum-dependent effective mass
\param r is the radial distance in fm
   */
double scat::eMassHF(double r)
{
  return  1. - HartreeFock->DerivativePotential(r);

}
//*******************************************************************
  /**
   * returns the frequency-dependent effective mass relative to the 
   * nucleon mass
   \param r is the radial distance in fm
  */
double scat::eMassBar(double r)
{

  return 1. - (Volume->DerivativeDispersive(r) 
    +          SurfaceLow->DerivativeDispersive(r)
	       + SurfaceHigh->DerivativeDispersive(r))/eMass(r);


}
//*******************************************************************
  /**
   * returns the derivative of the dispersive corrections
   \param r is the radial distance in fm
   */

double scat::DerDelta (double r)
{
  return SurfaceLow->DerivativeParticleHole(r) 
    + SurfaceHigh->DerivativeParticleHole(r) 
    + Volume->DerivativeParticleHole(r)
    + SpinOrbit->DerivativeParticleHole(r,LdotSigma);
}
//*********************************************************************
  /**
   * initialize integration of wavefunction
   */

void scat::initIntegration()
{
  rStart = .05;
  rStop = 12.;
  mWave = 60;
  nWave = mWave + 120;
  lMax = 100;
  llMax = lMax + 1;
// steps for wavefunction array in fm
  deltaR = (rStop-rStart)/(double)(mWave-1); 
  WaveArray = new double[nWave];
  WaveBar = new double[nWave];
  WaveMom = new double[100];

  Sigma = new double[llMax];
  eta2 = new double * [llMax];
  eta = new complex<double> * [llMax];
  phaseShift = new double * [llMax];
  for (int i=0;i<llMax;i++)
    {
      eta2[i] = new double [2];
      eta[i] = new complex<double> [2];
      phaseShift[i] = new double[2];
    }
 
  //set up array of Log derivatives.
  Nl = 7;
  LogDerMin = -123;
  LogDerMax = (int)EminRel;
  NLogDer = LogDerMax - LogDerMin + 1;
  
  LogDerArray = new double *[Nl+1];
  for (int i=0;i<=Nl;i++) LogDerArray[i] = new double [NLogDer];


  //cout << " Making LogDer array" << endl;
  for (l = 0;l<=Nl;l++)
    {

     for (int i=LogDerMin;i<=LogDerMax;i++)
      {
        energyCM = (double) i;
        Kwave2 = kconstant*mu*energyCM;
        muhbar = kconstant*mu;
        Kwave = sqrt(abs(Kwave2)); //wave number in fm**-1
        gamma = mu*Z*Zp/Kwave/28.820;  // Sommerfeld parameter

        LogDerArray[l][i-LogDerMin] = exteriorLogDer(l);
        
      }
    }
  //cout << " Finished LogDer array" << endl;

  //make array for absorption cross section
  SigmaAb = new double[lMax+1];
  SigmaAb_pos = new double [lMax+1];
  SigmaAb_neg = new double [lMax+1];
}
//*************************************************************************
  /**
   *Destructor
   */
scat::~scat()
{

  if (flag == 0) return;
  delete [] WaveArray;
  delete [] WaveBar;
  delete [] WaveMom;
  delete [] Sigma;

  for (int i=0;i<llMax;i++)
    {
      delete [] eta2[i];
      delete [] eta[i];
      delete [] phaseShift[i];
    }
  delete [] eta2;
  delete [] eta;
  delete [] phaseShift;

  for (int i=0;i<=Nl;i++) delete [] LogDerArray[i];

  delete [] LogDerArray;

  delete [] SigmaAb;
  delete [] SigmaAb_pos;
  delete [] SigmaAb_neg;
}
//************************************************************************
  /**
   * calculates the log derivative of the exterior wavefunction
   \param l is the orbital angular momentum of the nucleon
  */
double scat::exteriorLogDer(int l)
{
  if (proton == 0)
    {
      sphericalB sph;
      if (energyCM < 0.) return sph.LogDer_k(l,rStop*Kwave)*Kwave;
      else if (energyCM >0.) return sph.LogDer_y(l,rStop*Kwave)*Kwave;
      else return -1.e32;
    }
  else
    {

      if (energyCM <= 0.)//from whittaker function for bound states
	{
	  whit whittaker(16);
	  double outwave,DoutwaveDr;
	  if (energyCM < 0.)
	    {
	      outwave = whittaker.getWaveFunction(gamma,l,rStop*Kwave);
	      DoutwaveDr = whittaker.derivative*Kwave;
	    }
	  else 
	    {
	      if (proton == 0) return 1.e32;
	      outwave = whittaker.ZeroEnergy( 2.*mu*Z*Zp/28.820*rStop,l);
	      DoutwaveDr = whittaker.derivative*2.*mu*Z*Zp/28.820;
	    }
	  return DoutwaveDr/outwave;
	}
      else  // from irregular Coulomb wave function for quasi-bound states
	{
	  coul Coulomb;
	  Coulomb.init(l,gamma,rStop*Kwave);
	  return Coulomb.dG/Coulomb.G*Kwave;
	}
    }
}
//************************************************************************
  /**
   * adds the exterior part of the wavefunction to the array WaveArray
   * containing the interior part.
   \param l is the orbital angular momentum of the nucleon
  */
void scat::exteriorWaveFunct(int l)
{
  ANC=0.;
  if (proton == 0)
    {
      sphericalB sph;
      for (int i=0;i<=nWave-mWave;i++)
	{
	 double r =rStop + (double)i*deltaR;
	 double outwave;
         if (energyCM <= 0.) outwave = sph.k(l,r*Kwave);
	 else outwave = sph.y(l,r*Kwave);
	 if (i==0) ANC = WaveArray[mWave-1]/outwave;
	 else WaveArray[mWave+i-1] = outwave*ANC;
	}
    }
  else
    {
      if (energyCM <= 0.)//from whittaker function for bound states
	{
	  whit whittaker(16);
	  for (int i=0;i<=nWave-mWave;i++)
	    {
	      double r =rStop + (double)i*deltaR;
	      double outwave;
	      if (energyCM < 0.) outwave = 
              whittaker.getWaveFunction(gamma,l,r*Kwave);

	      else outwave = whittaker.ZeroEnergy( 2.*mu*Z*Zp/28.820*r,l);
	      if( i==0) ANC = WaveArray[mWave-1]/outwave;
	      else WaveArray[mWave+i-1] = outwave*ANC;
	    }
	}
      else  // from irregular Coulomb wave function for quasi-bound states
	{
	  coul Coulomb;
	  for (int i=0;i<=nWave-mWave;i++)
	    {
	      double r =rStop + (double)i*deltaR;
	      Coulomb.init(l,gamma,r*Kwave);
	      if( i==0) ANC = WaveArray[mWave-1]/Coulomb.G;
	      else WaveArray[mWave+i-1] = Coulomb.G*ANC;
	    }

	}
    }
}
//************************************************************
void scat::list()
{

  for (int i=-10;i<10;i++)
    {
      if (i==0) continue;
      energyCM = (double) i/4.;
      Kwave2 = kconstant*mu*energyCM;
      muhbar = kconstant*mu;
      Kwave = sqrt(abs(Kwave2)); //wave number in fm**-1
      gamma = mu*Z*Zp/Kwave/28.820;  // Sommerfeld parameter

    }
}
//****************************************************************
double scat::LogDerDifference( double j, int l)
{
  //this would be fine if you didn't need the wavefunction
  //valarray<double> Waveff = IntegrateBound(j,l);


  //if wavefunction is needed must use this subroutine
  double Waveff[2];
  GetWaveFunctionArray(j, l,Waveff[1]);
  Waveff[0] = WaveArray[mWave-1];

  //log derivative of interior wave function
  //double LOgDerInterior = Waveff[1]/Waveff[0];

  //interpolate the exterior LogDerative
  double en = floor(energyCM);
  double delta = energyCM-en;
  int i = (int)en- LogDerMin;
  if (i < 0 || i > NLogDer) cout << " outside of LogDerArray, Ecm = " 
				<< energyCM << " i= " << i << endl;
  double LogDerExterior = LogDerArray[l][i];

  if (delta != 0.)
    {
      if (proton == 0 && fabs(energyCM)< 1.) 
	LogDerExterior = exteriorLogDer(l);

      else LogDerExterior += delta*(LogDerArray[l][i+1]-LogDerExterior);
    }

  // if (l == 4) cout << "ext = " << LogDerExterior << endl;

  //tied the difference in interior and exterior logderivatives, 
  //but is not very very useful
  //as it only changes very close to a bound state andf hence they easily
  //missed.
  //return LogDerInterior - LogDerExterior;

  //instead
  //normalise interior and exterior wavefunctions to have the same slope,
  //then take the different in the magnitude of these wavefunctions - 
  // this crosses zero at a eigenfunction. one can easily find the zeros.

  return Waveff[0] - Waveff[1]/LogDerExterior;
  
}
//**********************************************************************
  /**
   *returns the expectation value of the imaginary part of the mean field
   *used for spectral functions
   */
double scat::spectralWidth()
{
  double sum = 0.;
  double r = rStart;
      for (int kk=0;kk<nWave;kk++)
	{
          double delta = deltaR;
	  if (kk == 0 || kk == nWave-1) delta /= 2.;
	  double pot = Volume->ImaginaryPot(r) 
          + SurfaceLow->ImaginaryPot(r)
          + SurfaceHigh->ImaginaryPot(r)
          + SpinOrbit->ImaginaryPotential(r,LdotSigma);
	  sum += pot*pow(WaveBar[kk],2)*delta; 
          r += deltaR;
	}
  return sum;
}
//**************************************************************
  /**
   *returns the transmission coefficients
   *the function integrateWave() must be executed first
   \param l is the orbital angular momentum of the nucleon
   \param j is the total angular momentum of the nucleon
   */
double scat::TransCoef(int l, double j)
{
  if (l > lStop) return 0.;
  int jj = (int)j;
  if (jj == l)return (1.-eta2[l][1]);
  else return 1. - eta2[l][0];
}
//********************************************************************
//********************************************************************
  /**
   * calculates the integrated potentials and the RMS values of the 
   * potentials
   */

void scat::VIntegrals()
{
  double deltaR = .1;
  double sumReal = 0.;
  double sumImag = 0.;
  double sumSO = 0.;
  double sumRealR2 = 0.;
  double sumImagR2 = 0.;
  double sumSOR2 = 0.;
  for (int i = 0;i<200;i++)
    {
      double r = ((double)i+0.5)*deltaR;
      double potReal = HartreeFock->RealPotential(r)
      + Volume->DispersiveCorrection(r)
	+ SurfaceLow->DispersiveCorrection(r)
	+ SurfaceHigh->DispersiveCorrection(r);
      double potImag = Volume->ImaginaryPot(r) 
        + SurfaceLow->ImaginaryPot(r)
        + SurfaceHigh->ImaginaryPot(r);
      double potSO = SpinOrbit->RealPotential(r,1.); 

      sumReal += potReal*pow(r,2)*deltaR;
      sumImag += potImag*pow(r,2)*deltaR;
      sumSO   += potSO*pow(r,2)*deltaR;

      sumRealR2 += potReal*pow(r,4)*deltaR;
      sumImagR2 += potImag*pow(r,4)*deltaR;
      sumSOR2 += potSO*pow(r,4)*deltaR;

    }

  JReal = sumReal*4.*pi/A;
  JImag = sumImag*4*pi/A;
  JSO = sumSO*4.*pi/pow(A,1./3.);

  RrmsReal = sqrt(sumRealR2/sumReal);
  RrmsImag = sqrt(sumImagR2/sumImag);  
  RrmsSO = sqrt(sumSOR2/sumSO);  
}
//******************************************************************
  /**
   * calculates the wavefunction at the specified momentum
   \param momentum is the momentum
   */
double scat::MomWave(double momentum)
{
  double r = rStart;
  double sum = 0.;
  for (int i = 0;i<nWave;i++)
    {
      double omega = momentum*r*.0050677;
      sum += WaveBar[i]*sin(omega)*deltaR;
      r += deltaR;
    }
  return 4./sqrt(pi)/momentum*sum;
}
//****************************************************************
  /**
   * calcuates the wavefunction in momentum space
   */

void scat::GetMomWave()
{
  double deltaMom = 3.;
  double sum = 0.;
  for (int i=0;i<100;i++)
    {
      double momentum = (double)i*deltaMom+.5;
      WaveMom[i] = MomWave(momentum);
      sum += WaveMom[i]*pow(momentum,2)*deltaMom;
    }
  sum *= 4.*pi;
  for (int i=0;i<100;i++)
    {
     WaveMom[i] /= sum;
     //cout << (double)i*deltaMom + .5 << " " << pow(WaveMom[i],2) << endl;
    }

}
//****************************************************
  /**
   *determines the number of nodes in the wavefunction
   */
int scat::Nodes()
{
  int nodes = 0;
  for (int i = 1;i<mWave-2;i++)
    {
 
      //cout << i << " " << WaveArray[i] << endl;
      if (WaveArray[i]*WaveArray[i-1] <= 0.) nodes++;

    }
  return nodes;
}
//******************************************************
//*****************************************************************
  /**
   *returns the absorption xsec associated with the isovector surface 
   *imaginary potential
   */
double scat::isovectorAbsorb()
{


  absorbAll = 0.;  //absorption xsec for all imaginary component
  absorbSurface = 0.;  //absorption xsec for surface imaginary component
  absorbVolume = 0.;
  absorbStandard = 0.;
  absorbSO = 0.; //spinOrbit

  // initialize the coupled differential equation solver
  initMerson(.001,.00001,.1,0);

  // find initial and matching wave functions
  //initialise wave function
  double rhoStop = rStop*Kwave; 


  waves outWave(rhoStop,gamma,lMax);
  for (int i=0;i<lMax;i++) Sigma[i] = outWave.Sigma[i];

  valarray<double> WaveFunct(4);
  double Wave2Array[mWave];  //array containing the modulus squared of wavefunction

  l = 0;
  lStop = -1;
  for (;;) // loop over orbital angular momentum
    {
      for (int plusMinus = -1;plusMinus<=1;plusMinus+=2)// spin up or spin down
	{
          j = (float)l + (float)plusMinus*0.5;
          if (j < 0.) continue;
          LdotSigma = j*(j+1.) - double(l*(l+1)) - 0.5*1.5 ;


          // potential at start
          double Vstart = RealPotential(rStart);
          double Wstart = ImaginaryPotential(rStart);
          //derivative of potential at start
          double dVstart = (RealPotential(rStart+.01) - Vstart)/0.01;
          double dWstart = (ImaginaryPotential(rStart+.01) - Wstart)/0.01;

	  // initialize wavefunctions
          double fact = pow(rStart,l+3)/2./(double)(2*l+3);

          WaveFunct[0] = pow(rStart,l+1) 
                                - muhbar*(energyCM-Vstart)*fact; // real part
          WaveFunct[2] =  Wstart*muhbar*fact;              //imaginary part

          // derivative of wavefunction
          fact = (double)(l+3)*pow(rStart,l+2)/2./double(2*l+3);
          WaveFunct[1] = (double)(l+1)*pow(rStart,l) 
                           - muhbar*(energyCM-Vstart)*fact; // real
          WaveFunct[3] = muhbar*Wstart*fact;          // imaginary

          fact = muhbar*pow(rStart,l+3)/2./(double)(2*l+3);
          WaveFunct[1] += dVstart*fact;
          WaveFunct[3] += dWstart*fact;


          Wave2Array[0] = pow(WaveFunct[0],2)+pow(WaveFunct[2],2);          
          //integrate wavefunctions out to matching radius, but storing value
          // at intervales of deltaR
          double r1 = rStart;
          double r2 = r1+deltaR;
          for (int i=1;i<mWave;i++)
             {
              solveMerson(&WaveFunct,r1,r2);
              if (WaveFunct.size() == 0)
                {
	          cout << "problem in GetWaveFunctionArray " << endl;
	          abort();
                }
              Wave2Array[i] = pow(WaveFunct[0],2)+pow(WaveFunct[2],2);
              WaveFunct = WaveFunct;
              r1 = r2;
              r2 += deltaR;
             }

          if (ok == 0) 
	    {
              cout << "j= " << j << " l= " << l << " Ecm= " <<
		     energyCM << endl;
              return 0;
	    }

          //outWave gives derivates with respect to rho = Kwave*r
	  //but at this point WaveFunctOut have derivative with respect to r
          // so make them with respect to rho=kwave*r

	  WaveFunct[1] /= Kwave;
          WaveFunct[3] /= Kwave;


	  // match wave functions 
	  //real WaveFunct = AA*F + BB*G
	  double  BB = outWave.dF[l]*WaveFunct[0] 
                     - outWave.F[l]*WaveFunct[1];
	  double AA = -outWave.dG[l]*WaveFunct[0] 
                      + outWave.G[l]*WaveFunct[1];

	  // imaginary part => Wavefunct  = CC*F + DD*G
	  double DD = outWave.dF[l]*WaveFunct[2] 
                    - outWave.F[l]*WaveFunct[3];
	  double CC = -outWave.dG[l]*WaveFunct[2] 
                    + outWave.G[l]*WaveFunct[3];


	  double denominator = pow(AA+DD,2) + pow(CC-BB,2);
	  double etaReal = (pow(AA,2)-pow(DD,2) + pow(CC,2) 
                    - pow(BB,2))/denominator;
	  double etaImag = 2.*(AA*BB+CC*DD)/denominator;
          double eta2 = pow(etaReal,2) + pow(etaImag,2);
         

          //normalization factor
	  //complex<double> norm ((BB-CC)/2.,(AA+DD)/2.); 
          double norm = pow(BB-CC,2)+pow(AA+DD,2);

          double absorb_all = 0.;
          double absorb_S = 0.;
          double absorb_V = 0.;
          double absorb_SO = 0.;
	  double r = rStart;
          for (int i=0;i<mWave;i++)
	    {
              absorb_all += Wave2Array[i]/norm*deltaR*
                     ImaginaryPotential(r);
              absorb_S += Wave2Array[i]/norm*deltaR*
                      SurfaceLow->ImaginaryPot(r)*gammaRel;
              absorb_V += Wave2Array[i]/norm*deltaR*
		(Volume->ImaginaryPot(r)+SurfaceHigh->ImaginaryPot(r))
                *gammaRel;
              if (r> 0)absorb_SO += Wave2Array[i]/norm*deltaR*
		SpinOrbit->ImaginaryPotential(r,LdotSigma)*gammaRel;


	      /*
	      if (energyLab == 20. && l == 4 && plusMinus ==1)
		cout << r << " " << RealPotential(r) << endl;
	       
		cout << r << " " << SurfaceLow->ImaginaryPot(r) << 
		  " " << Wave2Array[i]/norm << endl;
		*/

	      r+= deltaR;


	    }

	  if (plusMinus == 1)
	    {
	      absorbAll += absorb_all*(double)(l+1);
	      absorbSurface += absorb_S*(double)(l+1);
	      absorbVolume += absorb_V*(double)(l+1);
              absorbSO += absorb_SO*(double)(l+1);
              absorbStandard += (1.-eta2)*konst*(double)(l+1);

              //if (energyLab == 20.) cout << 1 << " " << l << " "
              // << absorb_S*(double)(l+1) << endl; 

	    }
	  else 
	    {
	      absorbAll += absorb_all*(double)(l);
	      absorbSurface += absorb_S*(double)(l);
	      absorbVolume += absorb_V*(double)(l);
              absorbSO += absorb_SO*(double)(l);
              absorbStandard += (1.-eta2)*konst*(double)(l);

	     
              //if (energyLab == 20.) cout << 2 << " " << l << " "
              //    << absorb_S*(double)(l+1) << endl; 
	      

	    }

	}

      l++;
      if (l > lMax) break;
    }
  //the following are true classically
  //absorbIsovector *= -6.046*mu/pow(Kwave,3);
  //absorbIsoscaler *= -6.046*mu/pow(Kwave,3);
  //absorbVolume *= -6.046*mu/pow(Kwave,3);
  //absorbAll *= -6.046*mu/pow(Kwave,3);
  absorbSurface *= -4.*pi/Kwave/energyCM*10;
  absorbVolume *= -4.*pi/Kwave/energyCM*10;
  absorbAll *= -4.*pi/Kwave/energyCM*10;
  absorbSO *= -4.*pi/Kwave/energyCM*10;
 
   
  return absorbSurface;
}
  //*********************************************************************

  double scat::getSmatrixEikonal(double b)
    {
      LdotSigma = 0.;
      double const deltaZ = .1;
      double z = 0.;
      double maxR = 0.;
      double maxI = 0.;
      double sumR = 0.;
      double sumI = 0.;
      for (;;)
        {
          double r = sqrt(pow(b,2)+pow(z,2)); 
          double addR = HartreeFock->RealPotential(r)
	  + Volume->DispersiveCorrection(r)
          + SurfaceLow->DispersiveCorrection(r)
	  + SurfaceHigh->DispersiveCorrection(r);
          maxR = max(maxR,abs(addR));
          if (z == 0.) addR /= 2;
          sumR += addR*deltaZ;

          double addI = ImaginaryPotential(r)/gammaRel;
          maxI = max(maxI,abs(addI));
          if (z == 0.) addI /= 2;
          sumI += addI*deltaZ;  

          if (abs(addR) < maxR/100. && abs(addI) < maxI/100.) break;
          if (z > 30.) break;
          z+= deltaZ;
        }

      double x = (pow(energyLab/m0+1,2)-1.);
      double beta = sqrt(x/(1.+x));
     


     double hvelocity = beta*197.326;
     double chiReal = -2.*sumR/hvelocity;
     double chiImag = -2.*sumI/hvelocity;


     //add the Coulomb term
     if (Zp == 1) chiReal += 2.*Z*e2/hvelocity*log(Kwave*b);


     //double phase = chiReal;
     double mag = exp(-chiImag);

     Sreal = mag*cos(chiReal);
     Simag = mag*sin(chiReal);


     return 1. - pow(mag,2);
    }




