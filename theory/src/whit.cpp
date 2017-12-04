#include "../include/whit.h"

double const whit::pi = acos(-1.);
double const whit::EulerGamma = 0.57721566490153;
//***********************************************************************
  /**
   *Constructor
   \param n0 is the largest order that can be called
   */

whit::whit(int n0)
{
    n0 = 40;
  //make sure n0 is even
  if (n0%2 != 0)
    {
      n0 += 1;
      cout << " n set to " << n0 << endl;
    }
  n = n0;
  m = n/2;
  //allocate 2d array
  try{
     pd = new double * [n+2];
     for (int i=0;i<n+2;i++)
       {
         pd[i] = new double [m+2];
       }
     array = new double[n+2];
     }

  catch (bad_alloc &memoryAllocationException)
    {
      cout << "could not allocate whit memory "
	   << memoryAllocationException.what() << endl;
    
    }
}
//***********************************************************************
  /**
   *Destructor
   */
whit::~whit()
{
  for (int i=0;i<n+2;i++)
    {
      delete [] pd[i];
    }
  delete [] pd;
  delete [] array;
}


//***********************************************************************
  /**
   * returns the wave function for a pure coulomb potential
   * \param nu = sommerfeld parameter =Z*Zp*e2*Kwave/energyCM/2.
   * \param l = angular momentum value
   * \param x = Kwave*r where Kwave is the imaginary asymptotic wave number 
   * and r is the radius
   */
double whit::getWaveFunction(double nu, int l, double rho)
{
  //old way
  //double out = AsymptoticExpansion(-nu,l,2.*rho);
  //derivative *= 2.;
  //return out;

  return whittakerW(nu, l , rho);

}

//************************************************************************
  /**
   *          Computation of the Pochammer symbols          
   */

double  whit::Pochammer (double a, int n)
{
  if (n == 0) return 1.;
  double out = 1.;
  for(int i=0;i < n; i++)
    {
      out *= a + (double)i;
    }
  return out;  
}
//************************************************************************
  /**
   *     Computation of the Kummer function truncated at N terms    
   */
double whit::KummerN (int N, double a, double b, double z, double & Derivative)
{
  Derivative = 0.;
  if (N == 0) return 1.;
  double output = 1.;
  
  for (int i=1;i<=N;i++)
    {
      double r1 = Pochammer(a,i);
      double r2 = Pochammer(b,i);
      double r3 = Pochammer(1,i);
      output += r1/r2/r3*pow(z,i);
      Derivative += r1/r2/r3*(double)i*pow(z,i-1);
  }  
  return output;
}

//****************************************************************************
/**          Computation of the Kummer function
 */
double whit::Kummer (double a, double b, double z, double &Derivative)
{
  Derivative = 0.;
  double output = 1.;
  for (int i=1;i<=25;i++)
    {
      double r1 = Pochammer(a,i);
      double r2 = Pochammer(b,i);
      double r3 = Pochammer(1,i);
      output += r1/r2/r3*pow(z,i);
      Derivative += r1/r2/r3*(double)i*pow(z,i-1);
  }  
  return output;
}
//***********************************************************************
  /**     
   * Computation of 0F1 function
   */
double whit::oF1(double b, double z)
{
  double output = 1.;
  double scale = z/b;
  output += scale;
  int n=1;
  for (;;)
    {
     n++;
     scale *= z/(b+(double)(n-1))/((double)n);
     output += scale;
     if (abs(scale) < 1e-14) break;
    }
  return output;
}
//****************************************************************************
/**
 *    COMPUTATION OF THE PSI DIGAMMA FUNCTION
 *  calculate for large x using 6.3.18, then use recurrence relation 6.3.6
 * to get our value (Abramowitz and stegan)
 *
 \param x is the independent variable
 */
double whit::psi (double x)
{

  double const Zstart=25;
  double  z = Zstart; // value of x to calculate 6.3.18
  int p = 0;
  if (x < z)
    {
       p = (int) (z-x);
       z = (double)p + x; // want xx and x separated by an interger value
    }
  else z = x;

  //calculate phi using 6.3.18 for z
  double  phiLarge = log(z) - 1./2./z - 1./12./z/z + 1./120./pow(z,4) 
           - 1./252./pow(z,6);

  double dif = 0.; //difference between wanted phi and phiLarge
  if (x < Zstart) for (int i=1;i<=p;i++) dif += 1./((double)i-1.+x);

  //cout << " psi " << x << " " << phiLarge - dif << endl;
  return phiLarge - dif;
}
//****************************************************************************
/**
 * Computation of the Euler Gamma function (Striling like expansion)
 \param z is the independent variable
 */
double  whit::gamma (double z)
{
  if (z == 0) return 1;
  double x = 1.;
  
  while (z<21.)
    {
      x = x * z;
      z++;
    }

 double som =  1/12./z - 1./360./pow(z,3) + 1./42./30./pow(z,5)
   - 1./30./56./pow(z,7) + 5./66./90./pow(z,9) - 691/12./11./2730./pow(z,11)
   + 7./6./13./14./pow(z,13) - 3617./510./15./16./pow(z,15)
   + 43867./798./17./18./pow(z,17) - 174611./330./19./20./pow(z,19)
   + 854513./138./21./22./pow(z,21);

 double y = (z - .5) * log(z) - z + .5 * log(2 * 3.14159265359) + som;
 return  exp(y) / x;
}
//****************************************************************************
  /**
   *                 Computation of the Pade Approximants     
   *                       with Wynn algorithm    
   */
void whit::paddie() 
{

  //initialize
  double pinfinity = 9.9999999999999E199;
  for (int i=0;i<=n;i++) pd[i][0] = pinfinity;
  for (int i=0;i<=m;i++) pd[0][i] = 0.;
  pd[1][0] = 1.;
  for (int i=1;i<=n+1;i++) pd[i][1] = pd[i-1][1] + array[i-1];


  for (int j=1;j<=m;j++)
    {
      for (int i=1;i<=n+1-j;i++)
       {
         if (i > n+2 || j > m+2) cout << "1*******" << endl;
	 double q = pd[i][j];
         if (j-1 > m+1 || j-1 < 0) cout << "2*************" << endl;
         double ww = pd[i][j-1] - q;
         if (ww == 0.) ww = 1.e-200;
         if (i-1 > n+1 || i-1 < 0) cout << "3************" << endl;
         double nn = pd[i-1][j]-q;
         if (nn == 0.) nn = 1.e-200;
         if (i+1 > n+1 || i+1 < 0) cout << "4************" << endl;
         double ss = pd[i+1][j]-q;
         if (ss == 0.) ss = 1.e-200;


         double denomi = nn*ww + ss*ww - ss*nn;
         if (denomi == 0.) denomi = 1.e-200;
         if (j+1 > m+1 || j+1 < 0) cout << "5************" << endl;
	 //check for overflow
          if (log(abs(ss))+log(abs(nn))+log(abs(ww))-log(abs(denomi)) > 460.)
             pd[i][j+1] = sign(ss)*sign(nn)*sign(ww)*sign(denomi)*1.e200;
	 else pd[i][j+1] = q + ss*nn*ww/denomi; 


	  if (i > 4 && i == j+1 && abs(1.-pd[i][i]/pd[i-1][i-1]) < 1.e-6)
	    {
	      //no need to continue 
	      pd[m][m] = pd[i][i];
              return;
	    }
	  
       }
    }

  /*
  if (pd[m][m] > 1.e10)
    {
      cout << pd[m][m] << " " << pd[m-1][m-1] << " " << pd[m-2][m-2] << endl;
      cout << pd[m-1][m-1]-pd[m-2][m-2] << endl;
    }
  */

}
//*********************************************************************
//*************************************************************************
//**********      Computation of the function  U(A,B,Z)           *********
//**********         using equation 13-1-6, page 504              *********
//**********             ABRAMOVITZ and STEGUN                    *********
//**********     "Handbook of mathematical functions", Dover.     *********
//**********       (Ascending series and Pade Approximants)       *********
//**********   Microsoft QUICK BASIC 4.5 (Compatible computers)   *********
//*************************************************************************
//
//Paper : Computation of the Whittaker functions W(z) with series
//        expansions and Pade Approximants
//
//by      Charles de IZARRA, Olivier VALLEE, Jacqueline PICART
//        and Nguyet TRAN MINH
//
//        Centre Universitaire de Bourges, GREMI, UFR SCIENCES
//        UNIVERSITE D'ORLEANS,
//        rue Gaston Berger, BP 4043
//        18028 BOURGES CEDEX FRANCE
//        e-mail : izarra@centre.univ-orleans.fr
//
//*********************************************************************
double  whit::AscendingSeries(double nu, int l ,double x)
{

  //double  ecartace[60], ecartpade[60];

  double a = (double)l + 1. - nu;
  int nn = 2*l + 1;

  double Dtt1Dx;
  double tt1 = Kummer(a,nn+1,x,Dtt1Dx);

  Dtt1Dx = Dtt1Dx*log(x) + tt1/x;
  tt1 *= log(x); // FUNCTION M(a,b,z)*Ln(z)


  // KUMMER LIMITED AT N TERMS
  double Dtt2Dx;
  double tt2 = KummerN(nn-1,a-(double)nn,1-nn,x,Dtt2Dx); 

  double t7 = Pochammer(1.,nn-1); //Corresponds to (N-1)!
  double t8 = gamma(a);

  tt2 *= t7/t8/pow(x,nn);
  Dtt2Dx = Dtt2Dx*t7/t8/pow(x,nn) - (double)nn*tt2/x;

  double denomi = gamma(a-(double)nn);

  double factorn = Pochammer(1.,nn);
  double facteur = pow(-1.,nn+1)/factorn/denomi;

  //====== Computation of coefficients
  for (int i=0;i<=n;i++)
    {
      double r = (double)i;
      double t1 = Pochammer(a,i);
      double t2 = Pochammer((double)nn+1.,i);
      double t3 = Pochammer(1.,i);
      double t4 = psi(a+r);
      double t5 = psi(1.+r);
      double t6 = psi(1.+(double)nn+r);
      array[i] = t1/t2/t3*(t4-t5-t6)*pow(x,i);

    }


  paddie();
  double output = exp(-x/2)*pow(x,l+1)*(facteur*(pd[m][m]+tt1)+tt2); 

  //now for the derivative
  array[0] = 0.;
  for (int i=1;i<=n;i++) array[i] *= (double)i/x;
  paddie();

  derivative = (-0.5+(double)(l+1)/x)*output + 
  exp(-x/2)*pow(x,l+1)*(facteur*(pd[m][m]+Dtt1Dx)+Dtt2Dx);

  return output;
}
//***************************************************
double whit::AsymptoticExpansion(double nu, int l ,double x)
{
  if (l == 0)
    {
     array[0] = 1.;
     for (int i=1;i<=n;i++)
       {
         double g = Pochammer(-nu+(double)l+1.,i);
         double h = Pochammer(-nu-(double)l,i);
         double f = Pochammer(1.,i);
         array[i] = h*g/f/pow(-x,i);
       }
 
     paddie();
    }
  else pd[m][m] = 1.;
     double output = exp(-x/2.)*pow(x,nu)*pd[m][m];
  
  //now for the derivative
  if (l != 0)
    {
     array[0] = 0.;
     for (int i=1;i<=n;i++) array[i] *= -(double)i/x;

     paddie();
    }
  else pd[m][m] = 0.;
  derivative = -0.5*output + nu/x*output 
                      + pd[m][m]*exp(-x/2.)*pow(x,nu);

  return output;
}

//***********************************************************************
  /**
   *gives the Whittaker function or equivalently, G (irregular Coulomb wave
   *function) at zero energy. The normalization is abitrary
   *if sommerfield paramter eta = beta/k, then c=2*beta*r where r is the radius
   */
double whit::ZeroEnergy(double c, int l)
{
  int m =2*l + 1;
  double factorial = 1.;
  for (int i=1;i<=m;i++) factorial *= (double)i;
  double factorial2 = factorial*(double)(m+1);
  double phi = -EulerGamma;
  for (int i=1;i<=m;i++) phi += 1./(double)i;
  double phi2 = phi + 1./(double)(m+1);

  double term = 1./factorial;
  double dterm = 1./factorial2;

  double C = term;
  double dC = dterm;
  double D = -term*phi;
  double dD = -dterm*phi2;

  for (int i=1;i<25;i++)
    {
      term *= c/(double)i/(double)(m+i);
      dterm *= c/(double)i/(double)(m+1+i);
      phi += 1./(double)(m+i);
      phi2 += 1./(double)(m+1+i);
      C += term;
      dC += dterm;
      D -= term*phi;
      dD -= dterm*phi2;
    }

  
  double Dminus = 0.;
  double dDminus = 0.;
  if (m > 0)
    {
      term = factorial/c;
      if (m%2 == 1) term *= -1.;
      for (int i=0;i<m;i++)
       {
	 term *= -c/(double)(m-i);
         if (i > 0) 
	   {
             term /= (double)i;
             dterm = term/c*(double)i;
             dDminus += dterm;
	   }
         Dminus += term;
       }
    }


  phi =- EulerGamma;
  term = pow(c,m)/factorial;
  Dminus -= term*phi;
  dDminus -= (double)m*term*phi/c;
  for (int i=m+1;i<25;i++)
    {
      phi += 1./(double)(i-m);
      term *= c/(double)i/(double)(i-m);
      Dminus -= term*phi;
      dDminus -= (double)i*term*phi/c;
    }
 


  double out = pow(c,l+1)*(log(c)*C + D + Dminus/pow(c,m));
  derivative = (double)(l+1)/c*out + pow(c,l+1)*(C/c + log(c)*dC +
      dD - (double)m*Dminus/pow(c,m+1) + dDminus/pow(c,m));

  return out;
}
//***********************************************************
double whit::sign(double x)
{
  if (x > 0.) return 1.;
  else if (x < 0.) return 1.;
  else return 0.;
}

//****************************************************************
double whit::whittakerW(double nu, int l ,double x)
{
  double U=hypergeometricU((double)(l+1)+nu,(double)(2*l+2),2.*x);
  double U2=hypergeometricU((double)(l+2)+nu,(double)(2*l+3),2.*x);
    
  double  front = exp(-x)*pow(2.*x,l+1);
  double  w = U*front;
  derivative = ((double)l/x-1.)*w - 2.*front*((double)(l+1)+nu)*U2;
  return w;
}


//*********************************************************************
  /**
   * returns the confluent hypergeometric function U for small arguement x
   * from COMPUTATION OF SPECIAL FUNCTIONS 
   * by Shanjie Zhang and Jianming Jin,  John Wiley and Sons
   */

double whit::hypergeometricU(double a, double b, double x)
{
  double out1= 0.;
  double out = 0.;
  double aa = a-b+1.0;
  bool IL1 = a == floor(a) && a <= 0.0;
  bool IL2 = aa == floor(aa) && aa <= 0.0;
  bool IL3 = abs(a*(a-b+1.0))/x <= 3.5; 

  bool BL1 = x <= 5.0 || (x <= 10.0 && a <= 2.0);
  bool BL2 = (x > 5.0 && x <= 12.5)&&(a >= 1.0 && b >= a+4.0);
  bool BL3 = x > 12.5 && a >= 5.0 && b >= a+5.0;
  bool BN = b == floor(b) && b != 0.0;

  int id1 = -100;
  if (b != floor(b))
    {
      out = chgus(a,b,x);
      id1 = id;
      method=1;
      if (id1 >= 6) return out;
      out1=out;
    }

  if (IL1 || IL2 || IL3)
    {
      out = chgul(a,b,x);
      method=2;
      if (id >= 6) return out;
      if (id1 > id)
	{
	  method=1;
          id=id1;
          out=out1;
	}
    }

 if (a >= 0.0)
   {
     if (BN && (BL1 || BL2 || BL3))
       {
         out = chgubi(a,b,x);
         method=3;
       }
     else
       {
         out = chguit(a,b,x);
         method=4;
       }
   }
 else
   {
     if (b <= a)
       {
	 double a00=a;
	 double b00=b;
         a=a-b+1.0;
         b=2.0-b;
         out = chguit(a,b,x);
         out=pow(x,(1.0-b00))*out;
         a=a00;
         b=b00;
	 method=4;
       }
     else if(BN && (! IL1))
       {
         out = cchgubi(a,b,x);
	 method=3;
       }
   }

 //if (id < 6) cout << "No accurate result obtained" << endl;

 return out;
}
//************************************************************
double whit::chgus(double a, double b, double x)
{
  id=-100;
  double gammaA = gamma2(a);
  double gammaB = gamma2(b);
  double xg1=1.0+a-b;
  double gammaAB = gamma2(xg1);
  double xg2=2.0-b;
  double gammaB2= gamma2(xg2);
  double HU0=pi/sin(pi*b);
  double R1=HU0/(gammaAB*gammaB);
  double R2=HU0*pow(x,(1.0-b))/(gammaA*gammaB2);
  double HU=R1-R2;
  double HMAX=0.0;
  double HMIN=1.0E300;
  double H0 = 0.;
  for (int j=1;j<=150;j++)
    {  
      double dj = (double)j;
      R1=R1*(a+dj-1.0)/(dj*(b+dj-1.0))*x;
      R2=R2*(a-b+dj)/(dj*(1.0-b+dj))*x;
      HU=HU+R1-R2;
      double HUA=abs(HU);
      if (HUA > HMAX) HMAX=HUA;
      if (HUA < HMIN) HMIN=HUA;
      if (abs(HU-H0) < abs(HU)*1.0E-15) break;
      H0=HU;
    }
 double d1=log10(HMAX);
 double d2=0.;
 if(HMIN != 0.0) d2=log10(HMIN);
 id=15-(int)abs(d1-d2);
 return HU;
}
//**************************************************
double whit::chgul(double a, double b, double x)
{
  int NM = 0;
  double out;
  id=-100;
  double  aa=a-b+1.0;
  bool IL1= a == floor(a) && a <= 0.0;
  bool IL2= aa== floor(aa) && aa <= 0.0;
  if (IL1) NM=(int)abs(a);
  if (IL2) NM=(int)abs(aa);
  if (IL1 || IL2)
    {
      out=1.0;
      double R=1.0;
      for (int k=1;k<=NM;k++)
	{
	  double dk = (double)k;
	  R=-R*(a+dk-1.0)*(a-b+dk)/(dk*x);
	  out=out+R;
	}
      out /= pow(x,a);
      id=10;
    }
  else
    {
      out=1.0;
      double R=1.0;
      double RA;
      double R0=0.;
      for (int k=1;k<=25;k++)
	{
	  double dk = (double)k;
          R=-R*(a+dk-1.0)*(a-b+dk)/(dk*x);
          RA=abs(R);
          if (k > 5 && RA >= R0 || RA < 1.0E-15) break;
          R0=RA;
          out=out+R;
	}
      id=(int)abs(log10(RA));
      out/=pow(x,a);
    }
  return out;
}
//*************************************************************************
double whit::chgubi(double a, double b, double x)
{
  double const EL = 0.5772156649015329;
  id=-100;
  int N=(int)abs(b-1);
  double RN1=1.0;
  double RN=1.0;
  for (int j=1;j<=N;j++)
    {
      RN *= (double)j;
      if (j == N-1) RN1=RN;
    }

  double PS = psi(a);
  double GA = gamma2(a);
  double A0,A1,A2,GA1,UA,UB;
  if (b > 0.0)
    {
      A0=a;
      A1=a-N;
      A2=A1;
      GA1 = gamma2(A1);
      UA=pow(-1.,N-1)/(RN*GA1);
      UB=RN1/GA/pow(x,N);
     }
  else
    {
      A0=a+N;
      A1=A0;
      A2=a;
      GA1= gamma2(A1);
      UA=pow(-1.,N-1)/(RN*GA)*pow(x,N);
      UB=RN1/GA1;
    }
  double HM1=1.0;
  double R=1.0;
  double HMAX=0.0;
  double HMIN=1.0E300;
  double H0 = 0.;
  for (int k=1;k<=150;k++)
    {
      double dk = (double)k;
      R=R*(A0+dk-1.0)*x/((N+dk)*dk);
      HM1=HM1+R;
      double HU1=abs(HM1);
      if (HU1 > HMAX) HMAX=HU1;
      if (HU1 < HMIN) HMIN=HU1;
      if (abs(HM1-H0) < abs(HM1)*1.0E-15) break;
      H0=HM1;
    }

  double DA2 = 0;
  double DA1=log10(HMAX);
  if (HMIN != 0.0) DA2=log10(HMIN);
  id=15-(int)abs(DA1-DA2);
  HM1=HM1*log(x);
  double S0=0.0;
  for (int m=1;m<N;m++)
    {
      double dm = (double)m;
      if (b >= 0.0) S0=S0-1.0/dm;
      if (b < 0.0) S0=S0+(1.0-a)/(dm*(a+dm-1.0));
    }

  double HM2=PS+2.0*EL+S0;
  R=1.0;
  HMAX=0.0;
  HMIN=1.0E300;
  for (int k=1;k<=150;k++)
    {
      double S1=0.0;
      double S2=0.0;

      if (N > 0.0)
	{
	  for (int m=1;m<=k;m++)
	    {
	      double dm = (double)m;
	      S1=S1-(dm+2.0*a-2.0)/(dm*(dm+a-1.0));
	    }
	  for (int m=1;m<=N;m++)S2 += 1.0/double(k+m);
	}
      else
	{
 	  for (int m=1;m<=k+N;m++)
	    {
	      double dm = (double)m;
              S1=S1+(1.0-a)/(dm*(dm+a-1.0));
	    }
	  for (int m=1;m<=k;m++) S2=S2+1.0/(double)m;
	}

      double HW=2.0*EL+PS+S1-S2;
      R *= (A0+(double)k-1.0)*x/(double)((N+k)*k);
      HM2=HM2+R*HW;
      double HU2=abs(HM2);
      if (HU2 > HMAX) HMAX=HU2;
      if (HU2 < HMIN) HMIN=HU2;
      if (abs((HM2-H0)/HM2) < 1.0E-15) break;
      H0=HM2;
    }

  double DB1=log10(HMAX);
  double DB2 = 0.;
  if (HMIN != 0.0) DB2=log10(HMIN);
  int ID1=15-(int)abs(DB1-DB2);
  if (ID1 < id) id=ID1;
  double HM3=1.0;
  if (N == 0) HM3=0.0;
  R=1.0;
  for (int k=1;k<=N-1;k++)
    {
      R=R*(A2+(double)k-1.0)/(double)((k-N)*k)*x;
      HM3 += R;
    }

  double SA=UA*(HM1+HM2);
  double SB=UB*HM3;
  double HU=SA+SB;
  int ID2 = 0;
  if (SA != 0.0) ID1=(int)(log10(abs(SA)));
  if (HU != 0.0) ID2=(int)(log10(abs(HU)));
  if (SA*SB < 0.0) id=id-abs(ID1-ID2);
  return HU;
}
//***************************************************************

double whit::chguit(double a, double b, double x)
{
  double t[30]={
   .259597723012478E-1, .778093339495366E-1, .129449135396945, 
   .180739964873425,.231543551376029, .281722937423262,
   .331142848268448, .379670056576798,.427173741583078, 
   .473525841761707,.518601400058570, .562278900753945,
   .604440597048510, .644972828489477,.683766327381356, 
   .720716513355730,.755723775306586, .788693739932264,
   .819537526162146, .848171984785930,.874519922646898,
   .898510310810046, .920078476177628, .939166276116423,
   .955722255839996, .969701788765053,.981067201752598, 
   .989787895222222, .995840525118838, .999210123227436};

  double  w[30]={
   .519078776312206E-1, .517679431749102E-1,
   .514884515009810E-1, .510701560698557E-1,
   .505141845325094E-1, .498220356905502E-1,
   .489955754557568E-1, .480370318199712E-1,
   .469489888489122E-1, .457343797161145E-1,
   .443964787957872E-1, .429388928359356E-1,
   .413655512355848E-1, .396806954523808E-1,
   .378888675692434E-1, .359948980510845E-1,
   .340038927249464E-1, .319212190192963E-1,
   .297524915007890E-1, .275035567499248E-1,
   .251804776215213E-1, .227895169439978E-1,
   .203371207294572E-1, .178299010142074E-1,
   .152746185967848E-1, .126781664768159E-1,
   .100475571822880E-1, .738993116334531E-2,
   .471272992695363E-2, .202681196887362E-2};

  id=7;
  double A1=a-1.0;
  double B1=b-a-1.0;
  double C=12.0/x;
  double out0 = 0.;
  double out1 = 0.;
  double out2 = 0.;
  for (int m=10;m<=100;m+=5)
    {
      out1=0.0;
      double G=0.5*C/(double)m;
      double D=G;
      for (int j=1;j<=m;j++)
	{
          double S=0.0;
          for (int k=1;k<=30;k++)
	    {
	      double T1=D+G*t[k-1];
              double T2=D-G*t[k-1];
              double F1=exp(-x*T1)*pow(T1,A1)*pow(1.0+T1,B1);
              double F2=exp(-x*T2)*pow(T2,A1)*pow(1.0+T2,B1);
              S += w[k-1]*(F1+F2);
	    }
          out1=out1+S*G;
          D=D+2.0*G;
	}
      if (abs(1.0-out0/out1) < 1.0E-7) break;
      out0=out1;
    }
  double GA = gamma2(a);
  out1=out1/GA;
  for (int m=2;m<=10;m+=2)
    {
      out2=0.0;
      double G=0.5/(double)m;
      double D=G;
      for (int j=1;j<=m;j++)
	{
	  double S=0.0;
          for ( int k=1;k<=30;k++)
	    {
	      double T1=D+G*t[k-1];
              double T2=D-G*t[k-1];
              double T3=C/(1.0-T1);
              double T4=C/(1.0-T2);
              double F1=T3*T3/C*exp(-x*T3)*pow(T3,A1)*pow(1.0+T3,B1);
              double F2=T4*T4/C*exp(-x*T4)*pow(T4,A1)*pow(1.0+T4,B1);
              S=S+w[k-1]*(F1+F2);
	    }
          out2=out2+S*G;
          D=D+2.0*G;
	}
      if (abs(1.0-out0/out2) < 1.0E-7) break;
      out0=out2;
    }
  GA = gamma2(a);
  out2=out2/GA;
  return out1+out2;
}

//*********************************************************************
double whit::cchgubi(double a, double b, double x)
{

  id=-100;
  int N=(int)(b-1.);
  double RN1=1.0;
  double RN=1.0;
  for (int j=1;j<=N;j++)
    {
      RN *= (double)j;
      if (j == N-1) RN1=RN;
    }
  double PS = psi(a);
  double GA = gamma2(a);
  double A0,A1,A2,GA1,UA,UB;
  if (b > 0.)
    {
      A0=a;
      A1=a-(double)N;
      A2=A1;
      GA1 = gamma2(A1);
      UA=pow(-1.,N-1)/(RN*GA1);
      UB=RN1/GA/pow(x,N);
    }
  else 
    {
      A0=a+(double)N;
      A1=A0;
      A2=a;
      GA1 = gamma2(A1);
      UA=pow(-1.,N-1)/(RN*GA)*pow(x,N);
      UB=RN1/GA1;
    }
  double HM1=1.0;
  double R=1.0;
  double HMAX=0.0;
  double HMIN=1.0E+300;
  double H0 = 0.;
  for (int k=1;k<=150;k++)
    {
      double dk = (double)k;
      R=R*(A0+dk-1.0)*x/(double)((N+k)*k);
      HM1=HM1+R;
      double HU1=abs(HM1);
      if (HU1 > HMAX) HMAX=HU1;
      if (HU1 < HMIN) HMIN=HU1;
      if (abs(HM1-H0) < abs(HM1)*1.0E-15) break;
      H0=HM1;
    }
  double DA1=log10(HMAX);
  double DA2 = 0.;
  if (HMIN != 0.0) DA2=log10(HMIN);
  id=(int)(15.-abs(DA1-DA2));
  HM1=HM1*log(x);
  double S0=0.0;
  for (int m=1;m<=N;m++)
    {
     double dm = (double)m;
     if (b >= 0.0) S0=S0-1.0/dm;
     if (b < 0.0) S0=S0+(1.0-a)/(dm*(a+dm-1.0));
    }
  double const EL = 0.5772156649015329;

  double HM2=PS+2.0*EL+S0;
  R=1.0;
  HMAX=0.0;
  HMIN=1.0E+300;
  for (int k=1;k<=300;k++)
    {
      double S1=0.0;
      double S2=0.0;
      if (b > 0.0)
	{
	  for (int m=1;m<=k;m++)
	    {
	      double dm = (double)m;
              S1=S1-(dm+2.0*a-2.0)/(dm*(dm+a-1.0));
	    }
          for (int m=1;m<=N;m++)S2=S2+1.0/(double)(k+m);
	}
      else
	{
          for (int m=1;m<=k+N;m++)
	    {
	      double dm = (double)m;
	      S1 += (1.0-a)/(dm*(dm+a-1.));
	    }
	  for (int m=1;m<=k;m++)S2 +=1.0/(double)m;
	}

      double HW=2.0*EL+PS+S1-S2;
      R=R*(A0+(double)k-1.0)*x/(double)((N+k)*k);
      HM2=HM2+R*HW;
      double HU2=abs(HM2);
      if (HU2 > HMAX) HMAX=HU2;
      if (HU2 < HMIN) HMIN=HU2;
      if (abs((HM2-H0)/HM2) < 1.0E-15) break;
      H0=HM2;
    }

  double DB1=log10(HMAX);
  double DB2 = 0.;
  if (HMIN != 0.0) DB2=log10(HMIN);
  int id1=(int)(15.-abs(DB1-DB2));
  if (id1 < id) id=id1;
  double HM3=1.0;
  if (N == 0) HM3=0.0;
  R=1.0;
  for (int k =1;k<N;k++)
    {
      R=R*(A2+(double)k-1.0)/(double)((k-N)*k)*x;
      HM3 += R;
    }
  double SA=UA*(HM1+HM2);
  double SB=UB*HM3;
  double out=SA+SB;
  int id2 = 0;
  if (SA != 0.0) id1=(int)log10(abs(SA));
  if (out != 0.0) id2=(int)log10(abs(out));
  if (SA*SB < 0.0) id -= abs(id1-id2);
  return out;
}

//*********************************************************************
double whit::gamma2(double x)
{ 
  double g[26]={
   1.0,0.5772156649015329,-0.6558780715202538, -0.420026350340952E-1,
   0.1665386113822915,-.421977345555443E-1,-.96219715278770E-2,
   .72189432466630E-2,-.11651675918591E-2, -.2152416741149E-3,
   .1280502823882E-3, -.201348547807E-4,-.12504934821E-5,
   .11330272320E-5,-.2056338417E-6, .61160950E-8,.50020075E-8,
   -.11812746E-8,.1043427E-9,.77823E-11,-.36968E-11, .51E-12,
   -.206E-13, -.54E-14, .14E-14, .1E-15};

  double GA,Z;
  double R = 0.;
  if (x == floor(x))
    {
     if (x > 0.0)
       {
         GA=1.0;
         int M1= (int)x-1;
         for (int k=2;k<=M1;k++) GA *= (double)M1;
       }
     else GA = 1.E300;
    }
  else
    {
      if (abs(x) > 1.)
	{
          Z=abs(x);
          int M=(int)Z;
          R=1.0;
	  for (int k=1;k<=M;k++) R *= (Z-(double)k);
          Z=Z-(double)M;
	}
      else Z = x;
    
     double GR = g[25];
     for (int k=25;k>=1;k--) GR = GR*Z+g[k-1];
     GA=1.0/(GR*Z);
     if(abs(x) > 1.0)
       {
         GA=GA*R;
         if (x < 0.0) GA=-pi/(x*GA*sin(pi*x));
       }
    }
  return GA;
}
