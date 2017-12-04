#include "../include/compound.h"
#include "../include/compoundException.h"

#include <iostream>
#include <fstream>
#include <algorithm>

using namespace std;

int const compound::dim = 32000;
compound::compound()
{}
//******************************************************
compound::compound(string* name)
{
  sigmaAbsorption = 0.;
  sigmaCompoundElastic  = 0.;
  string Name = *name+".CN";  
  ifstream CNfile(Name.c_str());
  if (CNfile.fail())
    {
      cout << "could not open file " << Name << endl;
      okay = 0;
      return;
    }
  else cout << Name << " opened" << endl;
  okay = 1;
  CNfile >> elasticState; //file name for levels with same 
                          //entrance and exit chennels
  CNfile >> elasticTrans; //file name for transmission coef.
  CNfile >> otherState;   //file name for levels from competing exit channel
  CNfile >> otherTrans;   //file name for transmission coeff.
  CNfile >> Qother;       // Qvalue to ground state of  competeing channel
  CNfile.close();  

  elastic = new channel(&elasticState,&elasticTrans,0.);

  other = new channel(&otherState,&otherTrans,Qother);

  Ntll = min(elastic->Ntll,other->Ntll);
  

  xsec = new double [Ntll];
  xsec0 = new double [Ntll];
  xsec1 = new double [Ntll];


  // normalization factors for Legendre poly
  norm0 = new double [Ntll];
  norm1 = new double [Ntll];

  for (int i=0;i<Ntll;i++)
    {
      double I = (double)i;
      norm0[i] = (I+0.5)/2./3.14159;
      if (i > 0) norm1[i] = norm0[i]/I/(I+1.);
    }

  T = new double[dim]; //array of transmission coefficients for each allowed channel
  nu = new double[dim]; //array of degrees of freedom for each allowed channel
  mult = new double[dim]; //multiplicity of channel - for continuum level density
  
}
//******************************************************
compound::~compound()
{
  cout << "destroying compound" << endl;
  if (okay == 0) return;
  delete elastic;
  delete other;
  delete [] xsec;
  delete [] xsec0;
  delete [] xsec1;
  delete [] norm0;
  delete [] norm1;
  delete [] T;
  delete [] nu;
  delete [] mult;
}
//*******************************************************
double compound::decayElastic(level CN)
{
  Nch = 0;
  return decay(CN,elastic);
}
//*********************************************************
//finds all possible decay channels 
double compound::decay(level CN,channel *Chan)
{

  if (okay ==0) return 0.;


  //loop over levels of the daughter nucleus
  for (int i=0;i<Chan->Nlevel;i++)
    {
      double decayE = CN.energy - Chan->Level[i].energy;




      if (decayE <= 0.) break;

      if (Chan->Level[i].energy >Chan->threshold) break;


      //loop over the decay AM
      //find maximum and minimum AM carried away by evaporated particle
      double decayImin = abs(CN.spin - Chan->Level[i].spin);
      int Imin = (int)decayImin;
      double decayImax = CN.spin + Chan->Level[i].spin;
      int Imax = (int)decayImax;

      int deltaI = Imax - Imin;

      for (int j=0;j<=deltaI;j++)
	{

	  double decayI = decayImin + (double)j;

          //loop over the orbital angular momentum
          double decayLmin = abs(decayI - 0.5);
          double decayLmax = decayI + 0.5;
          int Lmin = (int)decayLmin;
          int Lmax = (int)decayLmax;
          if (decayLmin - (double)Lmin > 0.5) Lmin++;
          if (decayLmax - (double)Lmax > 0.5) Lmax++;
          int deltaL = Lmax - Lmin; 


	   for (int k=0;k<=deltaL;k++)
	      {
		int decayL = Lmin + k;
                int parity = 1-2*(decayL%2);

                if (parity*Chan->Level[i].parity != CN.parity) continue;
                if (decayL > Ntll-1) break;



                int aligned;
	        if (decayL == (int)decayI) aligned = 1;// decayL and spin aligned
	         else  aligned = 0;

		double trans = Chan->TransCoef(decayE,decayL,decayI);

       
                if (trans <= 0.) continue;
                T[Nch] = trans;
                mult[Nch] = 1.;

		T_tot += trans;

		// save location of elastic channels
		if (Chan->Level[i].energy == 0.)
		  {
                   Nelastic = Nch;
		   Lelastic = decayL;
                   Aligned = aligned;
		  }

          
                Nch++;
                if (Nch == dim) 
		  {
                    cout << "increase dimension of T array " << dim << endl;
		    compoundException cmpEx;
                    throw cmpEx;
		  }
	      
	      }
	}

    }
  //****************************************************
  //continuum
  if (Chan->threshold == 1000.) return T_tot;
  double const deltaE= 1.;
  double E = Chan->threshold+ deltaE/2.;



  for (;;)
    {
      double gamE = 0.;
      double decayE = CN.energy - E - Chan->Qvalue;
      if (decayE <= 0.) break;
      for (int jj=0;jj<20;jj++) // loop over spin of CN
	{
          double gamJ = 0.;
          double spin = (double)jj + 0.5;
          for (int iparity = -1;iparity<2;iparity+=2) //loop over parity of CN
	    {
	      //double rho = Chan->FermiGas(E,jj);
	      double rho = Chan->ld->getLD(E,jj);
              if (rho <= 0.) break;
              //loop over the decay AM
              //find maximum and minimum AM carried away by evaporated particle
                double decayImin = abs(CN.spin - spin);
               int Imin = (int)decayImin;
               double decayImax = CN.spin + (double)spin;
               int Imax = (int)decayImax;

               int deltaI = Imax - Imin;

               for (int j=0;j<=deltaI;j++)
		 {
   	           double decayI = decayImin + (double)j;
                   //loop over the orbital angular momentum
                   double decayLmin = abs(decayI - 0.5);
                   double decayLmax = decayI + 0.5;
                   int Lmin = (int)decayLmin;
                   int Lmax = (int)decayLmax;
                   if (decayLmin - (double)Lmin > 0.5) Lmin++;
                   if (decayLmax - (double)Lmax > 0.5) Lmax++;
                   int deltaL = Lmax - Lmin; 

	           for (int k=0;k<=deltaL;k++)
	              {
		        int decayL = Lmin + k;
                        int parity = 1-2*(decayL%2);
                        if (parity*iparity != CN.parity) continue;
                        if (decayL > Ntll-1) break;
                        
		        double trans = Chan->TransCoef(decayE,decayL,decayI);
                        if (trans <= 0.) break;
                        T[Nch] = trans;
                        mult[Nch] = rho*deltaE;

                        double gam = trans*rho*deltaE;
		        T_tot += gam;
                        gamE += gam;
                        gamJ += gam;
                        Nch ++;

		      }

		 }
              
	    }
          if (gamJ < T_tot/10000.) break; 

	}

      E += deltaE;

      if (gamE < 0.001*T_tot) break;
    }

  return T_tot;
} 
//************************************************************
//************************************************************************
//returns the angular distribution for compound elastic 
double compound::DifferentialXsectionCE(double theta)
{
  if (okay == 0) return 0.;
  legendre Poly(Ntll-1);


  double sum = 0.;
  for (int decayL = 0;decayL<Ntll;decayL++)
    {

      if (xsec[decayL] <= 0.) break;
      sum += xsec0[decayL] * pow(Poly.LegendreP0(decayL,theta),2)*
	norm0[decayL];
      //spin flip contribution
       if (decayL > 0) sum += xsec1[decayL]* 
       pow(Poly.LegendreP1(decayL,theta),2)*norm1[decayL];

    }

  return sum;
}
//******************************************************
//returns the total compound elastic cross section and prepares 
// for other compound elastic quantites
//Ex is the excitation energy

double compound::statistical(double kconst,  double Ex)
{


  if (okay == 0) return 0.;
  sigmaCompoundElastic = 0.;
  sigmaAbsorption = 0.;
  //first zero xsec arrays
  for (int i=0;i<Ntll;i++)
    {
      xsec[i] = 0.;
      xsec0[i] = 0.;
      xsec1[i] = 0.;
    }





  //loop over the j
  for (int i=0;i<Ntll;i++) //rjc
    {

      double TcPos = 0.;
      double TcNeg = 0.;
      //positive parity
      T_tot = 0.;
      Nch = 0; // zero number of open channels
      level CN1((float)i+0.5,1,Ex);
      Nelastic = -1;
      Lelastic = -1;
      if (elastic->Qvalue<Ex)decay(CN1,elastic);

      //at this point if Nch == 0, then there are no channels available
      //for decay
      // if Nch> 0 && Nelastic == -1, then the channel can decay, 
      //  but the ealstic channel is not open


      if (Nch > 0 && Nelastic != -1)
	{

          TcPos = T[Nelastic];

          if (TcPos > 0.)
	    {

             if (other->Qvalue<Ex)decay(CN1,other);

             //fluctuation corrections and increment xsection
             corrections();

	    }
	 }


      //negative parity
      T_tot = 0.;
      Nch = 0; // zero number of open channels
      level CN2((float)i+0.5,-1,Ex);
      Nelastic = -1;
      Lelastic = -1;
      if (elastic->Qvalue<Ex)decay(CN2,elastic);
      if (Nch > 0 && Nelastic != -1)
	{
         TcNeg = T[Nelastic];
         if (TcNeg+TcPos <= 0.) break;
         if (TcNeg > 0.)
	   {
             if (other->Qvalue<Ex)decay(CN2,other);

             //fluctuation corrections and increment xsection
             corrections();
	    }
	}


    }


      for (int i=0;i<Ntll;i++)  
	{
	  xsec[i] *= kconst;
	  xsec0[i] *= kconst;
	  xsec1[i] *= kconst;

          sigmaCompoundElastic += xsec[i];
	}
      sigmaAbsorption*= kconst;


  return sigmaCompoundElastic;
}




void compound::corrections()
{


  if (okay == 0) return;
  //at this point we have found all the available channels for compound nucleus decay and their
  //transmission coeeficinets are saved in the array T
  //calculate the fluctuations corrections. These depend whether the nucleon came out in the same
  //state as it went in, i.e, whether its spin was flipped or not.
  //there are G and C corrections, first the G factors

  //find degrees of freedom for each channel
  for (int i=0;i<Nch;i++) nu[i] = 1.78 + (pow(T[i],1.212) - 0.78)
       *exp(-0.228*T_tot);
   



  //integrate over t
  double deltat= .05;
  double t = 0;
  double Gaa = 1.;
  double Saa = 1.;

  if (Nch > 500)
    Gaa = 0.;
    {
     for(;;)
       {
	double gaa = 1.;
        for (int i=0;i<Nch;i++)
           {
            double x = 1. + 2.*T[i]/nu[i]/T_tot*t;
            double y = pow(x,-nu[i]/2.);
            if (y == 1.) continue;
            gaa *= pow(y,mult[i]);
            if (i == Nelastic)  gaa *= pow(x,-2);
            }
        if (t == 0) Gaa = gaa*deltat/2.;
        else Gaa += gaa*deltat;

        if (gaa < .01) break;


        t += deltat;
       }


     //now for the S factors
     Saa = 1.+2./nu[Nelastic];

    }



  double factL;
  if (Aligned) 
    {
      double trans = TransCoef(Lelastic,(double)Lelastic+0.5);
      factL = (double)(Lelastic+1)*trans;
    }
  else 
    {
      double trans = TransCoef(Lelastic,(double)Lelastic-0.5);
     factL = (double)Lelastic *trans;
    }

  double fact = factL*T[Nelastic]/T_tot*Gaa*Saa;

  if (Lelastic < 0) 
    {
      cout << Lelastic << " " << Nelastic << endl;
      abort();
    }
  xsec[Lelastic] += fact;



  if (Aligned) 
    {
      xsec0[Lelastic] += (double)(Lelastic+1)/(double)(2*Lelastic+1)*fact;
      xsec1[Lelastic] += (double)(Lelastic)/(double)(2*Lelastic+1)*fact;
     }
  else 
    {
      xsec1[Lelastic] += (double)(Lelastic+1)/(double)(2*Lelastic+1)*fact;
      xsec0[Lelastic] += (double)(Lelastic)/(double)(2*Lelastic+1)*fact;
    }

  sigmaAbsorption += factL;



}
