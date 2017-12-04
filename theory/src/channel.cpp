#include "../include/channel.h"
#include <string>

channel::channel(string*StateFile, string*TransFile, double Q0)
{
  //read in states for target
  Qvalue = Q0;
  string filename = *StateFile+".state";
  ifstream File(filename.c_str());

  if (File.fail()) 
    {
     cout << filename << " (class channel) is not open " << endl;
     abort();
    }

  File >> Nlevel;
  //cout << "************ " << Nlevel << " " << Qvalue << endl;


  Level = new level [Nlevel];

  int i2;
  double e1,e2;

  for (int i=0;i<Nlevel;i++)
    {
      File >> e1 >> i2 >> e2;
      if (File.eof()) 
        {
	  cout << "eof when reading levels in " << filename << endl;
          abort();
	}
      if (File.bad())
	{
	  cout << "problems when reading levels in " <<filename << endl;
          abort();
	}
      Level[i].spin = e1;
      Level[i].parity = i2;
      Level[i].energy = e2 + Qvalue;


    }


  //------ level density info 
  string name;
  File >> name >> threshold;

  if (name != string("threshold")) 
    {
      cout << "threshold not found in " << filename << endl;
      cout << name << endl;
      abort();
     threshold = 1000.;
    }
  if (File.eof()) threshold = 1000.;
  if (File.bad()) threshold = 1000.;

  //cout << name << " " << threshold << endl;

  if (threshold < 1000.)
    {
      File >> name >> backshift;
      File >> name >> littlea ;
      File >> name >> A;
      sigma2t = 0.0888*littlea*pow((double)A,2./3.);
    }
  else
    {
      littlea = 1.;
      A = 1;
      sigma2t = 1.;
    }

  File.close();
  File.clear();

  
  ld = new levelD(*StateFile);

  //read in transmission coefficients
  //first for spin up
  filename = *TransFile+"_up.tl";
  File.open(filename.c_str(),ios::in);
  if (File.fail()) 
    {
    cout << filename << " (class channel) not open" << endl;
    abort();
    }
  File >> Ntll >> Ntle;
  tl_up = new double * [Ntll];
  tl_down = new double * [Ntll];
  Earray = new double [Ntle];

  for (int i=0;i<Ntll;i++) 
    {
     tl_up[i] = new double [Ntle];
     tl_down[i] = new double [Ntle];
    }
  
  for (int i=0;i<Ntle;i++)
     {
       File >> e1;
       Earray[i] = e1;
       for (int j=0;j<Ntll;j++) 
	 {
	   File >> e1;
           if (e1 < 0.) 
	     {
               cout << "Transmission coeff negative in " << filename << endl;
	       e1 = 0.;
	     }
           tl_up[j][i] = e1;
	   if (File.eof() == 1) 
	     {
              cout << " eof in tl_up read (class channel) " << filename <<endl;
	      abort();
	     }
	 } 
     }
   //now for down
   File.close();
   File.clear();
   filename = *TransFile+"_down.tl";
   File.open(filename.c_str(),ios::in);
   if (File.fail()) 
     {
      cout << filename << " is not open (class channel)" << endl;
      abort();
     }
   File >> Ntll >> Ntle;
   for (int i=0;i<Ntle;i++)
     {
       File >> e1;
       if (e1 < 0.) 
	 {
           cout << "Transmission coeff negative in " << filename << endl;
	 }
       Earray[i] = e1;
       for (int j=0;j<Ntll;j++) 
	 {
	   File >> e1;
           tl_down[j][i] = e1;
	   if (File.eof() == 1) 
	     {
            cout << " eof in tl_down read (class channel) " << filename <<endl;
            abort();
	     }
	 } 
     }
   File.close();
}
//*********************************************************
channel::~channel()
{
  //cout << "destroying channel" << endl;
  delete [] Level;
  delete [] Earray;
  
  for (int i=0;i<Ntll;i++)
    {
      delete [] tl_up[i];
      delete [] tl_down[i];
    }
  
  delete [] tl_up;
  delete [] tl_down;
}
//*******************************************************
// returns the transmission coefficient
double channel::TransCoef(double ek, int l, double j)
{

 // find channel for Transmission coefficent interpolation
 int ie = 0;
 for (;;)
   {
     if (Earray[ie] > ek) break;
     if (ie == Ntle-1) break;
     ie++;
   }
 if (ie == 0) ie = 1;

  double *tl;
  if ((int)j == l) tl = tl_up[l];
  else tl = tl_down[l];

  double trans = (tl[ie] - tl[ie-1])
   /(Earray[ie]-Earray[ie-1])*(ek-Earray[ie-1]) + tl[ie-1];

  return trans;
}
//***************************************************
//returns the level density
double channel::FermiGas(double E,int J)
{
  double U = E - backshift;
  if (U <= 0.) return 0.; 
  double temp = sqrt(U/littlea);
  double sigma2 = sigma2t*temp;

  double fact = (double)(2*J+1)/24./sqrt(2.)/pow(sigma2,1.5)/pow(littlea,.25);
  fact *= exp(2.*sqrt(littlea*U))/pow(U+temp,1.25);
  fact *= exp(-(double)(J*(J+1))/2./sigma2);


  return fact;
}
