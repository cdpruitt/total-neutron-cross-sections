#include "../include/reaction.h"

double const reaction::pi = acos(-1.);


//*********************************************************************
  /**
   *constructor reads in data from file
   * Set falg to zero is class used only to gives optical model potential
   \param title for input data
   \param flag indicates integrating of wavefunctions takes place 
  */
reaction::reaction(string *title0,bool flag0)
{
  bprint = 1;
  for (int i=0;i<NlevelMax;i++) LevelTh[i].SpectFunction = NULL;
  Esf = NULL;
  NXdata = 0;
  NTotXdata = 0;
  DOM = 1;
  flag = flag0;
  //directory = "/home/bob/DOMA/";
  directory = "";


  string line;
  title = *title0;
  string filename(title + ".inp");
  if (flag) cout << filename << endl;
  ifstream file (filename.c_str(),ios::in);
  // if one cannot open file quit
  if (file.fail()) 
    {
      cout << "couldn't open data file " << filename << " " << title << endl;
      abort();
    }

  file >> Zp >> Z >> A >> Efermi;
  file >> Nfermi;
  if (Nfermi < 1 || Nfermi > 2) cout << "Nfermi not possible" << endl;
  file >> ValenceHole.N >> ValenceHole.j >> ValenceHole.l >> 
          ValenceHole.energy;
  if (Nfermi == 2) 
     file >> ValenceParticle.N >> ValenceParticle.j >> ValenceParticle.l
          >> ValenceParticle.energy;
  file >> gap  >> gapOther;
  gapMin = min(gap,gapOther);
  Wstart = gap/2. + gapMin;

  //Wstart = 0.9*gap;


  asymmetry = (A-2.*Z)/A;
  if (Zp == 1.) sign = 1;
  else if (Zp == 0.) sign = -1.;
  else cout << " unknown Zp " << endl;
  if (file.bad()) cout << "problem with input file" << endl;
  if (file.eof()) cout << "eof with input file" << endl;

  file.close();
  file.clear();


  //asymmetry
  AsyVolume = 1;  // asymmetry for volume
  //alphaVolume = 1.65;
  //EaVolume = 140.;
  //EaVolume = 60.;
  

  //efective mass points

  scatter = new scat(Zp,Z,A,flag,title0); 
  if (flag == 0) return;
  //___________________________________________________________________
  //read elastic scattering data


  int iData = 0;
  Ndata = 0;
  DegreesOfFreedom = 0;
  xsecMax = 0.;
  xsecMin = 1.E32;

  //open data file
  filename = directory + title + ".data";
  cout << filename << endl;
  ifstream fileData (filename.c_str(),ios::in);
  // if one cannot open file quit
  if (fileData.fail()) 
    {
      cout << "couldn't open data file" << fileData << endl;
    }

  else
    {
      for(;;)
	{
	  int n;
	  string xa;
	  double Elab;
	  int nfit;
	  fileData >> Elab >> nfit;
	  // if nfit ==1 data is included in fit, =0 plotted but not fitted
	  //=-1 data is ignored
          //=2  data is ratio to ruth and in fit
	  if (Elab < 0.) break;


	  getline(fileData,line);
	  getline(fileData,line);
	  cout << line << endl;

	  if (nfit >= 0)
	    {
	      data[iData].energyLab = Elab;
              data[iData].name = line;
	      data[iData].fit = nfit;
	      double Ecm = energyLab2Cm(Elab);
	      data[iData].energyCM = Ecm;
              // if we need to calculate rutherford, load in the Ecm to 
	      // scatter
	      if (nfit == 2) scatter->loadPotential(NULL,NULL,NULL,
						    NULL,NULL,1.,Ecm,Elab); 
	    }

	  if (fileData.eof()) break;
	  if (fileData.bad())
	    {
	      cout << "trouble reading file" << endl;
	      break;
	    }

	  if (nfit >= 0) Ndata++;

	  //make room for cross xsection data
	  fileData >> n >> xa;
          //cout << "n = " << n << " xa = " << xa << endl;
	  if (xa != "X") cout << "X  problem for Elab = " << Elab << endl;

	  if (nfit >= 0) data[iData].nX = n;
	  if (n > 0)
	    {  

	      if (nfit >= 0)
		{
		  data[iData].Xtheta = new double[n];
		  data[iData].xsec = new double[n];
		  data[iData].Xsigma = new double[n];
		}

	      double theta, xsec, sigma;
	      for (int i=0;i<n;i++)
		{
		  fileData >> theta >> xsec >> sigma;
     
		  if (nfit >= 0)
		    {

		      data[iData].Xtheta[i] = theta;
		      data[iData].xsec[i] = xsec;
		      data[iData].Xsigma[i] = sigma;
                      


                      // input data is given as ratio to rutherford
                      if( nfit == 2) 
			{
                         double ruth = scatter->Rutherford(theta*pi/180.); 
                         data[iData].xsec[i]*= ruth;
			
			 if (sigma > 0) data[iData].Xsigma[i] *= ruth;
			}		    


		      if (data[iData].Xsigma[i] < 0.) 
			data[iData].Xsigma[i] *= -data[iData].xsec[i]/100.;

		      if (data[iData].Xsigma[i] < data[iData].xsec[i]*.1)
			data[iData].Xsigma[i] = data[iData].xsec[i]*.1;

                      if ( Zp == 0 && A == 92 && Elab > 10 ) 
                      data[iData].Xsigma[i] /= 5.;
                      if ( Zp == 1 && Elab > 14. && Elab < 100 && 
			   (A==48 || A==44))
                      data[iData].Xsigma[i] /= 5.;

                      if ( Elab > 14. && Elab < 50 && Z == 28)
                      data[iData].Xsigma[i] /= 5.;

                      if ( Zp == 1 && Elab > 100. && A==9)
                      data[iData].Xsigma[i] /= 5.;

                      if ( Zp == 1 && Elab > 40. && Elab < 100. && A==208)
                      data[iData].Xsigma[i] /= 5.;


                      //if ( Zp == 1 && Elab > 40. && Elab < 100. && A==92)
                      //data[iData].Xsigma[i] /= 5.;

                      if ( Zp == 0 && Elab > 20. && Elab < 40. && A==208)
                      data[iData].Xsigma[i] /= 5.;

                      if ( Zp == 0 && Elab > 20. && Elab < 40. && A==9)
                      data[iData].Xsigma[i] /= 5.;

		      if (nfit) DegreesOfFreedom++;

		      xsecMax = max(data[iData].xsec[i],xsecMax);
		      xsecMin = min(data[iData].xsec[i],xsecMin);
		    }
		}
	    }
	  //analysing power data
	  fileData >> n >> xa;

	  if (xa != "A") cout << "A problem for Elab = " << Elab << endl;
	  if (nfit >= 0) data[iData].nA = n;
	  if (n > 0)
	    {
	      if (nfit >= 0)
		{
		  data[iData].Atheta = new double[n];
		  data[iData].anal = new double[n];
		  data[iData].Asigma = new double[n];
		}
	      double theta, pol, sigma;
	      for (int i=0;i<n;i++)
		{
		  fileData >> theta >> pol >> sigma;

		  if (nfit >=0)
		    {
		      data[iData].Atheta[i] = theta;
		      data[iData].anal[i] = pol;
		      if (sigma < 0.05) sigma = 0.05; //****** ******
		      data[iData].Asigma[i] = sigma;
                      if (data[iData].energyLab > 50 && data[iData].energyLab
			  < 100) data[iData].Asigma[i] /= 4.;
		      DegreesOfFreedom++;
		    }
		}
	    }

 
	  if (nfit >=0 )
	    {
	      iData++;
	      if (Ndata > NdataMax) 
		{
		  cout << "increase NdataMax" << endl;
                  abort();
		  break;
		}
	    }
	} 

      cout << Ndata << " blocks of data read from file " << filename << endl;
      cout << " degrees of freedom= " << DegreesOfFreedom << endl;
      cout << " Max xsection = " << xsecMax << endl;
      cout << " Min xsection = " << xsecMin << endl;

      fileData.close();
    }
  //--------------------- read in levels
  //read in levels
  filename = directory + title + ".lev";
  cout << filename << endl;
  ifstream fileLevel (filename.c_str(),ios::in);
  // if one cannot open file quit
  if (fileLevel.fail()) 
    {
      cout << "couldn't open data file" << endl;
      NFitLevels = 0;
      Nlevel = 0;
     LevelDegreesOfFreedom = 0;
    }
  else
    {
      getline(fileLevel,line);
     cout << line << endl;
      getline(fileLevel,line);
     //cout << line << endl;
     int Ilevel=0;

     LevelDegreesOfFreedom = 0;

     NFitLevels = 0;
     for (;;)
       {
           fileLevel >> Level[Ilevel].Energy >> Level[Ilevel].SigmaEnergy >>
           Level[Ilevel].N >>
	   Level[Ilevel].j >> Level[Ilevel].l >> Level[Ilevel].color >>
   	   Level[Ilevel].Efit >>
           Level[Ilevel].Rrms >> Level[Ilevel].SigmaRrms >>
	   Level[Ilevel].Rfit >>
           Level[Ilevel].Delta >> Level[Ilevel].SigmaDelta >> 
	   Level[Ilevel].Dfit >>
           Level[Ilevel].SpectFactor >> Level[Ilevel].SigmaSpect >>
           Level[Ilevel].Sfit;

	   //if (A >= 40.) Level[Ilevel].SigmaSpect /= 2.;
        

         if (fileLevel.eof()) break;


         if (Level[Ilevel].Efit) 
	   {
            LevelDegreesOfFreedom++;
	    NFitLevels++;
	   }
         if (Level[Ilevel].Rfit) LevelDegreesOfFreedom++;
         if (Level[Ilevel].Dfit) LevelDegreesOfFreedom++;
         if (Level[Ilevel].Sfit) LevelDegreesOfFreedom++;


         if (Ilevel == NlevelMax-1) 
	   {
	     cout << "increase NlevelMax" << endl;
             Ilevel++;
	     break;
	   } 
         Ilevel++;
       }
     Nlevel = Ilevel;
     cout << Nlevel << " levels read in " << endl;
     fileLevel.close();
    }



  //prepare for spectral functions
  Nsf = 440;
  for (int i=0;i<NlevelMax;i++) LevelTh[i].SpectFunction = new double [Nsf];
  Elow = Efermi - 22.;
  Ehigh = Efermi;
  Esf = new double[Nsf];
  for (int i=0;i<Nsf;i++) Esf[i] = -(double)i/4.+30;



  XsecDegreesOfFreedom = 0;

  //readin reaction xsection data --------
 filename = directory+title + ".xsec";
  cout << filename << endl;
  file.open(filename.c_str(),ios::in);
  // if one cannot open file quit
  if (file.fail()) 
    {
      cout << "couldn't open reaction *.xsec file" << endl;
      cout << " no reaction xsec in fit" << endl;
      NXdata = 0;
    }
  else 
    {
      getline(file,line);
      getline(file,line);
     file >> NXdata;
     XsecDegreesOfFreedom += NXdata;
     Xdata = new xdata [NXdata];

     for (int i=0;i<NXdata;i++)
       {
         file >> Xdata[i].energyLab >> Xdata[i].xsec >> Xdata[i].sigma;
         double Elab = Xdata[i].energyLab;
         double Ecm = energyLab2Cm(Elab);
         Xdata[i].energyCM = Ecm;
       }
    }
  file.close();
  file.clear();

  cout << "xsec points = " << NXdata << endl;
  if (Zp == 0.) //for neutrons total cross section data
    {
      filename = directory+title + ".txsec";
      cout << filename << endl;
      file.open(filename.c_str(),ios::in);
      // if one cannot open file quit
      if (file.fail()) 
	{
	  cout << "couldn't open reaction *.txsec file" << endl;
	  cout << " no total xsec in fit" << endl;
	  NTotXdata = 0;
	}
      else 
	{
          getline(file,line);
          getline(file,line);

	  file >> NTotXdata;
	  XsecDegreesOfFreedom += NTotXdata;
	  TotXdata = new xdata [NTotXdata];

	  for (int i=0;i<NTotXdata;i++)
	    {
	      file >> TotXdata[i].energyLab >> 
              TotXdata[i].xsec >> TotXdata[i].sigma;

	      //make sure the error bars are not too small
              if (TotXdata[i].sigma < TotXdata[i].xsec*.02)
		TotXdata[i].sigma = 0.02*TotXdata[i].xsec;


             
	      double Elab = TotXdata[i].energyLab;
	      double Ecm = energyLab2Cm(Elab);
              //if (A == 208 && Elab > 25 ) TotXdata[i].sigma /= 5;
	      TotXdata[i].energyCM = Ecm;
              if (TotXdata[i].sigma < 0.) TotXdata[i].sigma *= 
					    -TotXdata[i].xsec/100.;
	    }
	}
    }
  else NTotXdata = 0;
  cout << "tot xsec points" << NTotXdata << endl;

  file.close();
  file.clear();
  //readin integrated moments
 filename = directory+title + ".Vint";
  cout << filename << endl;
  file.open(filename.c_str(),ios::in);
  // if one cannot open file quit
  if (file.fail()) 
    {
      cout << "couldn't open integtaed moment *.Vint file" << endl;
      cout << " no moments in fit" << endl;
      Nmoment = 0;
    }
  else 
    {
       getline(file,line);
      int flag;
      file >> Nmoment >> flag;
      cout << Nmoment << " " << flag << endl;
      if (flag == 0) Nmoment = 0;

      if (Nmoment > 0)
	{
        string ref;
	for (int i=0;i<Nmoment;i++)
	  {
	    file >> Moment[i].EnergyLab >> Moment[i].Jreal >>
	      Moment[i].Jimag >> Moment[i].RMSreal >> Moment[i].RMSimag >>
	      Moment[i].Jso >> Moment[i].RMSso;
	    getline(file,ref);
            Moment[i].EnergyCM = energyLab2Cm(Moment[i].EnergyLab);

	    //cout << Moment[i].EnergyLab << " " << Moment[i].EnergyCM << endl;
	  }
	}
    }
  file.close();
  file.clear();
}
//*********************************************************************
  /**
   * constructor to use when fitting single-energy data
   \param title0 is the name of input file without extension
   \param jdata th energy data in title0.data is fit
   */
reaction::reaction(string *title0, int jdata, bool btxsec,bool banal)
{
  bprint = 0;
  for (int i=0;i<NlevelMax;i++) LevelTh[i].SpectFunction = NULL;
  Esf = NULL;
  NXdata = 0;
  NTotXdata = 0;

  DOM = 0;
  flag = 1;
  //directory = "/home/bob/DOMA/";
  directory = "";


  string line;
  title = *title0;
  string filename(title + ".inp");
  if (flag) cout << filename << endl;
  ifstream file (filename.c_str(),ios::in);
  // if one cannot open file quit
  if (file.fail()) 
    {
      cout << "couldn't open data file " << filename << " " << title << endl;
      abort();
    }

  file >> Zp >> Z >> A >> Efermi;
  file >> Nfermi;
  if (Nfermi < 1 || Nfermi > 2) cout << "Nfermi not possible" << endl;
  file >> ValenceHole.N >> ValenceHole.j >> ValenceHole.l;
  if (Nfermi == 2) file >> 
           ValenceParticle.N >> ValenceParticle.j >> ValenceParticle.l;
  asymmetry = (A-2.*Z)/A;
  if (Zp == 1.) sign = 1;
  else if (Zp == 0.) sign = -1.;
  else cout << " unknown Zp " << endl;
  if (file.bad()) cout << "problem with input file" << endl;
  if (file.eof()) cout << "eof with input file" << endl;

  file.close();
  file.clear();


  //asymmetry
  AsyVolume = 1;  // asymmetry for volume
  //alphaVolume = 1.65;
  //EaVolume = 140.;
  //EaVolume = 60.;
  

  //effective mass points

  scatter = new scat(Zp,Z,A,flag,title0); 
  //___________________________________________________________________
  //read elastic scattering data



  Ndata = 0;
  DegreesOfFreedom = 0;
  xsecMax = 0.;
  xsecMin = 1.E32;

  //open data file
  filename = directory + title + ".data";
  cout << filename << endl;
  ifstream fileData (filename.c_str(),ios::in);
  // if one cannot open file quit
  if (fileData.fail()) 
    {
      cout << "couldn't open data file" << fileData << endl;
    }

  else
    {
      int counter = 0;
      for(;;)
	{
	  int n;
	  string xa;
	  double Elab;
	  int nfit;
	  fileData >> Elab >> nfit;
	  // if nfit ==1 data is included in fit, =0 plotted but not fitted
	  //=-1 data is ignored
          //=2  data is ratio to ruth and in fit
          if (nfit == 0) nfit = 1;
	  if (Elab < 0.) 
	    {
	      data[0].nX = -1;
              break;
	    }


	  getline(fileData,line);
	  getline(fileData,line);
	  cout << line << endl;

	  if (jdata ==  counter)
	    {
	      data[0].energyLab = Elab;
              data[0].name = line;
	      data[0].fit = nfit;
	      double Ecm = energyLab2Cm(Elab);
	      data[0].energyCM = Ecm;
              // if we need to calculate rutherford, load in the Ecm to 
	      // scatter
	      if (nfit == 2) scatter->loadPotential(NULL,NULL,NULL,
						    NULL,NULL,1.,Ecm,Elab); 
	    }

	  if (fileData.eof()) break;
	  if (fileData.bad())
	    {
	      cout << "trouble reading file" << endl;
	      break;
	    }

	  //make room for cross xsection data
	  fileData >> n >> xa;
          //cout << "n = " << n << " xa = " << xa << endl;
	  if (xa != "X") cout << "X  problem for Elab = " << Elab << endl;

	  if (jdata == counter ) data[0].nX = n;
          else data[0].nX = 0;
	  if (n > 0)
	    {  

	      if (jdata == counter)
		{
		  data[0].Xtheta = new double[n];
		  data[0].xsec = new double[n];
		  data[0].Xsigma = new double[n];
		}

	      double theta, xsec, sigma;
	      for (int i=0;i<n;i++)
		{
		  fileData >> theta >> xsec >> sigma;
     
		  if (jdata == counter)
		    {
		      data[0].Xtheta[i] = theta;
		      data[0].xsec[i] = xsec;
		      data[0].Xsigma[i] = sigma;

                      // input data is given as ratio to rutherford
                      if( nfit == 2) 
			{
                         double ruth = scatter->Rutherford(theta*pi/180.); 
                         data[0].xsec[i]*= ruth;
			
			 if (sigma > 0) data[0].Xsigma[i] *= ruth;
			}		    


		      if (data[0].Xsigma[i] < 0.) 
			data[0].Xsigma[i] *= -data[0].xsec[i]/100.;

		      if (data[0].Xsigma[i] < data[0].xsec[i]*.1)
			data[0].Xsigma[i] = data[0].xsec[i]*.1;

		      if (nfit) DegreesOfFreedom++;

		      xsecMax = max(data[0].xsec[i],xsecMax);
		      xsecMin = min(data[0].xsec[i],xsecMin);
		    }
		}
	    }
	  //analysing power data
	  fileData >> n >> xa;

	  if (xa != "A") cout << "A problem for Elab = " << Elab << endl;
	  if (jdata == counter)
	    {
              if (banal)data[0].nA = n;
              else data[0].nA = 0;
            }
	  if (n > 0)
	    {
	      if (jdata == counter && banal)
		{
		  data[0].Atheta = new double[n];
		  data[0].anal = new double[n];
		  data[0].Asigma = new double[n];
		}
	      double theta, pol, sigma;
	      for (int i=0;i<n;i++)
		{
		  fileData >> theta >> pol >> sigma;

		  if (jdata  == counter && banal)
		    {
		      data[0].Atheta[i] = theta;
		      data[0].anal[i] = pol;
		      if (sigma < 0.05) sigma = 0.05; //****** ******
		      data[0].Asigma[i] = sigma;
		      DegreesOfFreedom++;
                      
                      
		    }
		}
	    }
          cout << jdata << " " << counter << endl;
          if (jdata == counter) break;
	  counter++; 
	} 

      cout << " degrees of freedom= " << DegreesOfFreedom << endl;
      cout << " Max xsection = " << xsecMax << endl;
      cout << " Min xsection = " << xsecMin << endl;

      fileData.close();
    }
  Ndata = 1; //only one data set to fit
  //************************************
  if (btxsec)
    {
      if (Zp == 0.) //for neutrons total cross section data
    {
      filename = directory+title + ".txsec";
      cout << filename << endl;
      file.open(filename.c_str(),ios::in);
      // if one cannot open file quit
      if (file.fail()) 
	{
	  cout << "couldn't open reaction *.txsec file" << endl;
	  cout << " no total xsec in fit" << endl;
	  NTotXdata = 0;
	}
      else 
	{
          getline(file,line);
          getline(file,line);

	  file >> NTotXdata;
	  TotXdata = new xdata [1];
          double ElabOld = 0.;
	  double xsecOld = 0.;
          double xsec,sigma,Elab;
	  for (int i=0;i<NTotXdata;i++)
	    {
	      file >> Elab >> xsec >> sigma;
              if (Elab < data[0].energyLab) 
		{
                 ElabOld = Elab;
                 xsecOld = xsec;
		}
	      else
		{
		  xsec = xsecOld + (xsec-xsecOld)/(Elab-ElabOld)*
		    (data[0].energyLab-ElabOld);
                  TotXdata[0].xsec = xsec; 
	          TotXdata[0].energyCM = data[0].energyCM;
	          TotXdata[0].energyLab = data[0].energyLab;
	          TotXdata[0].sigma = sigma;
                  NTotXdata = 1;
                  NXdata = 0;
                  XsecDegreesOfFreedom = 1;
                  cout << "tot xsec fit = " << TotXdata[0].xsec << endl;;
		  break;
		}
	    }
	}
    }
  cout << "tot xsec points" << NTotXdata << endl;
    }
  else
    {
     NXdata = 0;
     NTotXdata = 0;
     XsecDegreesOfFreedom  = 0;
    }

  //--------------------- read in levels
  NFitLevels = 0;
  Nlevel = 0;
  LevelDegreesOfFreedom = 0;



  Nmoment = 0;
}


//******************************************************************
  /**
   * Destructor
   */
reaction::~reaction()
{

  delete scatter;
  cout << "destroying reaction" << endl;
  if (flag == 0) return;
  for (int i=0;i<Ndata;i++)
    {
      if (data[i].nX > 0)
	{
         delete [] data[i].Xtheta;
         delete [] data[i].xsec;
         delete [] data[i].Xsigma;
	}
      if(data[i].nA > 0)
	{
          delete [] data[i].Atheta;
	  delete [] data[i].anal;
	  delete [] data[i].Asigma;
	}
    }
  for (int i=0;i<NlevelMax;i++) 
    {
      if (LevelTh[i].SpectFunction == NULL) continue;
     delete  [] LevelTh[i].SpectFunction;
    }
  delete [] Esf;
  if (NXdata > 0) delete [] Xdata;
  if (NTotXdata > 0) delete [] TotXdata;
  cout << " reaction distroyed" << endl;
}
//******************************************************************
  /**
   *reinitializes all the energy-dependent potentials in the class
   *scatter
   \param Ecm is the energy in the center-of-mass frame [MeV]
   \param Elab is the energy in the laboratory frame [MeV]
   */
void reaction::InitializeForEcm(double Ecm,double ELab)
{
  if (DOM)
    {
      //find the energy dependent volume Imaginary potential and its 
      //correction to the real potential
      Volume.SetEnergy(Ecm);
      SurfaceLow.SetEnergy(Ecm);
      SurfaceHigh.SetEnergy(Ecm);
      HartreeFock.SetEnergy(Ecm);
      SpinOrbit.SetEnergy(Ecm);
    }

      // initialization to OM calculation
      scatter->loadPotential(&HartreeFock,&Volume,&SurfaceLow,&SurfaceHigh,
       &SpinOrbit,Rc,Ecm,ELab);



}


//*****************************************************************
  /**
   * loads in in new parameters and adjusts the depth of the Hartree
   * Fock potential to ge the Fermi energy correct
   */
bool reaction::prepare()
{
  load();
  return AdjustFermi();
}

//******************************************************************
  /**
   * returns chi squared per degree of freedom of the fit
   */
double reaction::ChiSquared()
{
  double sum = 0.;
  double sumXsec = 0.;
 

  
  //loop over number of energies in fit
  for(int i=0;i<Ndata;i++)
    {
      if (data[i].fit == 0) continue;
      double Ecm =data[i].energyCM;
      double Elab =data[i].energyLab;

 
      InitializeForEcm(Ecm,Elab);
       // integrate wavefunction and find phaseshifts
       if (scatter->integrateWave()==0) return 1.e6;
       //add prepare compound elastic
       if (Ecm < 15.) scatter->statistical(scatter->konst,Ecm);

       // loop over angles
       for (int j=0;j<data[i].nX;j++)
	 {
	   double angle = data[i].Xtheta[j]*pi/180.;
           double xsecTheory = scatter->DifferentialXsection(angle);
           // add in compound elastic
           if (Ecm < 15.) xsecTheory += scatter->DifferentialXsectionCE(angle);
           //add to chi squared
           sum += pow(xsecTheory-data[i].xsec[j],2)/
	     pow(data[i].Xsigma[j],2);


	 }
       for (int j=0;j<data[i].nA;j++)
	 {
	   double angle = data[i].Atheta[j]*pi/180.;
           double xsecTheory = scatter->DifferentialXsection(angle);
           //add to chi squared

           double A = scatter->AnalyzePower;
	   if (Ecm < 15.) 
	     {
               double shape = scatter->DifferentialXsectionCE(angle);
               A *= xsecTheory/(xsecTheory+shape);
	     }

           sum += pow(scatter->AnalyzePower-data[i].anal[j],2)/
	     pow(data[i].Asigma[j],2);
	 }

    }
  
  //sum *= 10.; //rjc used for fitting n+48Ca
  //-------------------------------------------------------------------


  // reaction xsections
  
  for(int i=0;i<NXdata;i++)
    {
      double Ecm = Xdata[i].energyCM;
      double Elab = Xdata[i].energyLab;
 
      InitializeForEcm(Ecm,Elab);

       // integrate wavefunction and find phaseshifts
       scatter->integrateWave();

       //fit absorption cross section for protons
       double xsec = scatter->AbsorptionXsection(); //protons

       //remove compound elastic part
       if (Ecm < 15.) xsec -= - scatter->statistical(scatter->konst,Ecm);

       sumXsec += pow((xsec-Xdata[i].xsec)/Xdata[i].sigma,2); //RJC
    }
  
  //total cross section
  //sumXsec = 0.;



  if (NTotXdata > 0 && Zp == 0.)
    {
      for(int i=0;i<NTotXdata;i++)
	{

	  double Ecm = TotXdata[i].energyCM;
	  double Elab = TotXdata[i].energyLab;
 
	  InitializeForEcm(Ecm,Elab);

	  // integrate wavefunction and find phaseshifts
	  scatter->integrateWave();


	  //fit total cross section for neutrons
	  double xsec = scatter->TotXsection(); //neutrons
	  sumXsec += pow((xsec-TotXdata[i].xsec)/TotXdata[i].sigma,2); //RJC

	}
    }
  



  //-------------------------------------------------
  



  // include levels in chi squared
  double sumLevel=0.;
  double sumR = 0.;
  double sumD = 0.;
  double sumS = 0.;
  double sumE = 0.;

  int FoundLevels = 0;
  for (int i=0;i<Nlevel;i++)

    {
      if (Level[i].Efit == 0) continue;
      double Elower;
      if (Z > 20.) Elower = -30.;
      else Elower = floor(Efermi - 25.);
      if (Level[i].N == 0 && Level[i].l == 0) Elower = -90.;
      if (Level[i].N == 0 && Level[i].l == 1) Elower = -55.;

      // rjc to get wavefunction information, Ifine = 1
      int const Ifine = 1;
      double Eupper = (float)scatter->LogDerMax;
      if (scatter->proton && Eupper > 10) Eupper = 10.;

      FoundLevels += bound(Elower,Eupper,Level[i].j,Level[i].l,Ifine);
      // check to see if the right N of level


      int nodes = scatter->Nodes();

      /*
      if (nodes > Level[i].N) cout << "warning level order too high " << nodes 
				   << " " << Level[i].N << " " << i << 
				" Energy= " << Elower << endl;
      */
      if (nodes< Level[i].N)
	{
          // need to search for highly lying level
          //cout << "searching for higher order level " << Level[i].N
	  //     << " " << Level[i].l << " " << Level[i].j << endl;
	  //cout << "energy was = " << Elower << endl;
	  FoundLevels--;
          Elower += .5;
          FoundLevels += bound(Elower,Eupper,Level[i].j,Level[i].l,Ifine);
          //cout << "energy is now " << Elower << endl;
	}
      if(Level[i].Rfit+ Level[i].Dfit + Level[i].Sfit > 0)
         LevelProperties(Elower,Level[i].j,Level[i].l);

      if (Level[i].Efit)
	sumE += pow(Elower-Level[i].Energy,2)/pow(Level[i].SigmaEnergy,2);


      if (Level[i].Rfit)
        sumR += pow(scatter->Rrms-Level[i].Rrms,2)
        /pow(Level[i].SigmaRrms,2)*4.;

      if (Level[i].Dfit)
        sumD += pow(scatter->Width-Level[i].Delta,2)
        /pow(Level[i].SigmaDelta,2)*40.;


      if (Level[i].Sfit)
	{
      sumS += pow(scatter->SpectFactor - Level[i].SpectFactor,2)
      /pow(Level[i].SigmaSpect,2)*20.;
	}
    }




  if (FoundLevels != NFitLevels) 
    {
      cout << " only " << FoundLevels << " levels found of " << NFitLevels <<
	" fit in " << title <<  endl;
      //Fout << " only " << FoundLevels << " found" << endl;
      sumLevel = 1.e9;
      sumE += 50.;
    }


  //--------------------------------------------------------
  //integrated moments
  double sumMoments = 0.;
  /*
  for (int i=0;i<Nmoment;i++)
    {
      InitializeForEcm(Moment[i].EnergyCM,Moment[i].EnergyLab);
      scatter->VIntegrals();
      if (Moment[i].Jreal != 0.) 
          sumMoments += pow((Moment[i].Jreal+scatter->JReal)/20.,2);
      if (Moment[i].Jimag != 0.) 
         sumMoments += pow((Moment[i].Jimag+scatter->JImag)/20.,2);
      //if (Moment[i].Jso != 0.) 
      //   sumMoments += pow((Moment[i].Jso+scatter->JSO)/20.,2);
      if (Moment[i].RMSreal != 0.) 
         sumMoments += pow((Moment[i].RMSreal-scatter->RrmsReal)/.02,2);
      if (Moment[i].RMSimag != 0.) 
         sumMoments += pow((Moment[i].RMSimag-scatter->RrmsImag)/.02,2);
      //if (Moment[i].RMSso != 0.) 
      //   sumMoments += pow((Moment[i].RMSso-scatter->RrmsSO)/.3,2);

    }
  
  */

  double total = 0.;
  double chiPoint[4] = {0.,0.,0.,0.};
  if (DegreesOfFreedom >0)
    {
     chiPoint[0] = sum/DegreesOfFreedom;
     total += chiPoint[0];
    }
  if (LevelDegreesOfFreedom > 0) 
    {
     
     sumLevel = sumE + sumR + sumD + sumS;
     chiPoint[1] = sumLevel/LevelDegreesOfFreedom;
     total += chiPoint[1];
     sumE /= LevelDegreesOfFreedom;
     sumR /= LevelDegreesOfFreedom;
     sumD /= LevelDegreesOfFreedom;
     sumS /= LevelDegreesOfFreedom;
    }
  if (XsecDegreesOfFreedom > 0) 
    {
     chiPoint[2] = sumXsec/XsecDegreesOfFreedom*10.;
     total += chiPoint[2];
    }
  if (Nmoment > 0)
    {
      chiPoint[3] = sumMoments/(6.*(double)Nmoment);
      total += chiPoint[3];
    }

  if (bprint)
    {
  cout << chiPoint[0] << " " << chiPoint[1] << " " << chiPoint[2] << 
    " " << chiPoint[3] << " E= " << sumE << " rms= " << sumR 
    << " delta=" << sumD << " S= " << sumS << endl;
    }


  return total;
}
//***************************************************************
  /** determines the bound state and quasibound states energies.
    * if one wants to get the wavefunction, you need to get the energy
    *  much more accurately than this, use Ifine=1
    * search is made in the energy range Elower to Eupper.
    * found level is returned as Elower.
    *function returns 1 is level is found, else 0.
    \param Elower is the lower bound for energy range of search  [MeV]
    \param Eupper is the upper bound for energy range of search  [MeV]
    \param j total angular momentum of the single-particle state
    \param l is the orbital angular momentum of the single-particle state
    \param Ifine fine search of energy - need if you want wavefunctions
  */
int reaction::bound(double& Elower,double Eupper, double j, int l, int Ifine)
{

  double const deltaE = 1.;
  double Ecm  = Elower;

  double Wave;
  double WaveOld = 0.;
  double WaveAbove = 0.;
  double EcmAbove = 0.;
  int tries = 0;
  int last = 0;
  valarray<double> WaveFunct;
  for (;;)
    {
      // find all energy dependent OM terms and initialize OM 
      InitializeForEcm(Ecm,0.);


       // integrate wavefunction and find phaseshifts
       WaveFunct = scatter->IntegrateBound(j,l);
       //Wave = WaveFunct[0];
       Wave = scatter->LogDerDifference(j,l);
       //if (l == 4) cout << Ecm << " " << Wave << endl;

       if (last) 
	 {
           if (Ifine)
	     {

	       if (Wave*WaveOld < 0.) 
		 {
		   Elower = Ecm - deltaE/2.;
                   Eupper = Ecm;
		 }
               else 		 {
                   Elower = Ecm;
                   Eupper = EcmAbove;
		 }
               return boundFine(Elower,Eupper,j,l);

	     }
	   else  // last iteration, use Ridder method for an interpolation
	 //see Numerical Recipes in Fortran77 Page351
	     {
               double sign = 1.;
               if (WaveOld < WaveAbove) sign = -1.;
               Elower = Ecm + deltaE/2.*sign*Wave/
	        sqrt(pow(Wave,2)-WaveOld*WaveAbove);
               return 1;
	     }
	 }

       if (tries > 0)
	 {
	   if (Wave*WaveOld < 0.)
	     {
               EcmAbove = Ecm;
               Ecm -= deltaE/2.;
               WaveAbove = Wave;
	       last = 1; // flag the last iteration
	       continue;
	     }
	 }

       tries++;
       WaveOld = Wave;
       Ecm += deltaE;
       if (Ecm > Eupper) break;
 
    }
  return 0;
}
//************************************************************************
  /**
   * finds the eigenstates accuracy enough so that the Wavefunction
   *can be determined to give physical properties of the state.
   *Ridder method for an interpolation,
   *see Numerical Recipes in Fortran77 Page351
    \param Elower is the lower bound for energy range of search  [MeV]
    \param Eupper is the upper bound for energy range of search  [MeV]
    \param j total angular momentum of the single-particle state
    \param l is the orbital angular momentum of the single-particle state
   */
int reaction::boundFine(double& Elower,double Eupper, double j, int l)
{
  int counter = 0;
  if (Elower >= Eupper) 
  cout << "Hey buddies, what gives? Elower >= Eupper" << endl;

  double Ecm;

  // start at lower limit, find difference in slope
  // of integrated wavefunction and matched whittaker function
  Ecm = Elower;
  InitializeForEcm(Ecm,0.);
  double ylower = scatter->LogDerDifference(j,l);

  //------------------------------------------------------
  // next the upper limit
  Ecm = Eupper;
  InitializeForEcm(Ecm,0.);
  double yupper = scatter->LogDerDifference(j,l);




  if (yupper*ylower > 0.)
    { 
    cout << "hey what gives, zero is not enclosed by limits" << endl;
    cout << ylower << " " << yupper << endl;
    cout << Elower << " " << Eupper << endl;
    }

  for (;;)
    {

      //cout << Elower << " " << ylower << " " << 
      // Eupper << " " << yupper << endl;

      double Emiddle  = (Elower+Eupper)/2.;
      double Ecm  = Emiddle;
 
      //return if precision is too small
      if (Ecm == Elower) return 1;
      if (Ecm == Eupper) return 1;

      InitializeForEcm(Ecm,0.);
      double ymiddle = scatter->LogDerDifference(j,l);
      //cout << "mid = " << Ecm << " " << ymiddle << endl;

       double sign = -1.;
       if (ylower > yupper) sign = 1.;
       
       Ecm = Emiddle + (Emiddle-Elower)*sign*ymiddle/
	 sqrt(pow(ymiddle,2)-ylower*yupper);

       InitializeForEcm(Ecm,0.);
       double y = scatter->LogDerDifference(j,l); 
       //cout << " new " << Ecm << " " << y << endl;

       if (abs(y) < .00001) 
	 {
	   Elower = Ecm;
           return 1;
	 }

       // find new bounding energies
       if (Ecm < Emiddle)
	 {
           if (y*ylower < 0.)
	     {
	       Eupper = Ecm;
               yupper = y;
	     }
	   else 
	     {
               Elower = Ecm;
               Eupper = Emiddle;
               ylower = y;
               yupper = ymiddle;
	     }
	 }
       else 
	 {
           if (y*yupper < 0.)
	     {
               Elower = Ecm;
               ylower = y;
	     }
	   else
	     {
               Elower = Emiddle;
               Eupper = Ecm;
               ylower = ymiddle;
               yupper = y;
	     }
	 }
       if (counter > 11) return 0;
       counter++;

    }
  return 0;
}

//*********************************************************************
/**
 *bound and quasi bound states are found when this function returns zero
 *Originally, I matched the integrated wave function to the Whittaker function
 *and returned the difference in slope. However at the matching radius I
 *generally use, it is just as accurate to find the solutions where integrated 
 *wave function is zero at the matching radius. this is quicker as one does 
 *not need to calculate the whittaker function.
 */

double reaction::SlopeDifference(double Ecm, double j, int l)
{

  // find all energy dependent OM terms and initialize OM 
  InitializeForEcm(Ecm,0.);
  // integrate wavefunction
  double derivative; 
  scatter->GetWaveFunctionArray(j,l,derivative);
  // find inertior log derivative at matching point
  double interiorLogDer = derivative/scatter->WaveArray[scatter->mWave-1];
  //find exterior log derivative at matching point

  double exteriorLogDer = scatter->exteriorLogDer(l);
  return interiorLogDer - exteriorLogDer; 


}
  //***************************************************************
    /**
     * finds the wavefunction and normalizes it
     * This function must be called first before many of the other
     * function associated with bound-state properties are called
     \param Ecm energy of bound state in MeV
     \param j total angular momentum of bound state
     \param l orbital angular momentum of boubd state
     */
  void reaction::LevelProperties(double Ecm,double j, int l) 
  {
      double derivative;
      InitializeForEcm(Ecm,0.);
      scatter->GetWaveFunctionArray(j,l,derivative);
      scatter->exteriorWaveFunct(l); // add wavefunction beyound matching region
      scatter->normalizeWaveFunction(Ecm,Efermi);
  }

//*********************************************************************
  void reaction::spectralFunction(double E0, int N, double* E, double * out)
{
  for (int i=0;i<N;i++)
    {
      double Ecm = E[i];

      //find the energy dependent volume Imaginary potential
      Volume.SetEnergy(Ecm);
      SurfaceLow.SetEnergy(Ecm); 
      SurfaceHigh.SetEnergy(Ecm); 

      double W = scatter->spectralWidth()/scatter->AverageMassEff;

     out[i] =  -scatter->SpectFactor*W/pi/(pow(Ecm-E0,2)+pow(W,2));

    }
}
//***************************************************************
  //gives the strength in an energy interval
  double reaction::strengthE(double* spectF, double Eres, double Width, 
       double SpectFactor)
{
  double sumIn = 0.;
  double sumOut = 0.;
  for (int i=0;i<Nsf;i++)
    {
      if (Esf[i] > Elow && Esf[i] < Ehigh) sumIn+= spectF[i];
      else sumOut += spectF[i];
    }

  double delta = abs(Esf[1]- Esf[0]);
  sumIn *= delta;
  sumOut *= delta;
  // if the width of the state is really narrow, then the step size for the 
  // Spectral function is probably not large enough to integrate over the peak,
  // thus we integrate outside of the energy inteval and subtract from the 
  // total (spectral factor) 
  if (Width < 1. && Eres > Elow && Eres < Ehigh) sumIn = SpectFactor - sumOut;
  return sumIn;
}
//**********************************************************
  /**
   * writes to files the tranmission coefficients
   * when s and l are parallel file is filename_up.tl
   * else filename_down.tl
    \param filename gives names of files
   */
void reaction::PrintTransCoef(string filename)
{
  int const NstepsE = 200;
  int const NstepsL = 11;

  struct tt
  {
    double down[NstepsE];
    double up[NstepsE];

  };
  tt transm[NstepsL+1];
  double Earray[NstepsE];
  string name = filename+"_up.tl";
  ofstream FileUp(name.c_str());
  name = filename+"_down.tl";
  ofstream FileDown(name.c_str());
  FileUp << " " << NstepsL +1 << " " << NstepsE << endl;
  FileDown << " " << NstepsL + 1 << " " << NstepsE << endl;
  for(int i=0;i<NstepsE;i++)
    {
      //double Ecm =((double)i+.5)/10.;
      double Ecm =((double)i+.5)/2.5; //rjc
      Earray[i] = Ecm;
      double Elab = Ecm*(A+1.)/A;
      InitializeForEcm(Ecm,Elab);

       // integrate wavefunction and find phaseshifts
       scatter->integrateWave();
       FileUp << Ecm << " ";
       FileDown << Ecm << " " ;
       for (int jj=0;jj<=NstepsL;jj++)
           {
	     double ttup =  scatter->TransCoef(jj,(double)jj+.5);
	     double ttdown =  scatter->TransCoef(jj,(double)jj-.5);

             // sometime the transmission coeff can be less than zero
             // due to numerical error
             if (ttup > -.0001 && ttup < 0.) ttup = 0.;
             if (ttdown > -.0001 && ttdown < 0.) ttdown = 0.;

	     if (ttup < 0. || ttdown < 0.)
	       cout << "warning transmission coef negative" << endl;
	     FileUp <<ttup << " " ;
	     FileDown <<ttdown << " ";  
             transm[jj].up[i] = ttup;
	     transm[jj].down[i] = ttdown;
           }
       FileUp << endl;
       FileDown << endl;
    }
  FileUp.close();
  FileDown.close();


#ifdef root
  string nnn = "trans_"+title;
  TCanvas Ctl(nnn.c_str());

 // TH2S histTr("histTr",title.c_str(),10,0,20,10,0,1); rjc
  TH2S histTr("histTr",title.c_str(),10,0,80,10,0,1);
  histTr.SetStats(kFALSE);
  histTr.GetXaxis()->SetTitle("Ek_{cm} [MeV]");
  histTr.GetYaxis()->SetTitle("Transmission Coefficient");
  histTr.Draw();
  
  TGraph *graphUP[NstepsL];
  TGraph *graphDOWN[NstepsL];
  for (int jj=0;jj<NstepsL;jj++)
    {
      
      graphUP[jj] = new TGraph(NstepsE,Earray,transm[jj].up);
      graphUP[jj]->SetLineStyle(1);
      graphUP[jj]->SetLineColor(jj+1);
      graphUP[jj]->SetLineWidth(2);
      graphUP[jj]->Draw("C");
      
      if (jj >= 1)
	{

      graphDOWN[jj] = new TGraph(NstepsE,Earray,transm[jj].down);
      graphDOWN[jj]->SetLineStyle(2);
      graphDOWN[jj]->SetLineColor(jj+1);
      graphDOWN[jj]->SetLineWidth(2);
      graphDOWN[jj]->Draw("C");

	}
      
    }
 
  Ctl.Write();
 
#endif
}
//*********************************************************************
  /**
   *relativistic conversion from Lab kinetic energyto total CM kinetic energy
   \param Elab energ in laboratory frame [MeV]
  */
double reaction::energyLab2Cm(double Elab)
{
  if (Elab < 0.) return Elab*A/(1.+A);
 // center of mass velocity in units of c
 double vcm = sqrt(Elab*(Elab+2.*scat::m0))/(Elab+(1.+A)*scat::m0);
 //gamma factor for this velocity
 double gam = 1./sqrt(1.-pow(vcm,2));
 double Ecm = (gam-1.)*(1.+A)*scat::m0 +
               gam*Elab*(A-1.)*scat::m0/((A+1.)*scat::m0+Elab);
 return Ecm;
}
//**********************************************************************
  /**
   *relativistic conversion from Lab kinetic energy to total CM kinetic energy
   \param Ecm is energy in the center-of-mass frame
  */
double reaction::energyCm2Lab(double Ecm)
{
  if (Ecm < 0.) return Ecm/A*(1.+A);

  //find momentum of projecile, also momentum of target
  double pc = sqrt(Ecm)*sqrt((Ecm+2.*scat::m0)*(Ecm+2.*A*scat::m0)*
              (Ecm+2.*(A+1.)*scat::m0))/2./(Ecm+(A+1)*scat::m0);
  //velocity of target in units of c
  double vtarget = pc/sqrt(pow(pc,2)+pow(A*scat::m0,2));
  //gamma factor for this velocity
  double gam = 1./sqrt(1.-pow(vtarget,2));
  // tot energy of projectile (neutron or proton in com frame)
  double Eproj = sqrt(pow(scat::m0,2)+pow(pc,2));
  double Elab = (Eproj + vtarget*pc)*gam;
  //this energy contains rest mass , so remove it 
  Elab -= scat::m0;
  return Elab;
}

//**********************************************************************
  /**
   * load all the paramter values into the classes associated with the 
   * optical model potential
   */
void reaction::load()
{

  HartreeFock.load(RHF,aHF,VHF,alpha,beta,gamma,Efermi,alphaS,betaS,gammaS);
  SpinOrbit.load(Rso,aso,Vso,VspinoE,Efermi,AWso,BWso);

  //initialize the energy dependent volume Imaginary potential and its 
  //correction to the real potential
  Volume.load(Rvolume,deltaRvolume,expRvolume,avolume,Avolume,Bvolume,
      Epvolume,Efermi,mvolume,AsyVolume,alphaVolume,EaVolume);  


  //initialize the energy dependent surface Imaginary potential and its 
  //correction to the real potential

  SurfaceLow.load(Rsurface,asurface,Asurface,Bsurface,Csurface,
	       Dsurface,Wstart*fGap,Efermi);

  SurfaceHigh.load(Rvolume,avolume,aHigh,bHigh,cHigh,
		   Epvolume,mvolume,Efermi);

}

//**********************************************************************
  /**
   * load all the paramter values into the classes associated with the 
   * optical model potential
   */
void reaction::loadOM()
{

  HartreeFock.load(VHF,0.,RHF,aHF);
  SpinOrbit.load(Vso,AWso,Rso,aso);

  Volume.load(Avolume,Rvolume,avolume);

  SurfaceLow.load(Asurface,Rsurface,asurface);
  SurfaceHigh.load(0.,Rsurface,asurface);

}


//*********************************************************************
  /**
   * Calculate the Fermi Energy in MeV
   */
double reaction::FindFermi()
{
      double Eupper = (float)scatter->LogDerMax;
      double Estart = Efermi-15.;
      int Ifine = 0;
      double Elower;
      int count = 0;
      double Efermi_calculate;
      for (;;)
	{
         Elower = Estart;


         bound(Elower,Eupper,ValenceHole.j,ValenceHole.l,Ifine);
         if (Elower > Estart) break;
         Elower -= 5.;
         count++;
         if (count > 4)
	   {
            cout << "hole level not found" << endl;
	    Elower = 15.; //set the hole level high so the chisq will be large
            break;
	   }
	}

      double Ehole = Elower;
      if (Nfermi > 1)
	{
         bound(Elower,Eupper,ValenceParticle.j,ValenceParticle.l,Ifine);
         if (Elower == Ehole) cout << "particle level not found" << endl;
         double Eparticle = Elower;
         Efermi_calculate = (Ehole+Eparticle)/2.;
	//cout << Ehole << " " << Eparticle << " " << Efermi_calculate << endl;
	}
      else 
	{
         Efermi_calculate = Ehole;
         //cout << Ehole << endl;
	}
      return Efermi_calculate;
}
//***************************************************************
  /**
   * adjusted the depth of the HF potential so that the calculated
   *Fermi energy is the same as the experimental Fermi energy
   */
bool reaction::AdjustFermi()
{
  int counts=0;
  double Efermi_calculate_old = 0.;
  double VHF_old = 0.;
  for (;;)
    {


      double Efermi_calculate = FindFermi();
      if (fabs(Efermi-Efermi_calculate) < .04) return 1;


      if (counts < 1) 
	{
         VHF_old = VHF;
         VHF -= (Efermi-Efermi_calculate)*1.2;
	}
      else
	{
          if (Efermi_calculate == Efermi_calculate_old)
	    {
	      cout << "Efermi_calculate == Efermi_calculate_old =" << 
		Efermi_calculate << endl;
              cout << "Efermi=" << Efermi << endl;
              cout << "VHF=" << VHF << " VHF_old= " << VHF_old << endl;
              cout << "counts=" << counts << endl;
              cout << "zp=" << Zp << endl;
	      //abort();
	      return 0;
	    }          

	  double VHF_new = VHF_old + (VHF-VHF_old)/
          (Efermi_calculate-Efermi_calculate_old)
           *(Efermi-Efermi_calculate_old);
          VHF_old = VHF;
          VHF = VHF_new;
          if (counts > 4) VHF = (VHF+VHF_old)/2.;
	}

      Efermi_calculate_old = Efermi_calculate;
      HartreeFock.load(RHF,aHF,VHF,alpha,beta,gamma,Efermi,alphaS,betaS,gammaS);

      counts++;
      if (counts > 10) 
	{
	  cout << "problems getting correct Fermi energy in AdjustFermi" 
	       << "VHF= "<< VHF << " VHF_old= " << VHF_old  
	       << " Efermi_calculate = " << Efermi_calculate 
	       << " Efermi_calculate_old= " << Efermi_calculate_old
	       << " Efermi=" << Efermi
               << endl;
          //abort();
         return 0;
         
	}
    }
}
//***************************************************
  /**
   *write out in a file title.out, the parameter AHF required to give the 
   *experimental fermi energy
   */
void reaction::WriteAHF()
{
  string filename(title + ".out");
  ofstream file (filename.c_str());
  double W = SurfaceLow.getMaxW();
  file << VHF << " " << W << endl;
  file.close();
  file.clear();
  string filename2(title + ".HF");
  file.open (filename2.c_str());
  double Ecoul = 0.;
  if (Zp == 1)
    {
     double R = 1.24*pow(A,(1./3.));
     Ecoul = Z*1.73/R;
     }

  HartreeFock.SetEnergy(Ecoul);
  scatter->VIntegrals();
  //file << HartreeFock.V << " " << scatter->JReal << endl;
  file << scatter->JReal <<  " " << HartreeFock.V << endl;
  file.close();


}
//*****************************************************
  /**
   *This subroutine is for cases where the Fermi energy is not known 
   *experimentally. Aan iterative procedure finds a consistant value of 
   *the Fermi energy
   */
double reaction::IterateFermi()
{
  //start with guess of Efermi from Jeukenne et al, PRC43 (1991) 2211
  if (Zp == 0) Efermi = -12.52 +31.3*asymmetry;
  else Efermi = -11.18 -57.5*asymmetry + 1.73*Z/(1.39+1.04*pow(A,.33333));
  cout << "Efermi from Jeukenne" << endl;
  double Efermi_old;
  load();
  int counter = 0;
  for(;;)
    {
     Efermi_old = Efermi;
     Efermi = FindFermi();
     if (abs(Efermi-Efermi_old) < 0.1) break;
     if (counter > 10) Efermi = (Efermi+Efermi_old)/2.;
     cout << "Efermi=" << Efermi << endl;
     load();      

     counter++;
    }
  string filename(title + ".fermi");
  ofstream file (filename.c_str());
  file << Efermi << endl;
  file.close();
  return Efermi;
}
//***************************************************
void reaction::getTransferWaveFct()
{
  //adjust the HF depth to get the valence hole energy level correct

     double Eupper = (float)scatter->LogDerMax;
     double Estart = Efermi-15.;
     double VHFstart = VHF;
     int const Ifine = 1;


     double Elower;
     double ElowerOld;
     double VHFold;
     int count = 0;

     for (;;)
       {
	count ++;
        ElowerOld = Elower;
        Elower = Estart;
        bound(Elower,Eupper,ValenceHole.j,ValenceHole.l,Ifine);
        double delta = Elower - ValenceHole.energy;
        if (abs(delta) < .001) break;

        if (count > 1)
          {
	    double slope = (Elower-ElowerOld)/(VHF-VHFold);
            VHFold = VHF;
            VHF -= delta/slope;
	  }
        else 
	  {
            VHFold = VHF;
            VHF += delta;
	  }


        HartreeFock.load(RHF,aHF,VHF,alpha,beta,gamma,Efermi,alphaS,betaS,gammaS);

       }

     /*
     cout << "shift in VHF is " << VHF-VHFstart << " MeV" << endl;

     LevelProperties(Elower,ValenceHole.j,ValenceHole.l);

     for (int i = 0;i<scatter->nWave;i++)
       {
	 cout << scatter->rStart + (double)i*scatter->deltaR 
              << " " << scatter->WaveBar[i] << endl;
       }
     */
}
