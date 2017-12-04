#ifndef _reaction
#define _reaction

#include "scat.h"
#include <string>
#include <algorithm>
#include <iomanip>
#include <sstream>
#include "volume.h"
#include "surfaceTF.h"
#include "surfaceSTD.h"
#include "hartreeFock.h"

#include <fstream>
#define root
#ifdef root
  #include "TFile.h"
  #include "TCanvas.h"
  #include "TGraphErrors.h"
  #include "TH2S.h"
  #include "TPad.h"
  #include "TLine.h"
  #include "TLatex.h"
#endif

struct xdata
{
  double energyLab;
  double energyCM;
  double xsec;
  double sigma;
};


struct datas
{
  double energyLab;
  double energyCM;
  string name;
  int nX;
  int nA;
  int fit;
  double* Xtheta;
  double* xsec;
  double* Xsigma;
  double* Atheta;
  double* anal;
  double* Asigma;
};

struct lev
{
  int N;
  double j;
  int l;
  double energy;
};

struct levels
{
  double Energy; // level energy
  double SigmaEnergy;  // experimetal uncertainty
  int N;
  double j;
  int l;
  int color;
  int Efit; // fit level energy
  int Rfit; // fit level rms radius
  int Dfit; // fit level width
  int Sfit; // fit level spectroscopic factor
  double Rrms;     //rms radius
  double SigmaRrms;  //exp uncertainty
  double Delta; // level width (MeV)
  double SigmaDelta; // experimental uncertainty
  double SpectFactor; //spectroscopic Factor
  double SigmaSpect; // error in Spect Factor
};

//calculated level properties
/*
 *\brief stores information on calculated levels
 */
struct levelTh
{
  double N; //!< number of nodes in the radial wavefunction
  double j; //!< total angular momentum of level
  int l; //!< orbital angular momentum of level
  double energy; //!< level energy in MeV
  double Rrms; //!< rms radius of wavefunction
  double SpectFactor; //!< Spectroscopic factor
  double Occupation;  //!<occupation probability of level
  double Width;  //!< width of level in MeV
  int color;
  string name; //!< spectscopic name of level in latex
  double *SpectFunction;
  double SpectE;
  double ANC; //!< asymptotic normalization coefficient
};
int const NdataMax = 40;
int const NlevelMax = 30;


//integrated moment
struct moment
{
  double EnergyLab;
  double EnergyCM;
  double Jreal;
  double Jimag;
  double Jso;
  double RMSreal;
  double RMSimag;
  double RMSso;
};


  /**
   *\brief deal with a single reaction, ie, p+40Ca
   *
   * performs dispersive optical model calculations of elastic scattering, 
   * reaction and toal sections, bound state energies, RMS radii, widths, and 
   * spectroscopic factors for a single reaction. 
   */


class reaction
{
 public:
  bool DOM; //!< indicates that a DOM fit is taking place
  static double const pi;
  reaction () {};
  reaction (string *title0,bool flag0);
  reaction (string *title0,int jdata,bool btxsec, bool banal);
  ~reaction();
  int Ndata; //!< number of data sets
  double xsecMax;
  double xsecMin;
  int DegreesOfFreedom;
  int LevelDegreesOfFreedom;
  int XsecDegreesOfFreedom;
  bool flag;

  int Nfermi; //!< number of levels defining the Fermi energy
  lev ValenceHole;
  lev ValenceParticle;

  //elastic scattering data
  datas data[NdataMax];
  string title;
  string directory;

  //reaction xsection data
  xdata *Xdata;
  int NXdata;

  //total xsection data
  xdata *TotXdata;
  int NTotXdata;

  bool prepare();
  double ChiSquared(); // returns the chi squared.
  void PlotFit(); // plots data and fitted curve
  void PlotPotentialEcm(); //plots potential as function of Ecm energy 
  void plotSmatrix();
  void printSmatrix();
  void printSmatrix2();
  void PrintTransCoef(string);
  int bound(double&,double,double,int,int);
  int boundFine(double&,double,double,int);
  void OpenRootFile();
  void CloseRootFile();
  void InitializeForEcm(double,double); // find all energy-dependent OM terms 
  double SlopeDifference(double,double,int);
  void LevelProperties(double,double,int);
  double FindFermi();
  bool AdjustFermi();
  void WriteAHF();
  double IterateFermi();
  void spectralFunction(double,int,double*,double*);
  double strengthE(double*,double,double,double);
  double energyLab2Cm(double);
  double energyCm2Lab(double);
  void load();
  void loadOM();
  void getTransferWaveFct();

  scat *scatter;

  volume Volume; 
  surfaceTF SurfaceLow;
  surfaceSTD SurfaceHigh;
  
  hartreeFock HartreeFock;
  spinOrbit SpinOrbit;

  double Zp; //charge of projectile 
  double Z;
  double A;
  double asymmetry; // (N-Z)/A
  double sign;
  double Rc;
  double Efermi;
  double gap; //!< gap of particle sepcies under consideration
  double gapOther; //!< gap for other particle species (proton or neutron)
  double gapMin; //!< minimum of proton and neutron gaps
  double Wstart; //!< energy above and below Efermi where Imaginary pot starts
  double Epvolume;

  int mvolume;

  int Ndim;

  double VHF;
  double alpha;
  double beta;
  double gamma;
  double alphaS;
  double betaS;
  double gammaS;
  double coulShift;
  double RHF;
  double aHF;
  double fGap;
  double Asurface;
  double Bsurface;
  double Csurface;
  double Dsurface;

  double Rsurface;
  double asurface;

 
  double Avolume;
  double Bvolume;
  double deltaRvolume;
  double expRvolume;
  double Rvolume;
  double evolume;
  double avolume;
  int AsyVolume;
  double alphaVolume;
  double EaVolume;

  double aHigh;
  double bHigh;
  double cHigh;

  double Vso;
  double VspinoE; // energy dependent term
  double Rso;
  double aso;
  double AWso;
  double BWso;

  levels Level[NlevelMax];
  int Nlevel;
  int NFitLevels; // number of levels in fit

  levelTh LevelTh[NlevelMax];
  levelTh LevelThSort[NlevelMax];
  int NlevelTh;
  double Ef; // calculated Fermi energy
  moment Moment[100];
  int Nmoment;


  // spectrospopic function
  int Nsf; // number of energy in array
  double Elow; // for electon scattering lower limits for integrated strength
  double Ehigh; // upper limit
  double *Esf; //array of energies


  ofstream Fout;
#ifdef root
  TFile *f;
#endif


 private:
  bool bprint;

};

#endif
