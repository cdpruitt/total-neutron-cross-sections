#ifndef _channel
#define _channel

#include "level.h"
#include <string>
#include <fstream>
#include <iostream>
#include <cmath>
#include <cstdlib>
#include "levelD.h"

using namespace std;


/**
 *\brief levels and transmission coefficients for a decay branch
 *
 *reads in and store the daughter levels and the transmission coefficients
 *for a particular decay of a compound nucleus
 */

class channel
{
 public:
  channel(string*,string*,double);
  ~channel();
  double TransCoef(double,int,double);
  double FermiGas(double,int);
  int Ntll;

  // target states
  level *Level;
  int Nlevel;
  double Qvalue; // separation energy for this channel
  double threshold;
  levelD * ld;

 private:
  //transmission coeff
  int Ntle;
  double **tl_up;
  double **tl_down;
  double *Earray;

  double littlea;
  double backshift;

  double sigma2t; //sigma**2 divided by temp
  int A;

};

#endif
