#ifndef _level 
#define _level

/**
 *\brief level information structure, energy,spin,parity
 *
 * Structure to store level information: energy, spin, particle of excited
 * states in the compound and daughter nuclei.
 */

class level
{
 public:
  double spin;
  int parity;
  double energy;
  level();
  level(double,int,double);
};
#endif
