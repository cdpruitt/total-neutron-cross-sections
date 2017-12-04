#include <string>
#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
using namespace std;

/**
 *!\brief level densities from the HF-BCS model
 *
 * this class gives an interpolation of spin-dependent level
 * densities from the model of S. Goriely 
 *References
 *[1] S. Goriely, F. Tondeur, J. M. Pearson (2001), Atomic Data Nuclear Data
 *   Tables 77, 311.
 *[2] P. Demetriou and S. Goriely (2001), Nucl. Phys. A695, 95.
 *[3] M. Arnould and F. Tondeur (1981) Proc. Conf. on Nuclei far from stability,
 *   Helsingor, CERN, 81-09, vol.1, p.229.
 *[4] S. Goriely (1996) Nucl. Phys. A605, 28.
 *[5] S. Goriely (1997) Proc. Int. Conf. on Nuclear Data for Science and 
 *   Technology (eds. G. Reffo et al., Italian Physical Society, Italy), 811.
 * the input data files for each nucleus can be cut and pasted from the 
 * RIPL2 site http://www-nds.iaea.org/RIPL-2/
 *
 */

class levelD
{
 private:
  static int const nE = 55; //!< number of excitation energy bins
  static int const nJ = 30; //!< number of spin bins
  float ld[nE][nJ];//!< level densities as a function of excitation energy and spin
  float Ex0[nE]; // !< array of excitation energies
  bool hasLevelD = false;

 public:
  levelD(string filename);
  double getLD(double Ex, int J);
};
