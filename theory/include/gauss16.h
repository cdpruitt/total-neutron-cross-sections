#include "gaussInteg.h"

class gauss16 : public gaussInteg
{
 public:
  gauss16();
 private:
  static double const x16[16];
  static double const w16[16];
};

