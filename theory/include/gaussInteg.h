

//*********************************************************************
//Numerical integration by Gaussian Quadrature
class gaussInteg
{
 public:
   gaussInteg(int, const double*, const double*);
   gaussInteg(){}
   gaussInteg( const gaussInteg&);
   double x (int);
   double w (int);
   int nPoints;
   const  double *pos;
   const  double *weight;
};
