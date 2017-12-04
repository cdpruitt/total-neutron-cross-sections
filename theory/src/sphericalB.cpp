#include "../include/sphericalB.h"

double const sphericalB::pi = acos(-1.);


/**
 * returns the modified spherical Bessel function of the second kind needed
 * for bound states.
 \param l = order of the function (orbital angular momentum)
 \param rho = independent variable (rho = k * r)
 */
double sphericalB::k(int l, double rho)
{
  switch (l)
    {
    case 0:
      return k0(rho);
    case 1:
      return k1(rho);
    case 2: 
      return k2(rho);
    case 3: 
      return k3(rho);
    case 4:
      return k4(rho);
    case 5:
      return k5(rho);
    case 6:
      return k6(rho);
    case 7:
      return k7(rho);
    default:
      cout << "no l>6 programed in sphericalB" << endl;
      return 0.;
    }
}
//****************************************
  /**
   * returns modified spherical Bessel function of the second kind
   * order = 0 , i.e. l=0
   *
   \param rho = independent variable (rho = k * r)
   */

double sphericalB::k0(double rho)
{
  double out = exp(-rho);
  derivative = -out;
  return out;
}
//***************************************
  /**
   * returns modified spherical Bessel function of the second kind
   * order = 1 , i.e. l=1
   *
   \param rho = independent variable (rho = k * r)
   */
double sphericalB::k1(double rho)
{
 double  out = exp(-rho)*(1.+1./rho);
 derivative = -out - exp(-rho)/pow(rho,2);
 return out;
}
//***************************************
  /**
   * returns modified spherical Bessel function of the second kind
   * order = 2 , i.e. l=2
   *
   \param rho = independent variable (rho = k * r)
   */
double sphericalB::k2(double rho)
{
 double  out = exp(-rho)*(1.+3./rho+3./pow(rho,2));
 derivative = -out + exp(-rho)*(-3./pow(rho,2)-6./pow(rho,3));
 return out;
}
//***************************************
  /**
   * returns modified spherical Bessel function of the second kind
   * order = 3 , i.e. l=3
   *
   \param rho = independent variable (rho = k * r)
   */
double sphericalB::k3(double rho)
{
 double  out = exp(-rho)*(1.+6./rho+15./pow(rho,2)+15./pow(rho,3));
 derivative = -out + exp(-rho)*(-6./pow(rho,2)-30./pow(rho,3)-
                                      45./pow(rho,4));
 return out;
}
//***************************************
  /**
   * returns modified spherical Bessel function of the second kind
   * order = 4 , i.e. l=4
   *
   \param rho = independent variable (rho = k * r)
   */
double sphericalB::k4(double rho)
{
 double  out = exp(-rho)*(1.+10./rho+45./pow(rho,2)+105./pow(rho,3)
+105./pow(rho,4));
 derivative = -out + exp(-rho)*(-10./pow(rho,2)-90./pow(rho,3)-
                                      315./pow(rho,4)-420./pow(rho,5));
 return out;
}
//***************************************
  /**
   * returns modified spherical Bessel function of the second kind
   * order = 5 , i.e. l=5
   *
   \param rho = independent variable (rho = k * r)
   */
double sphericalB::k5(double rho)
{
 double  out = exp(-rho)*(1.+15./rho+105./pow(rho,2)+420./pow(rho,3)
+945./pow(rho,4)+945./pow(rho,5));
 derivative = -out + exp(-rho)*(-15./pow(rho,2)-210./pow(rho,3)-
                                      1260./pow(rho,4)-3780./pow(rho,5)
                                      -4725./pow(rho,6));
 return out;
}
//***************************************
  /**
   * returns modified spherical Bessel function of the second kind
   * order = 6 , i.e. l=6
   *
   \param rho = independent variable (rho = k * r)
   */
double sphericalB::k6(double rho)
{
 double  out = exp(-rho)*(1.+21./rho+210./pow(rho,2)+1260./pow(rho,3)
+4725./pow(rho,4)+10395./pow(rho,5)+10395./pow(rho,6));
 derivative = -out + exp(-rho)*(-21./pow(rho,2)-420./pow(rho,3)-
 3780./pow(rho,4)-18900./pow(rho,5)-51975./pow(rho,6) - 62370./pow(rho,7));
 return out;
}
//***************************************
  /**
   * returns modified spherical Bessel function of the second kind
   * order = 7 , i.e. l=7
   *
   \param rho = independent variable (rho = k * r)
   */
double sphericalB::k7(double rho)
{
 double  out = exp(-rho)*(1.+28./rho+378./pow(rho,2)+3150./pow(rho,3)
+17325./pow(rho,4)+62370./pow(rho,5)+135135./pow(rho,6)+135135/pow(rho,7));
 derivative = -out + exp(-rho)*(-28./pow(rho,2)-756./pow(rho,3)-
 9450./pow(rho,4)-69300./pow(rho,5)-311850./pow(rho,6) - 810810./pow(rho,7) -
 945945./pow(rho,8));
 return out;
}
//****************************************************************
/**
 * returns the spherical Bessel function of the second kind
 * needed for continuum states
 \param l = order of the function (orbital angular momentum)
 \param rho = independent variable (rho = k * r)
 */
double sphericalB::y(int l, double rho)
{
  switch (l)
    {
    case 0:
      return y0(rho);
    case 1:
      return y1(rho);
    case 2: 
      return y2(rho);
    case 3: 
      return y3(rho);
    case 4:
      return y4(rho);
    case 5:
      return y5(rho);
    case 6:
      return y6(rho);
    case 7:
      return y7(rho);
    default:
      cout << "no l>6 programed in sphericalB" << endl;
      return 0.;
    }
}
//****************************************
  /**
   * returns modified spherical Bessel function of the second kind
   * order = 0 , i.e. l=0
   *
   \param rho = independent variable (rho = k * r)
   */
double sphericalB::y0(double rho)
{
  double out = -cos(rho);
  derivative = sin(rho);
  return out;
}
//****************************************
  /**
   * returns modified spherical Bessel function of the second kind
   * order = 1 , i.e. l=1
   *
   \param rho = independent variable (rho = k * r)
   */
double sphericalB::y1(double rho)
{
  double out = -cos(rho)/rho - sin(rho);
  derivative = (cos(rho)-pow(rho,2)*cos(rho)+rho*sin(rho))/pow(rho,2);
  return out;
}
//****************************************
  /**
   * returns modified spherical Bessel function of the second kind
   * order = 2 , i.e. l=2
   *
   \param rho = independent variable (rho = k * r)
   */
double sphericalB::y2(double rho)
{
  double out = cos(rho)*(1.-3./pow(rho,2)) - sin(rho)*3./rho;
  derivative = -(3.*(pow(rho,2)-2.)*cos(rho)+rho*(pow(rho,2)-6.)*sin(rho))
  /pow(rho,3);
  return out;
}
//****************************************
  /**
   * returns modified spherical Bessel function of the second kind
   * order = 3 , i.e. l=3
   *
   \param rho = independent variable (rho = k * r)
   */
double sphericalB::y3(double rho)
{
  double out = cos(rho)*(6./rho-15./pow(rho,3)) + sin(rho)*(1.-15./pow(rho,2));
  derivative = ((45.-21.*pow(rho,2)+pow(rho,4))*cos(rho) + 
		3.*rho*(15.-2.*pow(rho,2))*sin(rho))/pow(rho,4);

  return out;
}
//****************************************
  /**
   * returns modified spherical Bessel function of the second kind
   * order = 4 , i.e. l=4
   *
   \param rho = independent variable (rho = k * r)
   */
double sphericalB::y4(double rho)
{
  double out = cos(rho)*(-1.+45./pow(rho,2)-105./pow(rho,4)) 
    + sin(rho)*(10./rho-105./pow(rho,3));
  derivative = (5.*(84.-39.*pow(rho,2)+2.*pow(rho,4))*cos(rho)+
		rho*(420.-55.*pow(rho,2)+pow(rho,4))*sin(rho))/pow(rho,5);

  return out;
}
//****************************************
  /**
   * returns modified spherical Bessel function of the second kind
   * order = 5 , i.e. l=5
   *
   \param rho = independent variable (rho = k * r)
   */
double sphericalB::y5(double rho)
{
  double out = cos(rho)*(-15./rho+420./pow(rho,3)-945./pow(rho,5))
    + sin(rho)*(-1.+105./pow(rho,2)-945./pow(rho,4));
  derivative = (-(-4725.+2205*pow(rho,2)-120.*pow(rho,4)+pow(rho,6))*cos(rho)
		+15.*rho*(315.-42.*pow(rho,2)+pow(rho,4))*sin(rho))/pow(rho,6);

  return out;
}
//****************************************
  /**
   * returns modified spherical Bessel function of the second kind
   * order = 6 , i.e. l=6
   *
   \param rho = independent variable (rho = k * r)
   */
double sphericalB::y6(double rho)
{
  double out = cos(rho)*(1.-210./pow(rho,2)+4725/pow(rho,4)-10395./pow(rho,6))
    +sin(rho)*(-21./rho+1260./pow(rho,3)-10395./pow(rho,5));
  derivative = -(21.*(-2970.+1395*pow(rho,2)-80.*pow(rho,4)+pow(rho,6))
		 *cos(rho)+rho*(-62370.+8505.*pow(rho,2)-231.*pow(rho,4)+
				pow(rho,6))*sin(rho))/pow(rho,7);

  return out;
}
//****************************************
  /**
   * returns modified spherical Bessel function of the second kind
   * order = 7 , i.e. l=7
   *
   \param rho = independent variable (rho = k * r)
   */
double sphericalB::y7(double rho)
{
  double out = cos(rho)*(28./rho-3150./pow(rho,3)+62370./pow(rho,5)-
    135135./pow(rho,7)) + sin(rho)*(1.-378./pow(rho,2)
				    +17325./pow(rho,4)-135135./pow(rho,6));
  derivative = ((945945. - 446985*pow(rho,2) + 26775*pow(rho,4) -
   406.*pow(rho,6)+pow(rho,8))*cos(rho)+7.*rho*(135135.-18810*pow(rho,2)+
   558.*pow(rho,4)-4.*pow(rho,6))*sin(rho))/pow(rho,8);

  return out;
}
//**********************************************
  /**
   * returns the logarithmic derivative of the k speherical Bessel function
   \param l order of function (orbital angular momentum)
   \param rho independent variable
   */
  double sphericalB::LogDer_k(int l ,double rho)
{
  if (l > 0 && rho == 0.) return -1.e32;
  double out = k(l,rho);
  return derivative/out;
}
//**********************************************
  /**
   * returns the logarithmic derivative of the y speherical Bessel function
   \param l order of function (orbital angular momentum)
   \param rho independent variable
   */
  double sphericalB::LogDer_y(int l ,double rho)
{
  if (l > 0 && rho == 0.) return -1.e32;
  double out = y(l,rho);
  return derivative/out;
}
