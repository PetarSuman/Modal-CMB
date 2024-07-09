#include <math.h>

#define pi M_PI
#define gamma_one 2.6789385347
#define gamma_two 1.3541179394 

double spherical_bessel ( int l , double z )
{
  double  nu = l + 0.5;
  double  nu2 = nu*nu;
  double  nu3 = nu2*nu;
  double  nu4 = nu2*nu2;
  double  mu = pow(l,0.325);
  double  lower_limit=nu - 1.31*mu;
  double  upper_limit=nu + 1.48*mu;
  double value;

// *************** small z case ************** 
  if (z < 1.0e-3) {
    if (l==0) {
      return 1.0;
    }
    return 0.0;
// *************** small l case **************
  } else if (l < 6) {
    double  sin_z = sin (z);
    double  cos_z = cos (z);
    double  z1 = 1/z;
    double  z2 = z1*z1;
    double  z3 = z2*z1;
    double  z4 = z2*z2;
    double  z5 = z3*z2;
    double  z6 = z3*z3;
    
    switch (l) {
    case 0:
      value = z1*sin_z;
      break;
    case 1:
      value = z2*sin_z - z1*cos_z;
      break;
    case 2:
      value = (3.0*z3 - z1)*sin_z - 3.0*z2*cos_z;
      break;
    case 3:
      value = (15.0*z4 - 6.0*z2)*sin_z - (15.0*z3  - z1)*cos_z;
      break;
    case 4:
      value = (105.0*z5 - 45.0*z3 + z1)*sin_z - (105.0*z4 - 10.0*z2)*cos_z;
      break;
    case 5:
      value = (945.0*z6 - 420.0*z4 + 15.0*z2)*sin_z - (945.0*z5 - 105.0*z3 + z1)*cos_z;
      break;
    }
// *************** z << l case **************
  } else if (z < lower_limit) {
    double z2 = z*z;
    double beta = log(nu/z + sqrt(pow(nu/z,2) - 1.0) );
    double coth = 1.0  / sqrt(1.0 - z2/nu2);
    double coth3 = coth*coth*coth;
    double coth6 = coth3*coth3;
    double sech = z/nu;
    double sech2 = sech*sech;
    double sech4 = sech2*sech2;
    double sech6 = sech2*sech4;
    double prefactor = sqrt(coth/sech)/(2.0*nu);

    double sum1 = 2.0+3.0*sech2;
    double sum2 = 4.0+sech2;
    double sum3 = 16.0-1512.0*sech2-3654.0*sech4-375.0*sech6;
    double sum4 = 32.0+288.0*sech2+232.0*sech4+13.0*sech6;

    double exp_term;
    exp_term = sum1*coth3/(24.0*nu);
    exp_term = exp_term - sum2*sech2*coth6/(16.0*nu2);
    exp_term = exp_term - sum3*coth3*coth6/(5760.0*nu3);
    exp_term = exp_term - sum4*sech2*coth6*coth6/(128.0*nu4);
    exp_term = exp(-nu*beta+nu/coth-exp_term);

    value = prefactor*exp_term;
// *************** z << l case **************
  } else if (z > upper_limit) {
    double z2 = z*z;
    double beta = acos(nu/z);
    double coth = nu / (z*sqrt(1.0 - nu2/z2));
    double coth3 = coth*coth*coth;
    double coth6 = coth3*coth3;
    double sech = z/nu;
    double sech2 = sech*sech;
    double sech4 = sech2*sech2;
    double sech6 = sech2*sech4;
    double prefactor = sqrt(coth/sech)/nu;

    double sum1 = 2.0+3.0*sech2;
    double sum2 = 4.0+sech2;
    double sum3 = 16.0-1512.0*sech2-3654.0*sech4-375.0*sech6;
    double sum4 = 32.0+288.0*sech2+232.0*sech4+13.0*sech6;

    double exp_term;
    exp_term = sum2*sech2*coth6/(16.0*nu2);
    exp_term = exp_term - sum4*sech2*coth6*coth6/(128.0*nu4);
    exp_term=exp(-exp_term);

    double trig_arg;
    trig_arg = nu/coth - nu*beta - pi/4.0;
    trig_arg = trig_arg - sum1*coth3/(24.0*nu);
    trig_arg = trig_arg - sum3*coth3*coth6/(5760.0*nu3);
    double trig_cos;
    trig_cos = cos(trig_arg);

    value = prefactor*exp_term*trig_cos;
    
// ***************  x near l   **************

  } else {
    double beta = z-nu;
    double beta2 = beta*beta;   
    double beta4 = beta2*beta2;
    double beta6 = beta4*beta2;
    double sz = 6.0/z;
    double sz2 = sz*sz;
    double sz3 = sz2*sz;
    double sech = pow(sz,1.0/3);
    double sech2 = sech*sech;
    double prefactor = sqrt(sz/pi) / 12.0;


    double sum1 = (beta2/6.0 - 1.0/15.0)*beta;
    double sum2 = beta4/24.0 - beta2/24.0 + 1.0/280.0;
    double sum3 = beta6/720.0 - 7.0*beta4/1440.0 + beta2/288.0 - 1.0/3600.0;
    double sum4 = (beta6/5040.0 - beta4/900.0 + 19.0*beta2/12600.0 - 13.0/31500.0)*beta;
    double sum5 = (beta4*beta4/362880.0 - beta6/30240.0 + 71.0*beta4/604800.0 - 121.0*beta2/907200.0 + 7939.0/232848000.0)*beta;

    double deriv;
    deriv = gamma_one*sech;
    deriv = deriv+ beta*gamma_two*sech2;
    deriv = deriv - sum1*sz*sech*gamma_one/3.0;
    deriv = deriv - 2.0*sum2*sz*sech2*gamma_two/3.0;
    deriv = deriv + 4.0*sum3*sz2*sech*gamma_one/9.0;
    deriv = deriv + 10.0*sum4*sz2*sech2*gamma_two/9.0;
    deriv = deriv - 28.0*sum5*sz3*sech*gamma_one/27.0;

    value = prefactor*deriv;
  }
  return value;
}
