/*==========================================================================*/
/* Version: 0.10           date:19/01/12                                    */
/*                                                                          */
/* Author: Francesco Montesano, MPE, Garching                               */
/*                                                                          */
/* Purpose: calculate P22(k) as in equation 29-30                           */
/* in Crocce & Scoccimarro, Phys.Rev.D,77(2008),023533                      */
/*                                                                          */
/*==========================================================================*/

#define P22_CLASS_CPP
#include "p22_class.h"

/*==========================================================================
 * Initialise the points and weights for the gauleg integration
 * Parameters
 * ----------
 * rkmax, rkmin: maximum and minimum k to use in the radial integration
 * rbins: number of bins for the radial integration. Positive: 
 *   approximate linear bins; negative: approximate logarithmic binning
 * abins: number of bins for the angular integration
 *==========================================================================*/
void p22::init_gauleg(double rkmax, double rkmin, int rbins, size_t abins){
  //save the number of bins
  radbins=abs(rbins);
  angbins=abins;

  kr = new double[radbins];   //allocate the values of k for the radial integration
  wr = new double[radbins];   //allocate the values of w for the radial integration
  ct = new double[abins];   //allocate the values of cos(theta) for the angular integration
  wa = new double[abins];   //allocate the values of w for the angular integration
  for(int i=0; i<radbins; ++i){
    kr[i]=0.;
    wr[i]=0.;
  }
  for(size_t i=0; i<angbins; ++i){
    ct[i]=0.;
    wa[i]=0.;
  }
  
  //initialise the abscissa and weights for the radial integration
  if(rbins > 0.){
    gauleg(rkmin, rkmax, kr, wr, radbins);   //linear binning
    exponent = 2.;
  }
  else{
    gauleg(log(rkmin), log(rkmax), kr, wr, radbins);   //logarithmic binning
    for(int i=0; i<radbins; ++i) kr[i] = exp(kr[i]);  //recovert k to the correct values
    exponent = 3.;   //if logarithmic integration, an extra factor kr needed
  }
  gauleg(-1., +1., ct, wa, angbins);  //initialise the abscissa and weights for the angular integration
}

/*==========================================================================
 * Gauss-Legendre decomposition 
 * Parameters
 * ----------
 * x1, x2: minimum and maximum values for the decomposition
 * n: number of bins of the decomposition
 * output
 * ------
 * x: abscissas of the Gauss-Legendre n-point quadrature formula in [x1,x2]
 * w: weights of the Gauss-Legendre n-point quadrature formula
 *==========================================================================*/
void p22::gauleg(double x1, double x2, double *x, double *w, int n){
  int m,j,i;
  double z1,z,xm,xl,pp,p3,p2,p1;

  m=(n+1)/2;
  xm=0.5*(x2+x1);
  xl=0.5*(x2-x1);
  for(i=1;i<=m;i++) {
    z=cos(M_PI*(i-0.25)/(n+0.5));
    do{
      p1=1.0;
      p2=0.0;
      for(j=1;j<=n;j++) {
	p3=p2;
	p2=p1;
	p1=((2.0*j-1.0)*z*p2-(j-1.0)*p3)/j;
      }
      pp=n*(z*p1-p2)/(z*z-1.0);
      z1=z;
      z=z1-p1/pp;
    } while(fabs(z-z1) > EPS);
    x[i-1]=xm-xl*z;
    x[n-i]=xm+xl*z;
    w[i-1]=2.0*xl/((1.0-z*z)*pp*pp);
    w[n-i]=w[i-1];
  }
}

/*==========================================================================
 * integrand of the radial part                                             
 * q**exponent * P(q) * integrand_ang. exponent=2,3 if lin or log binning   
 * Parameters
 * ----------
 * q: wavenumber at which evaluate the radial integrand
 * k: wavenumber at which evaluate the power spectrum
 * output
 * ------
 * angintegr: q**exponent * P(q) * integrand_ang
 *==========================================================================*/
double p22::integrand_rad(double q, double k){
  double angintegr=0;   //integral of the angular part
  double plin=0;   //Plin from the spline evaluation

  for(size_t i=0; i<angbins; ++i) 
    angintegr += wa[i] * integrand_ang(ct[i], q, k);  //do the angular integration

  int err_cod = gsl_spline_eval_e(spll, q, accl, &plin);   //evaluate the interpolation and save the error
  if(err_cod != 0){  //if the interpolation didn't work
    if(err_cod == GSL_EDOM) plin=0.;   //if q is out of the range set plin to 0
    else{  //else exit
      std::cerr << "Error to spline the linear power spectrum in k="<< q << std::endl;
      exit(20);
    }
  }
  return( angintegr * pow(q, exponent) * plin );  // return the integrand
}
/*==========================================================================
 * integrand of the radial part                                             
 * q**exponent * P(q) * integrand_ang. exponent=2,3 if lin or log binning   
 * Parameters
 * ----------
 * q: wavenumber at which evaluate the radial integrand
 * k: wavenumber at which evaluate the power spectrum
 * kstar: gaussian damping scale to apply to the linear power spectrum
 *   when integrating
 * output
 * ------
 * angintegr: q**exponent * P(q) * integrand_ang
 *==========================================================================*/
double p22::integrand_rad(double q, double k, double kstar){
  double angintegr=0;   //integral of the angular part
  double plin=0;   //Plin from the spline evaluation

  for(size_t i=0; i<angbins; ++i) 
    angintegr += wa[i] * integrand_ang(ct[i], q, k, kstar);  //do the angular integration

  int err_cod = gsl_spline_eval_e(spll, q, accl, &plin);   //evaluate the interpolation and save the error
  if(err_cod != 0){  //if the interpolation didn't work
    if(err_cod == GSL_EDOM) plin=0.;   //if q is out of the range set plin to 0
    else{  //else exit
      std::cerr << "Error to spline the linear power spectrum in k="<< q << std::endl;
      exit(20);
    }
  }
  return( angintegr * pow(q, exponent) * plin /* exp(-(q*q/(kstar*kstar)))*/ );  // return the integrand
}

/*==========================================================================
 * integrand of the angular part                                            
 * F2(k,q,t)**2 * P(sqrt(k**2+q**2-2*k*q*t)), with t=cos(theta)             
 * Parameters
 * ----------
 * t: cos of the angle at which evaluate the angular integrand
 * q: wavenumber at which evaluate the radial integrand
 * k: wavenumber at which evaluate the power spectrum
 * output
 * ------
 * angintegr: F2(k,q,t)**2 * P(sqrt(k**2+q**2-2*k*q*t))
 *==========================================================================*/
double p22::integrand_ang(double t, double q, double k){
  double plin=0;   //Plin from the spline evaluation
  double qk = sqrt( q*q + k*k - 2*q*k*t );   //|q-k| 

  int err_cod = gsl_spline_eval_e(spll, qk, accl, &plin);   //evaluate the interpolation and save the error
  if(err_cod != 0){  //if the interpolation didn't work
    if(err_cod == GSL_EDOM) plin=0.;   //if q is out of the range set plin to 0
    else{  //else exit
      std::cerr << "Error to spline the linear power spectrum in k="<< q << std::endl;
      exit(20);
    }
  }
  return( F_2(k,q,t) * plin );
}
/*==========================================================================
 * integrand of the angular part                                            
 * F2(k,q,t)**2 * P(sqrt(k**2+q**2-2*k*q*t)), with t=cos(theta)             
 * Parameters
 * ----------
 * t: cos of the angle at which evaluate the angular integrand
 * q: wavenumber at which evaluate the radial integrand
 * k: wavenumber at which evaluate the power spectrum
 * kstar: gaussian damping scale to apply to the linear power spectrum
 *   when integrating
 * output
 * ------
 * angintegr: F2(k,q,t)**2 * P(sqrt(k**2+q**2-2*k*q*t))
 *==========================================================================*/
double p22::integrand_ang(double t, double q, double k, double kstar){
  double plin=0;   //Plin from the spline evaluation
  double qk = sqrt( q*q + k*k - 2*q*k*t );   //|q-k| 

  int err_cod = gsl_spline_eval_e(spll, qk, accl, &plin);   //evaluate the interpolation and save the error
  if(err_cod != 0){  //if the interpolation didn't work
    if(err_cod == GSL_EDOM) plin=0.;   //if q is out of the range set plin to 0
    else{  //else exit
      std::cerr << "Error to spline the linear power spectrum in k="<< q << std::endl;
      exit(20);
    }
  }
  return( F_2(k,q,t) * plin /* exp(-(qk*qk/(kstar*kstar)))*/ );
}

