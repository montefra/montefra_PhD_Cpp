/*==========================================================================*/
/* Version: 0.10           date:19/01/12                                    */
/*                                                                          */
/* Author: Francesco Montesano, MPE, Garching                               */
/*                                                                          */
/* Purpose: calculate P22(k) as in equation 29-30                           */
/* in Crocce & Scoccimarro, Phys.Rev.D,77(2008),023533                      */
/*                                                                          */
/*==========================================================================*/

#include <iostream>

#include <gsl/gsl_spline.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>

//gauleg precision
#define EPS 3.0e-11

class p22
{
  private:
    double *kr, *ct;   //values of k and of cos(theta) used for the gauleg integration
    double *wr, *wa;   //radial and angular weights used in the gauleg integration
    double exponent;   //exponent of kr in the radial integration
    int radbins;   //number of radial bins
    size_t angbins;  //number of angular bins

    gsl_spline *spll;   //gsl spline of the linear and the 1loop ps
    gsl_interp_accel *accl;   //gsl accellation of the linear and the 1loop ps

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
    void gauleg(double x1, double x2, double *x, double *w, int n);

    /*==========================================================================
     * integrand of the radial part                                             
     * q**exponent * P(q) * integrand_ang. exponent=2,3 if lin or log binning   
     * Parameters
     * ----------
     * q: wavenumber at which evaluate the radial integrand
     * k: wavenumber at which evaluate the power spectrum
     * output
     * ------
     * radintegr: q**exponent * P(q) * integrand_ang
     *==========================================================================*/
    double integrand_rad(double q, double k);
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
     * radintegr: q**exponent * P(q) * integrand_ang
     *==========================================================================*/
    double integrand_rad(double q, double k, double kstar);

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
    double integrand_ang(double t, double q, double k);
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
    double integrand_ang(double t, double q, double k, double kstar);

    /*==========================================================================
     * standard F_2 kernel with k=(0,0,k)
     * Parameters
     * ----------
     * k: wavenumber at which evaluate the power spectrum
     * q: wavenumber at which evaluate the radial integrand
     * t: cos of the angle at which evaluate the angular integrand
     * output
     * ------
     * f: (F_2)^2                                                          
     *==========================================================================*/
    double F_2(double k, double q, double t){
      double f = 4.*k*k*t*t - k*q*t - 3.*q*q;
      f /= 14.*(k*k + q*q - 2.*k*q*t);
      f += 3./14. + k*t/(2.*q);
      return(f*f);
    }

  public:
    /*==========================================================================
     * p22 class constructor: copy spline of the power spectrum to integrate
     * Parameters
     * ----------
     * spl: gsl_spline object containing the power spectrum to integrate
     * acc: gsl_interp_accel of the spline
     *==========================================================================*/
    p22( gsl_spline *spl, gsl_interp_accel *acc ){
      spll = spl;
      accl = acc; 
    }
    /*==========================================================================
     * class destructor
     *==========================================================================*/
    ~p22(){
      delete[] kr;
      delete[] wr;
      delete[] ct;
      delete[] wa;
    }

    /*==========================================================================
     * Initialise the points and weights for the gauleg integration
     * Parameters
     * ----------
     * rkmax, rkmin: maximum and minimum k to use in the radial integration
     * rbins: number of bins for the radial integration. Positive: 
     *   approximate linear bins; negative: approximate logarithmic binning
     * abins: number of bins for the angular integration
     *==========================================================================*/
    void init_gauleg(double rkmax, double rkmin, int rbins, size_t abins);

    /*==========================================================================
     * evaluate P22(k) 
     * Parameters
     * ----------
     * k: wavenumber where to evaluate P22
     * output
     * ------
     * P22: evaluated in k
     *==========================================================================*/
    double computep22(double k){
      double p22=0;
      for(int i=0; i<radbins; ++i) 
        p22 += wr[i] * integrand_rad(kr[i], k);  //do the radial integration
      return(p22/(2.*M_PI*M_PI));
    }
    /*==========================================================================
     * evaluate P22(k) 
     * Parameters
     * ----------
     * k: wavenumber where to evaluate P22
     * kstar: gaussian damping scale to apply to the linear power spectrum
     *   when integrating
     * output
     * ------
     * P22: evaluated in k
     *==========================================================================*/
    double computep22(double k, double kstar){
      double p22=0;
      for(int i=0; i<radbins; ++i) 
        p22 += wr[i] * integrand_rad(kr[i], k, kstar);  //do the radial integration
      return(p22 * exp(-k*k/(kstar*kstar)) / (2.*M_PI*M_PI));
    }

};
