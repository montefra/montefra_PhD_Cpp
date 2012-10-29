/*==========================================================================*/
/* Author: Francesco Montesano, montefra@mpe.mpg.de                         */
/* 06/02/2012,    version: 0.01                                             */
/*                                                                          */
/* Corrections to the MAS to apply to the square of the density in order    */
/* to correct for the FFT alias                                             */
/*==========================================================================*/

#ifndef MAS_CORRECTION_H
#define MAS_CORRECTION_H

#include <cstdlib>
#include <cmath>

#include <gsl/gsl_sf_trig.h>

/*==========================================================================
 * Virtual base class to create a pointer to the needed class
 * containing the MAS correction
 *==========================================================================*/
class correction_base{
  public:
    virtual double corr(double kx, double ky, double kz) = 0;  //prototipe of the correction
};

/*==========================================================================
 * No correction wanted, 'corr' returns 1
 *==========================================================================*/
class no_correction: public correction_base{
  public:
    /*==========================================================================
     * Correction to the MAS
     * Parameters
     * ----------
     * kx, ky, kz: wavenumbers where to evaluate the correction
     * output
     * ------
     * corr: 1.
     *==========================================================================*/
    double corr(double kx, double ky, double kz){return 1.;}
};

/*==========================================================================
 * class with the definition of the correction 
 * [prod_i W^2(k_i)]^power, with W(k_i)=sinc(pi*k_i/[2*kN]), i=x,y,z
 * and power=1,2,3 for NGP, CIC and TSC
 *==========================================================================*/
class W_k_correction: public correction_base{
  private:
    double kN;   //Nyquist frequency * 2
    int power;   //1: NGP; 2: CIC; 3: TSC
  public:
    /*==========================================================================
     * Class constructor
     * Parameters
     * ----------
     * kNf: Nyquist wavenumber
     * p: power of the correction. p=1,2,3 for NGP, CIC and TSC
     *==========================================================================*/
    W_k_correction(double kNf, int p){ 
      kN = kNf*2.; 
      power = p;
    }
    /*==========================================================================
     * Correction to the MAS
     * Parameters
     * ----------
     * kx, ky, kz: wavenumbers where to evaluate the correction
     * output
     * ------
     * corr: [prod_i W^2(k_i)]^power, with W(k_i)=sinc(pi*k_i/[2*kN]), i=x,y,z
     *   and power=1,2,3 for NGP, CIC and TSC
     *==========================================================================*/
    double corr(double kx, double ky, double kz){
      double wk = gsl_sf_sinc(kx/kN) * gsl_sf_sinc(ky/kN) * gsl_sf_sinc(kz/kN);
      return( pow(wk, 2.*power) );
    }
};

/*==========================================================================
 * class with the definition of the correction for CIC
 * W^2(K+2NkN) = [prod_i W^2(k_i)], with W(k_i)=1-2/3*sin^2(pi*k_i/[2*kN])
 * and i=x,y,z
 *==========================================================================*/
class W_k_2nkN_correction_CIC: public correction_base{
  private:
    double kN;   //Nyquist frequency*2/pi

    /*==========================================================================
     * approximate one dimensional correction
     * Parameters
     * ----------
     * k: wavenumber where to evaluate the correction
     * output
     * ------
     * w_k: 1-2/3*sin^2(k)
     *==========================================================================*/
    double w_k(double k) { return( 1. - (2./3.) * pow(sin(k), 2.) ); }
  public:
    /*==========================================================================
     * Class constructor
     * Parameters
     * ----------
     * kNf: Nyquist wavenumber
     *==========================================================================*/
    W_k_2nkN_correction_CIC(double kNf){ kN = kNf*2./M_PI; }
    /*==========================================================================
     * Correction to the MAS
     * Parameters
     * ----------
     * kx, ky, kz: wavenumbers where to evaluate the correction
     * output
     * ------
     * corr: [prod_i W^2(k_i)], with W(k_i)=1-2/3*sin^2(pi*k_i/[2*kN]), 
     *   for i=x,y,z
     *==========================================================================*/
    double corr(double kx, double ky, double kz){
      return( w_k(kx/kN) * w_k(ky/kN) * w_k(kz/kN) );
  }
};
/*==========================================================================
 * class with the definition of the correction for TSC
 * W^2(K+2NkN) = [prod_i W^2(k_i)], with 
 * W(k_i)=1-sin^2(pi*k_i/[2*kN]) + 2/15*sin^4(pi*k_i/[2*kN])
 * and i=x,y,z
 *==========================================================================*/
class W_k_2nkN_correction_TSC: public correction_base{
  private:
    double kN;   //Nyquist frequency*2/pi

    /*==========================================================================
     * approximate one dimensional correction
     * Parameters
     * ----------
     * k: wavenumber where to evaluate the correction
     * output
     * ------
     * w_k: 1-2/3*sin^2(k)
     * w_k: 1-sin^2(k) + 2/15*sin^4(k)
     *==========================================================================*/
    double w_k(double k) { 
      return( 1. - pow(sin(k), 2.) + (2./15.) * pow(sin(k), 4.) );
    }
  public:
    /*==========================================================================
     * Class constructor
     * Parameters
     * ----------
     * kNf: Nyquist wavenumber
     *==========================================================================*/
    W_k_2nkN_correction_TSC(double kNf){ kN = kNf*2./M_PI; }
    /*==========================================================================
     * Correction to the MAS
     * Parameters
     * ----------
     * kx, ky, kz: wavenumbers where to evaluate the correction
     * output
     * ------
     * corr: [prod_i W^2(k_i)], with W(k_i)=1-2/3*sin^2(pi*k_i/[2*kN]), 
 * W(k_i)=1-sin^2(pi*k_i/[2*kN]) + 2/15*sin^4(pi*k_i/[2*kN])
     *   for i=x,y,z
     *==========================================================================*/
    double corr(double kx, double ky, double kz){
      return( w_k(kx/kN) * w_k(ky/kN) * w_k(kz/kN) );
    }
};

#endif

