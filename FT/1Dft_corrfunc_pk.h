/*==========================================================================*/
/* version: 0.01, 13/07/2011                                                */       
/* author: Francesco Montesano (montefra@mpe.mpg.de)                        */
/*                                                                          */
/* Purpose: compute the 1D fourier transform of the power spectrum into the */
/* correlation function or viceversa using gsl integration routines         */
/* P(k) = 4\pi \int_0_\infty dr r^2 xi(r) j_0(rk)                           */
/* xi(r) = 1/(2\pi^2) \int_0_\infty dk k^2 P(k) j_0(rk)                     */
/*==========================================================================*/

#include <iostream>
#include <stdlib.h>
#include <string>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_spline.h>

using namespace std;

struct gslfparams{
  gsl_spline *spl;   // declaration of the interpolation scheme
  gsl_interp_accel *acc;   // declaration of the accelarator for the interpolation;
  double a,b;     //extrema of the integration
};

/*==========================================================================*/
/* Fourier transform the power spectrum into the correlation function       */
/* 'k' and 'pk' are the input wavenumber and P(k) with dimention 'n' or 'nk'*/
/* 'xi' contains the output correlation function.                           */
/* all pointers have to point to already allocated vectors                  */
/*==========================================================================*/
double* pk2corfunc(double* k, double* pk, double* xi, size_t n);  //xi(r) must have 'n' elements and is evaluated in r=\pi/k
void pk2corfunc(double* k, double* pk, size_t nk, double* r, double* xi, size_t nr);   //'xi' is evaluated at the given 'r' scales; 'xi' and 'r' must have 'nr' elements

/*==========================================================================*/
/* Fourier transform the correlation function into the power spectrum       */
/* 'r' and 'xi' are the input scale and xi(r) with dimension 'n' or 'nr'    */
/* 'pk' contains the output correlation function.                           */
/* all pointers have to point to already allocated vectors                  */
/*==========================================================================*/
double* corfunc2pk(double* r, double* xi, double* pk, size_t n);  //'pk' must have 'n' elements and is evaluated in k=\pi/r
void corfunc2pk(double* r, double* xi, size_t nr, double* k, double* pk, size_t nk);   //'pk' is evaluated at the given 'k' wavenumbers; 'Pk' and 'k' must have 'nk' elements

/*==========================================================================*/
/* compute f(x) = int_0^\infty dk k^2 fk(k) j_0(k*x)                        */
/* j_0(k*x) is the spherical bessel function of order 0                     */
/* it uses the gsl function 'gsl_integration_qawo' to compute the integral  */
/* f(x) = int_k[0]^k[nk-1] dk k/x fk(k) sin(k*x)                            */
/* 'nx' is the size of 'x' and 'fx' and 'nk' is the size of 'k' and 'fk'    */
/* 'fx' needs to be already allocated                                       */
/*==========================================================================*/
void ftQawf(double* k, double* fk, size_t nk, double* x, double* fx, size_t nx);

/*==========================================================================*/
/* function to be passed to the gsl integration routine                     */
/*==========================================================================*/
double function(double x, void* params);


/*==========================================================================*/
/* misc functions                                                           */
/*==========================================================================*/

/* allocate a qawo table and check for errors                               */
/* for the defination of the first four quantities see gsl manual           */
/* the last two are name of the file and the number of the line where the   */
/* function is called: it's use for the error message                       */
/* can be given with the __FILE__ and __LINE__ preprocessor definitions     */
gsl_integration_qawo_table* alloc_qawo_table(double omega, double L, enum gsl_integration_qawo_enum sine,
                                             size_t n, string filename, int linenumber);

/* allocate a integration workspace and check for errors                    */
/* for the defination of the first quantity see gsl manual                  */
/* the last two are name of the file and the number of the line where the   */
/* function is called: it's use for the error message                       */
/* can be given with the __FILE__ and __LINE__ preprocessor definitions     */
gsl_integration_workspace* alloc_integration_workspace(size_t n, string filename, int linenumber);
