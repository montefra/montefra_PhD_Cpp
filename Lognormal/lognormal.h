/*==========================================================================*/
/* version: 1.00, 05/09/2011                                                */       
/* author: Francesco Montesano (montefra@mpe.mpg.de)                        */
/*                                                                          */
/* Purpose: create lognormal mocks in cubic grids                           */
/*==========================================================================*/

#include <cmath>
#include <cstdlib>
#include <ctime>
#include <limits>
#include <iostream>
#include <sys/stat.h>
#include <stdlib.h>
#include <string>
#include <vector>

//gsl libraries
#include <gsl/gsl_math.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>  
#include <gsl/gsl_sf_log.h>  
#include <gsl/gsl_spline.h>

//header of the functions to compute the 1D FT transform of the power spectrum and correlation function
#include <1Dft_corrfunc_pk.h>
//header of the class for the power spectrum which implement FFTW 3D complex to complex FFT
//#include <powspec.h>
#include <fftw_r2c_c2r.h>
//head of the input-output functions
#include <io.h>
//header of functions to convert vectors to arrays and viceversa
#include <vector_array.h>
//standard operations which enclud strings
#include <mystrings.h>
//command line parsing 
#include <tclap/CmdLine.h> 

using namespace std;

/*==========================================================================*/
/* default rmin, rmax, kmin, kmax for the 1D FT xi(r)<->P(k)                */
/* used to create the r and k arrays the first time the FT is used          */
/*==========================================================================*/

#ifdef MAIN_CPP
double drmin=1.e-3, drmax=200.;  //maximum and minimum valued of r
double dkmin=1.e-4, dkmax=100.;  //maximum and minimum valued of k
#else
extern double drmin, drmax;  //maximum and minimum valued of r
extern double dkmin, dkmax;  //maximum and minimum valued of k
#endif

/*==========================================================================*/
/* Error handling classes                                                   */
/*==========================================================================*/
/* custom class derived by exception to throw error reading file 'ifile'    */       
//class inputexception : public exception
//{
//  string s;  //output string
//  public:
//    inputexception(string ifile) : s(ifile) {}  //constructor that takes a single file name
//    ~inputexception () throw() {};   //destructor
//    virtual const char* what() const throw() { return s.c_str(); }
//};
//
//class outputexception : public exception
//{
//  string s;  //output string
//  public:
//    outputexception(string ifile) : s(ifile) {}  //constructor that takes a single file name
//    ~outputexception () throw() {};   //destructor
//    virtual const char* what() const throw() { return s.c_str(); }
//};

/*==========================================================================*/
/* function declarations                                                    */
/*==========================================================================*/

/* get the number of digits of a size_t integer (unsigned)                  */
int getNumberOfDigits(size_t i);

/* reads 'infile' with 'nc' columns and returns 'x' and 'y' arrays of double*/
/* if there are problems it throws and exception                            */
/* the size of the arrays returned                                          */
size_t readandconvert(string ifile, int nc, double* &x, double* &y);

/* compute xi_G(r) = ln(1+xi(r)). The calculation is done inplace           */
/* xi_G(r) has 'n' elements                                                 */
void logxi(double* xi, size_t n);

/* compare double 'a' and 'b' requiring either relative or absolute         */
/* precision 'epsrel' or 'epsabs'. return 'true' if a == b with the         */
/* required precision                                                       */
bool almostEqual(float a, float b, float epsrel=1.e-4, 
                 float epsabs=numeric_limits<double>::epsilon()*1.e+2);
/* check the binning of array 'k' of size 'n'                               */
/* and returns if it is 'lin', 'log' or 'irr'                               */
string check_binning(double* k, size_t n);
/* Compute the width of the k bins in 'k' and returns it in 'dk'.           */
/* Both arrays must be allocated and of size 'n'                            */
/* if 'isKnown' == True the 'binning' is assumed (options: 'lin', 'log'     */
/* or 'irr'). Otherwise the function checks whether the binning is          */
/* linear, logarithmic or irregular. In the last case the bins are computed */
/* as [(k[i]-k[i-1])/2, (k[i+1]-k[i])/2] and the first and last bins are    */
/* simmetric                                                                */
void k2dk(double* k, size_t n, double* dk, bool isKnown, string binning="None");

/* Compute the extremes of the bins in 'k' and returns them in 'extrk'      */
/* Both arrays must be allocated and of size 'n' and 'n'+1                  */
/* if 'isKnown' == True the 'binning' is assumed (options: 'lin', 'log'     */
/* or 'irr'). Otherwise the function checks whether the binning is          */
/* linear, logarithmic or irregular. In the last case the bins are computed */
/* as [(k[i]-k[i-1])/2, (k[i+1]-k[i])/2] and the first and last bins are    */
/* simmetric                                                                */
void k2extr(double* k, size_t n, double* extrk, bool isKnown, string binning="None");


/* when computin the velocity field in each direction n=x,y,z, the k_n is   */
/* is needed. In order to avoid if statements pointer to function that      */
/* returns the index in the nth direction is set before filling the grid    */
/* defined inline in order to avoid function calls                          */
#ifdef MAIN_CPP
ptrdiff_t (*_return_index)(ptrdiff_t, ptrdiff_t, ptrdiff_t);
ptrdiff_t return_index(ptrdiff_t ii, ptrdiff_t jj, ptrdiff_t ll){ return (*_return_index)(ii, jj, ll); }
ptrdiff_t _return_index_ii(ptrdiff_t ii, ptrdiff_t jj, ptrdiff_t ll){ return ii; }
ptrdiff_t _return_index_jj(ptrdiff_t ii, ptrdiff_t jj, ptrdiff_t ll){ return jj; }
ptrdiff_t _return_index_ll(ptrdiff_t ii, ptrdiff_t jj, ptrdiff_t ll){ return ll; }
#endif

