/*==========================================================================*/
/* Version: 4.00           date:17/01/12                                    */
/*                                                                          */
/* Author: Francesco Montesano, MPE, Garching                               */
/*                                                                          */
/* Purpose : includes and definitions for the main in 'mcmc.cpp'            */
/*==========================================================================*/

#include <cmath>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <limits>
#include <iostream>
#include <sys/stat.h>
#include <sstream>
#include <stdlib.h>
#include <string>
#include <vector>

/* use of machine time in second as random seed */
#include <time.h>   
#include <ulimit.h>
// GSL libraries
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_errno.h>
//head of the input-output functions
#include <io.h>
#include <mystrings.h>

//command line parsing 
#include <tclap/CmdLine.h> 

//include the classes for the mode coupling power spectrum
#include "pMC.h"

/*==========================================================================*/
/* structures                                                               */
/*==========================================================================*/
struct in_data{   //structure used to read part of the input files
  double x, y, z;
};
struct cosmology{  //this structure contains parameters which depends only on cosmology
  double f;   //dln(D(z))/dln(a)
};

/*==========================================================================*/
/* Common variables                                                         */
/*==========================================================================*/

#ifdef MCMC_CPP
gsl_vector *data, *kj, *G2, *W0j;     //P(k), kj, G^2(ki), W0j
double G20;    //G^2(0)
gsl_matrix *invcov, *window; //inversed of the covariance matrix and window matrix Wij
gsl_spline *spllin;   //gsl spline of the linear and the 1loop ps
gsl_interp_accel *acclin;   //gsl accellation of the linear and the 1loop ps
gsl_vector *theory, *wtheory;  // contains the theoretical PS and the theoretical PS convolved with the window function
gsl_vector *tempv;   // temporary vector
gsl_matrix *tempm, *tempm2;   // temporary matrices
cosmology cosmo;   //cosmology structure
#else
extern gsl_vector *data, *kj, *G2, *W0j;     //P(k), kj, G^2(ki), W0j
extern double G20;    //G^2(0)
extern gsl_matrix *invcov, *window; //inversed of the covariance matrix and window matrix Wij
extern gsl_spline *spllin;   //gsl spline of the linear and the 1loop ps
extern gsl_interp_accel *acclin;   //gsl accellation of the linear and the 1loop ps
extern gsl_vector *theory, *wtheory;  // contains the theoretical PS and the theoretical PS convolved with the window function
extern gsl_vector *tempv;   // temporary vector
extern gsl_matrix *tempm, *tempm2;   // temporary matrices
extern cosmology cosmo;   //cosmology structure
#endif

/*==========================================================================*/
/* Functions declaration                                                    */
/*==========================================================================*/

int input_data(std::string input_file, std::vector<in_data> &indata, int which);
gsl_matrix * winmat(std::string winfile, int kimin, int kjmin, int nki, int nkj, int sizei, int sizej);
gsl_spline * th2spl(std::string fname);

/*==========================================================================*/
int paramnames(std::string ofile, std::vector<std::string> &sname);
void readini(std::string inifile, std::vector<std::string> sname, 
    double *params, double *upper, double *lower, double *sigma);
/*==========================================================================
 * Substitute the a value in the input array with the desired one
 * Parameters
 * ----------
 * array: pointer to double 
 *   array where to modify
 * n: size_t
 *   size of array
 * i: size_t
 *   position to substitute
 * val: double
 *   array[i] = val
 * output
 * ------
 * array: pointer to double
 *   pointer to the array with the substituted value
 *==========================================================================*/
double * substitute_value( double *array, size_t n, size_t i, double val );

/*==========================================================================*/
gsl_vector * vec2gslvec(std::vector<in_data> &indata, double *kmin, double
    kmax, int x_or_y, bool mod_kmin);
gsl_matrix * invert(std::vector<in_data> &indata, int diminv, int mininv, 
    int n_mocks, bool diag);

double likelihood(double *par, pmc *p22);
double p_obj(double k, double *par, pmc *p22);
/*==========================================================================
 * angular average redshift space distortions: analytical expression
 * from Ariel.
 * it includes the Kaiser boost factor and finger of god.
 * Parameters
 * ----------
 * k: double
 *   wave number 
 * par: array of doubles
 *   parameters for the step of the mcmc
 * cosmo: cosmology struct
 *   cosmology related parameters. In this case f = ln(D)/ln(a)   
 * output
 * ------
 * rsd: double
 *   value of the correction for the given set of parameters and k
 *==========================================================================*/
double rsd( double k, double *par, cosmology c);

void dealloc_comm(); // deallocate common variables

