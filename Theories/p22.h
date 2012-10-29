/*==========================================================================*/
/* Version: 0.10           date:19/01/12                                    */
/*                                                                          */
/* Author: Francesco Montesano, MPE, Garching                               */
/*                                                                          */
/* Purpose: calculate P22(k) as in equation 29-30                           */
/* in Crocce & Scoccimarro, Phys.Rev.D,77(2008),023533                      */
/*                                                                          */
/* WARNING:distances expressed in unit of h                                 */
/*==========================================================================*/


#include <cmath>
#include <cstdlib>
#include <cstring> 
#include <string>

//gsl libraries
#include <gsl/gsl_spline.h>

//command line parsing 
#include <tclap/CmdLine.h> 
//head of the input-output functions
#include <io.h>
//standard operations which enclude strings
#include <mystrings.h>

#include "p22_class.h"

#ifdef p22_CPP
gsl_spline *spllin;   //gsl spline of the linear and the 1loop ps
gsl_interp_accel *acclin;   //gsl accellation of the linear and the 1loop ps
/*==========================================================================*/
/* read file 'fname', and save the first two columns into a gsl spline      */
/* object that is returned                                                  */
/*==========================================================================*/
inline gsl_spline * th2spl(std::string fname)
{
  std::vector<double> x,y;  //temporary vector to read the file
  gsl_spline * spl;   //output spline
  
  if( readfirst2columns(fname, 2, x, y) == 1 ) 
    throw inputexception(fname);   //read the input file and save in two vectors 

  spl = gsl_spline_alloc(gsl_interp_linear, x.size());
  if(gsl_spline_init(spl, &x[0], &y[0], x.size()) !=0){
    cerr << "Initialization of the spline " << gsl_spline_name(spl);
    cerr << " for the data in file " << fname << " failed" << endl;
    exit(41);
  }

  x.clear();
  y.clear();
  return spl;  //return the spline
}
#else
extern gsl_spline *spllin;   //gsl spline of the linear and the 1loop ps
extern gsl_interp_accel *acclin;   //gsl accellation of the linear and the 1loop ps
gsl_spline * th2spl(std::string fname)
#endif

