/*==========================================================================*/
/* version: 0.01, 19/07/2011                                                */       
/* author: Francesco Montesano (montefra@mpe.mpg.de)                        */
/*                                                                          */
/* Purpose: collection of function to transform from vectors to arrays      */
/* or viceversa                                                             */
/*==========================================================================*/

#include <cmath>
#include <fstream>
#include <iostream>
#include <stdlib.h>
#include <vector>

//gsl libraries
#include <gsl/gsl_sf_log.h>  

using namespace std;

/* allocate 'array' with dimension 'n'                                      */
template <typename T>  void allocarray(T* &array, size_t n, string filename, int linenumber);
#ifndef VECTOR_ARRAY_H
#define VECTOR_ARRAY_H
template <typename T>  //set template name
void allocarray(T* &array, size_t n, string filename, int linenumber)
{
  try {array = new T[n];}
  catch(bad_alloc& e) {
    cerr << filename << ": " << linenumber << ". Error allocating the array" << endl;
    exit(33);
  }
}
#endif

/* fill 'array' of size 'n' with values between 'amin' and 'amax'           */
/* if 'isLin' == true linear binning, otherwise logarithmic                 */
/* if 'useExtremes==true' the extrems of the 'n'-1 bins returned            */
void fillarray(double* array, size_t n, double amin, double amax, bool isLin, bool useExtremes=false);

/* convert vector of doubles 'in' into a corresponding array, allocating it */
/* the size of the array is returned                                        */
size_t vector2array(vector<double> &in, double* &out);
size_t vector2array(vector<double> &in, double* &out, bool clear);   //same as before, but 'in' is cleared if 'clear'==true

/* shift the content of 'array' of size 'n' of 'shift' elements             */
/* if 'shift>0' rightward shift and the last element goes in front          */
/* if 'shift<0' leftward shift and the first element goes in the end        */
void shiftarray(double* array, size_t n, ptrdiff_t shift);
