/*==========================================================================*/
/* version: 1.00, 05/09/2011                                                */       
/* author: Francesco Montesano (montefra@mpe.mpg.de)                        */
/*                                                                          */
/* Purpose: functions used in lognormal.cpp                                 */
/*==========================================================================*/

#include "lognormal.h"

/*==========================================================================*/
/* get the number of digits of a size_t integer (unsigned)                  */
/*==========================================================================*/
int getNumberOfDigits(size_t i)
{
      return( i > 0 ? (int) log10((double)i) + 1 : 1 );
}

/*==========================================================================*/
/* reads 'infile' with 'nc' columns and returns 'x' and 'y' arrays of double*/
/* if there are problems it throws and exception                            */
/* the size of the arrays returned                                          */
/*==========================================================================*/
size_t readandconvert(string ifile, int nc, double* &x, double* &y)
{
  size_t insize=0;  //size of the input arrays
  vector<double> tempx, tempy;   //temporary vectors to store the input
  if( readfirst2columns(ifile, nc, tempx, tempy) !=0 )
    throw inputexception(ifile);
  insize = vector2array(tempx, x, true);   //convert the input k from vector to array and clear the vector
  insize = vector2array(tempy, y, true);   //convert the input P(k) from vector to array and clear the vector
  return(insize);
}

/*==========================================================================*/
/* compute xi_G(r) = ln(1+xi(r)). The calculation is done inplace           */
/* xi_G(r) has 'n' elements                                                 */
/*==========================================================================*/
void logxi(double* xi, size_t n)
{
  size_t i;  //loop integer
  for(i=0; i<n; ++i) xi[i] = gsl_sf_log_1plusx(xi[i]);
}

/*==========================================================================*/
/* compare double 'a' and 'b' requiring either relative or absolute         */
/* precision 'epsrel' or 'epsabs'. return 'true' if a == b with the         */
/* required precision                                                       */
/*==========================================================================*/
bool almostEqual(float a, float b, float epsrel, float epsabs)
{
  float relativeError;
  if (fabs(a - b) < epsabs)   //if difference less than absolute error
      return true;
  relativeError = (fabs(b) > fabs(a)) ? fabs((a - b) / b) : fabs((a - b) / a);
  if (relativeError <= epsrel)
      return true;
  return false;
}

/*==========================================================================*/
/* check the binning of array 'k' of size 'n'                               */
/* and returns if it is 'lin', 'log' or 'irr'                               */
/*==========================================================================*/
string check_binning(double* k, size_t n)
{
  size_t i;   //loop integer
  string binning;   //output string
  double dkcomp;   //dk of the first bin to use as comparison

  binning = "lin";  //first guess

  dkcomp = k[1]-k[0];  //binsize for the linear binning
  for(i=1; i<n-1; ++i) if( almostEqual(dkcomp, k[i+1]-k[i]) == false ){   //check if it is linear binning
      binning = "log";
      break;
    }
  if(binning.compare("log") == 0){   //if the binning wasn't linear check if it is logarithmic
    dkcomp = gsl_sf_log(k[1]/k[0]);  //binsize for the logarithmic binning
    for(i=1; i<n-1; ++i) if( almostEqual(dkcomp, gsl_sf_log(k[i+1]/k[i])) == false ){   //check if it is logarithmic binning
	binning = "irr";
	break;
      }
  }
  return(binning);
}
/*==========================================================================*/
/* Compute the width of the k bins in 'k' and returns it in 'dk'            */
/* Both arrays must be allocated and of size 'n'                            */
/* if 'isKnown' == True the 'binning' is assumed (options: 'lin', 'log'     */
/* or 'irr'). Otherwise the function checks whether the binning is          */
/* linear, logarithmic or irregular. In the last case the bins are computed */
/* as [(k[i]-k[i-1])/2, (k[i+1]-k[i])/2] and the first and last bins are    */
/* simmetric                                                                */
/*==========================================================================*/
void k2dk(double* k, size_t n, double* dk, bool isKnown, string binning)
{
  size_t i;   //loop integer
  double* klims;   //contains the limits of the kbins

  if(isKnown == false) binning = check_binning(k, n);   //check the binning if not given

  if( binning.compare("lin") == 0 ){
    allocarray(klims, 1, __FILE__, __LINE__);    //allocate 'klims' with one element
    klims[0] = k[1]-k[0];   //linear binning
    for(i=0; i<n; ++i) dk[i] = klims[0];   //fill the delta k array
    delete[] klims;  //deallocate 'klims'
  }
  else if( binning.compare("log") == 0 ){
    allocarray(klims, n+1, __FILE__, __LINE__);    //allocate 'klims' with one element
    fillarray(klims, n+1, sqrt(k[0]*k[0]*k[0]/k[1]), sqrt(k[n-1]*k[n-1]*k[n-1]/k[n-2]), false);
    for(i=0; i<n; ++i) dk[i] = klims[i+1]-klims[i];   //fill the delta k array
    delete[] klims;  //deallocate 'klims'
  }
  else if( binning.compare("irr") == 0 ){
    dk[0] = k[1] - k[0];   //first bin simmetrix around k[0] with width 2*((k[1]-k[0])/2)
    for(i=1; i<n-1; ++i) dk[i] = (k[i+1]-k[i-1])/2.;   //fill the delta k array between 1 and n-2
    dk[n-1] = k[n-1] - k[n-2];    //last bin simmetrix around k[n-1] with width 2*((k[n-1]-k[n-2])/2)
  }
  else{
    cerr << "Only 'lin', 'log' and 'irr' keywords accepted." << endl;  
    exit(40);
  }
}

/*==========================================================================*/
/* Compute the extremes of the bins in 'k' and returns them in 'extrk'      */
/* Both arrays must be allocated and of size 'n' and 'n'+1                  */
/* if 'isKnown' == True the 'binning' is assumed (options: 'lin', 'log'     */
/* or 'irr'). Otherwise the function checks whether the binning is          */
/* linear, logarithmic or irregular. In the last case the bins are computed */
/* as [(k[i]-k[i-1])/2, (k[i+1]-k[i])/2] and the first and last bins are    */
/* simmetric                                                                */
/*==========================================================================*/
void k2extr(double* k, size_t n, double* extrk, bool isKnown, string binning)
{
  size_t i;   //loop integer

  if(isKnown == false) binning = check_binning(k, n);   //check the binning if not given

  if( binning.compare("lin") == 0 )
    fillarray(extrk, n+1, (3.*k[0]-k[1])/2., (3.*k[n-1]-k[n-2])/2., true);   //linear binning: fill the array with extremes
  else if( binning.compare("log") == 0 )
    fillarray(extrk, n+1, sqrt(k[0]*k[0]*k[0]/k[1]), sqrt(k[n-1]*k[n-1]*k[n-1]/k[n-2]), false);   //logarithmic binning: fill the array with extremes
  else if( binning.compare("irr") == 0 ){
    extrk[0] = (3.*k[0] - k[1])/2.;   //first bin simmetrix around k[0] with width 2*((k[1]-k[0])/2)
    for(i=1; i<n; ++i) extrk[i] = (3.*k[i]-k[i-1])/2.;   //fill the bin edges array between 1 and n-1
    extrk[n] =(3.*k[n-1] - k[n-2])/2.;    //last bin simmetrix around k[n-1] with width 2*((k[n-1]-k[n-2])/2)
  }
  else{
    cerr << "Only 'lin', 'log' and 'irr' keywords accepted." << endl;  
    exit(40);
  }
}
