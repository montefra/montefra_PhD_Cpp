/*==========================================================================*/
/* version: 0.01, 19/07/2011                                                */       
/* author: Francesco Montesano (montefra@mpe.mpg.de)                        */
/*                                                                          */
/* Purpose: collection of function to transform from vectors to arrays      */
/* or viceversa                                                             */
/*==========================================================================*/

#include "vector_array.h"

/*==========================================================================*/
/* fill 'array' of size 'n' with values between 'amin' and 'amax'           */
/* the bin value is associated to the middle of the bin                     */
/* if 'isLin' == true linear binning, otherwise logarithmic                 */
/* if 'useExtremes==true' the extrems of the 'n'-1 bins returned            */
/*==========================================================================*/
void fillarray(double* array, size_t n, double amin, double amax, bool isLin, bool useExtremes)
{
  size_t i; // loop integer
  double binsize;   //bin size
  double shift;  //shift of the bin value 
  size_t n_bins;  //number of bins

  if(useExtremes == true){   //if the 'n' extremes of the 'n-1' bins required
    shift = 0.;    //set the shift to 0
    n_bins = n-1;   //reduce by 1 the number of bins
  }
  else{
    shift = 0.5;    //set the shift to 0.5
    n_bins = n;   //reduce by 1 the number of bins
  }

  if(isLin == true){
    binsize = (amax-amin)/n_bins;
    for(i=0; i<n; ++i) array[i] = binsize * (double(i)+shift) + amin;
  }
  else{
    binsize = gsl_sf_log(amax/amin)/n_bins;
    amin = gsl_sf_log(amin);
    for(i=0; i<n; ++i) array[i] = exp(binsize * (double(i)+shift) + amin);
  }
}

/*==========================================================================*/
/* convert vector of doubles 'in' into a corresponding array, allocating it */
/* the size of the array is returned                                        */
/*==========================================================================*/
size_t vector2array(vector<double> &in, double* &out)
{
  size_t i; //loop integer
  size_t insize;  //size of the input vector

  insize = in.size();  //get the size
  try{
    out = new double[insize];   //allocate space of the 'out' array
  }
  catch(bad_alloc& e){
    cerr << "Error allocating the array 'out'. Message: " << e.what() << endl;
  }
  for(i=0; i<insize; ++i) out[i] = in[i];  //copy 'in' in 'out'
  return(insize);   //return the size
}

size_t vector2array(vector<double> &in, double* &out, bool clear)
//same as before, but 'in' is cleared if 'clear'==true
{
  size_t insize = vector2array(in, out);   //convert the vector in array and gets back the size
  if(clear == true) in.clear();   //clear the vector
  return(insize);   //return the size
}

/*==========================================================================*/
/* shift the content of 'array' of size 'n' of 'shift' elements             */
/* if 'shift>0' rightward shift and the last element goes in front          */
/* if 'shift<0' leftward shift and the first element goes in the end        */
/*==========================================================================*/
void shiftarray(double* array, size_t n, ptrdiff_t shift)
{
  size_t i;   //loop integer
  double* temp;   //temporary array

  if(shift > 0){
    allocarray(temp, shift, __FILE__, __LINE__); // allocate the temporary array
    n -= shift;   //take out the shifted components from 'n'
    for(i=0; i<shift; ++i) temp[i] = array[i+n];   //copy the last 'shift' elements into the temporary array
    for(i=0; i<n; ++i) array[n-i+shift-1] = array[n-i-1];   //shift forward the elements in the array
    for(i=0; i<shift; ++i) array[i] = temp[i];    //copy the last 'shift' elments at the beginning of the array
  }
  else{
    allocarray(temp, -shift, __FILE__, __LINE__); // allocate the temporary array
    for(i=0; i<-shift; ++i) temp[i] = array[i];   //copy the first 'shift' elements into the temporary array
    n += shift;   //take out the shifted components from 'n'
    for(i=0; i<n; ++i) array[i] = array[i-shift];   //shift backwards the elements in the array
    for(i=0; i<-shift; ++i) array[i+n] = temp[i];    //copy the first 'shift' elments at the end of the array
  }
  delete[] temp;   //release the memory
}




