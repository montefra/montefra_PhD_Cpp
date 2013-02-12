/*==========================================================================*/
/* Author: Francesco Montesano, montefra@mpe.mpg.de                         */
/* 06/02/2012,    version: 0.01                                             */
/*                                                                          */
/* includes and definitions for power spectra calculation codes             */
/*==========================================================================*/

#ifndef COMMON_FUNCTIONS_H
#define COMMON_FUNCTIONS_H

#include <string>
#include <sstream>

/*==========================================================================
 * Convert 'value' to string. Templatized, so accept all standard types     
 * Parameters
 * ----------
 * value: any type that can be assigned to a stringstream
 *   value to convert to string
 * output
 * ------
 * to_string: string
 *   std::string of the input value
 *==========================================================================*/
//#ifdef MAIN_CPP
template <typename T>  //set template name
std::string to_string(T const& value){
    std::stringstream sstr;
    sstr << value;
    return sstr.str();
}

/*==========================================================================
 * Convert double to string with required floating point accuracy
 * Parameters
 * ----------
 * value: double 
 *   value to convert to string
 * p: int
 *   number of digits after floating point
 * output
 * ------
 * to_string: string
 *   std::string of the input value
 *==========================================================================*/
//#ifdef MAIN_CPP
inline std::string to_string(double value, int p){
    std::stringstream sstr;
    sstr.setf(std::ios_base::fixed);
    sstr.precision(p);
    sstr << value;
    return sstr.str();
}


//#else
//template <typename T> std::string to_string(T const& value);
//#endif
/*=======================================================================
 * Template function that checks if 'val' is positive                    
 * Parameters
 * ----------
 * val: any type that support 'greater than' comparison
 *   value to check
 * output
 * ------
 * positive: bool
 *   true if 'val' is positive, false otherwise
 *=======================================================================*/
//#ifdef MAIN_CPP
template<class T>
bool positive(T val){
  if( val> 0 ) return true;
  else return false;
}
//#else
//template<class T> bool positive(T val);
//#endif

/*==========================================================================
 * Check if file exists
 * Parameters
 * ----------
 * fname: string
 *   file name to check
 * output
 * ------
 * fileexists: bool
 *   true if the file exists, false otherwise
 *==========================================================================*/
//#ifdef MAIN_CPP
inline bool fileexists(std::string fname)
{
  std::ifstream check(fname.c_str(), std::ifstream::in);
  if (check.is_open()){
    check.close();
    return(true);
  }
  else return(false);
}
//#else
//bool fileexists(std::string fname);
//#endif


#endif
