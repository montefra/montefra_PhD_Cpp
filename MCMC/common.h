/*==========================================================================
 * namespace containing some common function and variables that are shared
 * between all the files of the mcmc
 *==========================================================================*/

#include <sstream>
#include <string>

#pragma once

namespace common{

  /*==========================================================================
   * enable verbose or not. Take care to initialise it in 'main.cpp'
   *==========================================================================*/
#ifdef MAIN
  bool verbose=false;
#else
  extern bool verbose;
#endif

  /*==========================================================================
   * Convert 'value' to string. Templatized, so accept all standard types     
   * Parameters
   * ----------
   *  value: value to convert to string
   * output
   * ------
   *  str: string 
   *==========================================================================*/
  template <typename T>  //set template name
  std::string to_string(T const& value) {
    std::stringstream sstr;
    sstr << value;
    return sstr.str();
  }
  /*==========================================================================
   * Convert string 'str' to number. Templatized, so accept all standard types     
   * Parameters
   * ----------
   *  str: string to convert
   *  err: int
   *    0 if conversion worked, >=1 if failed
   * output
   * ------
   *  result: T
   *==========================================================================*/
  template <typename T>  //set template name
  T to_number(std::string const& str, int *err) {
    std::stringstream sstr(str);
    T result;
    *err = 0;
    if(!(sstr >> result)) //if it fails set err to 1
      *err =1;
    return(result);
  }
}
