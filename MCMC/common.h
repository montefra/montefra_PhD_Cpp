/*==========================================================================
 * namespace containing some common function and variables that are shared
 * between all the files of the mcmc
 *==========================================================================*/

#pragma once

#include <fstream>
#include <sstream>
#include <string>

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
   * convert string 'str' to number. templatized, so accept all standard types     
   * parameters
   * ----------
   *  str: string to convert
   *  err: int
   *    0 if conversion worked, >=1 if failed
   * output
   * ------
   *  result: t
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

}


