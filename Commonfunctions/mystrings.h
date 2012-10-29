/*==========================================================================*/
/* version: 0.01, 25/07/2011                                                */       
/* author: Francesco Montesano (montefra@mpe.mpg.de)                        */
/*                                                                          */
/* Purpose: collection of string functions                                  */
/*==========================================================================*/

#include <cstring>
#include <iomanip>
#include <sstream>
#include <string>
#include <vector>


using namespace std;

/*==========================================================================*/
/* Convert 'value' to string. Templatized, so accept all standard types     */
/*==========================================================================*/
template <typename T> string to_string(T const& value);
#ifndef MYSTRINGS_H
#define MYSTRINGS_H
template <typename T>  //set template name
string to_string(T const& value) {
    stringstream sstr;
    sstr << value;
    return sstr.str();
}
#endif


/* padd the integer 'num' with zeroes in order to have 'pads' digits        */
string padwithzero(int num, int pads);


/* return the path and file name from string 'str'. the codes does not      */
/* check if 'str' is a pure file name or a pure path. Unix dividers         */
/* the path is given without the trailing '/'                               */
string extractPath(const string& str);
string extractFilename(const string& str);

/*==========================================================================
 * Convert a vector of strings in a char**
 * Parameters
 * ----------
 * vec: vector of strings
 * output
 * ------
 * chars: array of char arrays
 *==========================================================================*/
char ** vector2char(std::vector<std::string> vec);

