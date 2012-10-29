/*==========================================================================*/
/* version: 0.01, 25/07/2011                                                */       
/* author: Francesco Montesano (montefra@mpe.mpg.de)                        */
/*                                                                          */
/* Purpose: collection of string functions                                  */
/*==========================================================================*/

#include "mystrings.h"

/*==========================================================================*/
/* padd the integer 'num' with zeroes in order to have 'pads' digits        */
/*==========================================================================*/
string padwithzero(int num, int pads)
{
  ostringstream ss;
  ss << setw( pads ) << setfill( '0' ) << num;
  return ss.str();
}

/*==========================================================================*/
/* return the path and file name from string 'str'. the codes does not      */
/* check if 'str' is a pure file name or a pure path. Unix dividers         */
/* the path is given without the trailing '/'                               */
/*==========================================================================*/
string extractPath(const string& str)
{
return str.substr( 0, str.find_last_of( '/' ) );
}
string extractFilename(const string& str)
{
return str.substr( str.find_last_of( '/' ) +1 );
}

/*==========================================================================
 * Convert a vector of strings in a char**
 * Parameters
 * ----------
 * vec: vector of strings
 * output
 * ------
 * chars: array of char arrays
 *==========================================================================*/
char** vector2char(std::vector<std::string> vec){
  size_t svec = vec.size();   //get the number of elements in the vector
  char **chars = new char*[svec];   //allocate the array of char pointers
  for(size_t i=0; i<svec; ++i){     //loop over the vector elements
    chars[i] = new char[vec[i].size()+1];   //allocate enought space to copy each string
    std::strcpy( chars[i], vec[i].c_str() );  //copy the string to the char
  }
  return(chars);  //returns the array of chars
}
