/*==========================================================================
 * Version: 5.00           date:xx/xx/xx                                    
 *                                                                          
 * Author: Francesco Montesano, MPE, Garching                               
 *                                                                          
 * Purpose: Markov Chain Monte Carlo
 *==========================================================================*/

#include<iostream>
#include<string>
#include<vector>

#include "parse_ini.h"
#include "mcmc_like.h"

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
string to_string(T const& value) {
    stringstream sstr;
    sstr << value;
    return sstr.str();
}
