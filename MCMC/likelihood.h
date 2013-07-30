/*==========================================================================
 * Version: 5.00           date:xx/xx/xx                                    
 *                                                                          
 * Author: Francesco Montesano, MPE, Garching                               
 *                                                                          
 * Purpose: likelihood class. Basically contains a data object and a theory
 * object and take care of connecting those two to the mcmc engine
 *==========================================================================*/

#pragma once

#include <cmath>
#include <map>
#include <string>
#include <vector>

#include <gsl/gsl_vector.h>

#include "common.h"
#include "parse_ini.h"
#include "read_file.h"
#include "theory.h"

class Likelihood{

  private:
    int dataset;    //dataset number
    ParseIni *ini;   //inifile
    Theory *theory;  //theory objec
    Dataset *data;

    //short and long parameter names to pass to the mcmc chains
    std::vector<std::string> paramnames, long_names;
    //maps the paramnames in the theory object to the ones passed to the mcmc engine
    //the second are the first plus the dataset number appended
    std::map<std::string, std::string> mcmc2theory;

    //gsl vectors containing the value of k for the model power spectrum
    gsl_vector *k;
    size_t n_k; //size of k
    // theory power spectrum and convolved power spectrum
    gsl_vector *pktheory, *pkconvolved; 

    /*==========================================================================
     * get the paramnames from the theory object, append the dataset number 
     * to the paramnames and save them
     *==========================================================================*/
    void retrieve_paramnames();

  public:
    /*==========================================================================
     * Constructor
     * Parameters
     * ----------
     *  dataset: int
     *    dataset number
     *  ini: ParseIni object
     *    container of the inifile
     *  theory: Theory object
     *    used to get the model
     *==========================================================================*/
    Likelihood(int dataset, ParseIni &ini, Theory &theory);

    /*==========================================================================
     * return the long and short names of the variables
     *==========================================================================*/
    std::vector<std::string> get_paramnames() {return(paramnames);}
    std::vector<std::string> get_long_names() {return(long_names);}

    /*==========================================================================
     * return the likelihood with the given parameters
     * Parameters
     * ----------
     *  params: map
     *    map with parameter name as key and its value as value
     *==========================================================================*/
    double get_like(std::map<std::string, double> params);

};
