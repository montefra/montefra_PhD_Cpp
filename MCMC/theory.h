/*==========================================================================
 * Version: 5.00           date:xx/xx/xx                                    
 *                                                                          
 * Author: Francesco Montesano, MPE, Garching                               
 *                                                                          
 * Purpose: read the linear and 1loop power spectrum and save them into 
 * gsl spline objects. 
 * Provide method(s) to that return the model
 *==========================================================================*/

#pragma once

#include <cmath>
#include <map>
#include <string>
#include <vector>

#include <gsl/gsl_spline.h>

#include "parse_ini.h"
#include "gsl_funcs.h"


class Theory{
  private:

    gsl_interp_accel *acclin, *acc1l;  //interpolation accelerators
    gsl_spline *spllin, *spl1l; //spline objects

    //vectors with teh short and long name of the parameters used in the model
    //this are passed to the likelihood object that creates the proper
    //paramnames attaching the number of the dataset and then passed to the
    //mcmc engine that print out the short and long name and read the parameter
    //limits
    std::vector<std::string> paramnames, long_names;
    /*==========================================================================
     * fill the paramnames with the names of the parameters used here
     * the order in the paramnames vector is the same as the one accepted by
     * 'get_model'
     *==========================================================================*/
    void fill_paramnames();

  public:

    /*==========================================================================
     * Constructor
     *
     * Parameters
     * ----------
     *  ini: ParseIni object
     *    ini parser objects containing. This code need that ini contains 'plin'
     *    and 'p1loop' files
     *==========================================================================*/
    Theory(ParseIni ini);
    /*==========================================================================
     * Destructor
     * deallocate stuff
     *==========================================================================*/
    ~Theory();

    /*==========================================================================
     * return the long and short names of the variables
     *==========================================================================*/
    std::vector<std::string> get_paramnames() {return(paramnames);}
    std::vector<std::string> get_long_names() {return(long_names);}

    /*==========================================================================
     * compute the model in wavenumber k given the model parameter in double
     * array. The parameter names are defined in 'fill_paramnames()'
     * Parameters
     * ----------
     *  k: double
     *    wavenumber where to evaluate the model
     *  params: map
     *    map containing the parameter name as key and its value as value
     * output
     * ------
     *  pk: double
     *    model power spectrum evaluated in k
     *==========================================================================*/
    double get_model(double k, std::map<std::string, double> params);

};

