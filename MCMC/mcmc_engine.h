/*==========================================================================
 * Version: 5.00           date:xx/xx/xx                                    
 *                                                                          
 * Author: Francesco Montesano, MPE, Garching                               
 *                                                                          
 * Purpose: initialise the parameter and run the mcmc 
 *==========================================================================*/

#pragma once

#include <algorithm>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <iostream>
#include <map>
#include <string>
#include <vector>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include "common.h"
#include "likelihood.h"
#include "parse_ini.h"

class MCMC{
  private:
    std::vector<Likelihood> likes; //list of likelihoods
    size_t n_likes;

    std::string file_root;  //root of the ouput files
    size_t n_steps;  //number of steps in the chain
    gsl_rng *r, *rg;   // gsl random generators
    //parameter names from the likelihoods
    std::vector<std::string> paramnames, long_names;
    size_t n_parameters;  //number of parameters

    //maps containing the starting point, the mimimum and the maximum values
    //possible and the step size
    std::map<std::string, double> start, min, max, step;
    /*==========================================================================
     * initialiase the random generators
     *==========================================================================*/
    void init_random(){
      r = gsl_rng_alloc(gsl_rng_taus2);
      rg = gsl_rng_alloc(gsl_rng_taus2);
      gsl_rng_set(r, time(NULL));
      gsl_rng_set(rg, time(NULL)+1);
    }

    /*==========================================================================
     * retrieve the parameter names from the likelihoods
     *==========================================================================*/
    void get_paramnames();
    /*==========================================================================
     * save them in the paramnames file
     *==========================================================================*/
    void save_paramnames();
    /*==========================================================================
     * read the paramnames limits and steps from the inifile and save in 
     * maps
     * Parameters
     * ----------
     *  ini: ParseIni object
     *    inifile containing the parameters 
     *==========================================================================*/
    void read_params(ParseIni &ini);

    /*==========================================================================
     * initialise the parameter names (mostly call the three above functions)
     *==========================================================================*/
    void initialise_paramnames(ParseIni &ini){
      get_paramnames();
      save_paramnames();
      long_names.clear();  //clear the long names
      read_params(ini);  //get the paramnames limits and steps from the inifile
    }

    /*==========================================================================
     * randomly generate a new set of parameters
     * Parameters
     * ----------
     *  params: map<string, double>
     *    map containing the parameters to substitute
     * output
     * ------
     *  newparams: map<string, double>   
     *    map with the new parameters
     *==========================================================================*/
    std::map<std::string, double> new_parameters(std::map<std::string, double>
        &params);

  public:
    /*==========================================================================
     * Constructor
     * get the number of datasets used and the inifile objects containing the
     * parameter limits
     * Parameters
     * ----------
     *  n_datasets: int
     *    number of datasets
     *  ini: ParseIni object
     *    ini object containing the the parametes limits and step
     *==========================================================================*/
    MCMC(std::vector<Likelihood> &likes, ParseIni &ini);

    /*==========================================================================
     * Destructor
     *==========================================================================*/
    ~MCMC(){};

    /*==========================================================================
     * run the mcmc chain
     *==========================================================================*/
    void run();

};
