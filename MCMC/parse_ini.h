/*==========================================================================
 * Version: xxxx           date:xxxxxxxx                                    
 *                                                                          
 * Author: Francesco Montesano, MPE, Garching                               
 *                                                                          
 * Purpose: parse the inifile and provide methods to access its content
 *                                                                          
 *==========================================================================*/

#pragma once

#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "common.h"

class ParseIni
{

  private:
    const std::string ininame;  //name of the ini file 
    std::vector<std::string> inis;   //vector containing the lines of the inifile
    size_t inisize;   //size of the inifile vector
    std::string comment;   //string preceding comments

    /*==========================================================================
     * read file 'inifile', skipping empty lines,
     * and stripping comments
     * Parameters
     * ----------
     * fname: string
     *   ini file name
     * output
     * ------
     * lines: vector of strings
     *   inifile stripped all all the empty lines and comments
     *==========================================================================*/
    std::vector<std::string> readlines(std::string fname);

  public:

    /*==========================================================================
     * Construnctor
     * Parameters
     * ----------
     * inifile: string
     *   name of the input file contining the input parameters
     *==========================================================================*/
    ParseIni(std::string inifile);
    /*==========================================================================
     * Destructor
     *==========================================================================*/
    ~ParseIni(){this->inis.clear();}

    /*==========================================================================
     * return the name of the ini file
     * output
     * ------
     *  fname: string
     *==========================================================================*/
    std::string get_fname(){return(ininame);}

    /*==========================================================================
     * standard error when retreaving a value from a inifile
     * parameters
     * ----------
     *  param: string
     *    name of the searched parameter
     *  error: int
     *    error code to give to exit
     *==========================================================================*/
    void ini_error(std::string param, int error);


    /*==========================================================================
     * Read MCMC paramter starting guess, boundaries and sigma
     * Parameters
     * ----------
     * parname: string
     *   name of the parameter in the inifile
     * output
     * ------
     * guess: double
     *   guess of the initial position
     * lower: double
     *   lower limit of the flat prior
     * upper: double
     *   upper limit of the flat prior
     * sigma: double
     *   standard deviation of the gaussian used for the random walk
     * return
     * ------
     * error: int
     *   0: if no error occurs; 1: if the parameter do not exists, 
     *   -1: if less than 4 arguments are found
     *==========================================================================*/
    int get_mcmc_params(std::string parname, double *guess, double* lower,
        double *upper, double *sigma);

    /*==========================================================================
     * Read a single parameter following a given string and '='
     * Parameters
     * ----------
     * parname: string
     *   parameter to search
     * output
     * ------
     * value: std::string, int, double, bool
     *   first value following the parameter name and '='
     *   if bool:
     *     if the parameter after '=' is T,t,True,true retruns true
     *     if the parameter after '=' is F,f,False,false retruns false
     * return
     * ------
     * error: int
     *   0: if no error occurs; 1: if the parameter do not exists, 
     *   -1: if the parameter has no argument
     *   if bool:
     *     -2 if the parameter is not among the choises of above
     *   otherwise:
     *     -2 the values could not be converted to the desired type
     *==========================================================================*/
    int get_param(std::string parname, std::string *value);
    int get_param(std::string parname, int *value);
    int get_param(std::string parname, size_t *value);
    int get_param(std::string parname, double *value);
    int get_param(std::string parname, bool *value);
    /*==========================================================================
     * Read n parameters and return them in an vector
     * Parameters
     * ----------
     * parname: string
     *   parameter to search
     * n: int
     *   number of parameters to read
     * output
     * ------
     * values: std vector of strint, int, double
     *   first n values following the parameter name and '='
     * return
     * ------
     * error: int
     *   0: if no error occurs; 1: if the parameter do not exists, 
     *   -1: if the parameter has less than n arguments
     *   -2: the values could not be converted to the desired type
     *==========================================================================*/
    int get_param_list(std::string parname, size_t n, std::vector<std::string> &values);
    int get_param_list(std::string parname, size_t n, std::vector<int> &values);
    int get_param_list(std::string parname, size_t n, std::vector<size_t> &values);
    int get_param_list(std::string parname, size_t n, std::vector<double> &values);
    //int get_param_list(std::string parname, size_t n, std::vector<bool> &values);
};
/*==========================================================================
 * END OF THE CLASS
 *==========================================================================*/

