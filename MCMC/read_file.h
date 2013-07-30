/*==========================================================================
 * Version: 5.00           date:xx/xx/xx                                    
 *                                                                          
 * Author: Francesco Montesano, MPE, Garching                               
 *                                                                          
 * Purpose: read the pk data to use in the mcmc chain
 *==========================================================================*/

#pragma once

#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>

#include "common.h"
#include "parse_ini.h"
#include "gsl_funcs.h"

class Dataset{
  private:
    ParseIni *dataset;

    struct nbins{
      int kitot, kimin, kimax, kiuse;  //bins in the measured power spectrum
      //bins in the kj file: used for the evaluation of the model power
      //spectrum before convolving with the window function
      int kjtot, kjmin, kjmax, kjuse;  
    } n_bins;
    //gsl vectors and matrices containing the content of files
    gsl_vector *data; // measured P(k)
    gsl_matrix *invcov; //inverse covariance
    gsl_matrix *Wij;  //window matrix
    gsl_vector *kj, *W0j, *G2i;     //kj, W0j, G^2(ki)
    double G20;  // G^2(ki=0), stored in the first element of the G20i file
    gsl_vector *convolved_model;  //model convolved with the window matrix
    gsl_vector *tempv; //store the matrix multiplication with the inverse covariance

    /*==========================================================================
     * read the number of bins to read from the files and to use for the
     * likelihood
     * Parameter
     * ---------
     *  ini: inifile
     *==========================================================================*/
    void get_bin_numbers(ParseIni ini);
    /*==========================================================================
     * read measured power
     * Parameter
     * ---------
     *  ini: inifile
     *==========================================================================*/
    void read_pk(ParseIni ini);
    /*==========================================================================
     * read and invert the covariance matrix
     * Parameter
     * ---------
     *  ini: inifile
     *==========================================================================*/
    void read_invert_cov(ParseIni ini);
    /*==========================================================================
     * read the files for the window matrix
     * Parameter
     * ---------
     *  ini: inifile
     *==========================================================================*/
    void read_window(ParseIni ini);

    /*==========================================================================
     * Allocate gls vectors and matrices used when convolving and computing the
     * chi^2
     *==========================================================================*/
    void alloc_gsl();

  public:
    /*==========================================================================
     * Constructor of the object
     * Parameter
     * ---------
     *  pk_file: string with the file name containing the name of the file to
     *  read in and the bins to consider
     *==========================================================================*/
    Dataset(std::string &pk_dataset);
    /*==========================================================================
     * Destructor of the object
     *==========================================================================*/
    ~Dataset();

    /*==========================================================================
     * return the gsl vector containing the wavenumber where to evaluate the
     * model power spectrum
     * output
     * ------
     *  k: gsl vector
     *==========================================================================*/
    gsl_vector *get_k(){return(kj);}

    /*==========================================================================
     * return the input gsl vector convolved with the window matrix
     * Parameters
     * ----------
     *  model: gsl vector
     *    power spectrum to be convolved
     * output
     * ------
     *  convolved_model: gsl vector
     *==========================================================================*/
    gsl_vector *convolve(gsl_vector *model);

    /*==========================================================================
     * compute the chi^2 give the model
     * Parameters
     * ----------
     *  model: gsl vector
     *    model power spectrum
     * output
     * ------
     *  chi2: double
     *    chi^2 from the data stored here and a model
     *==========================================================================*/
    double get_chi2(gsl_vector *model);

};  
