/*==========================================================================
 * Version: 5.00           date:xx/xx/xx                                    
 *                                                                          
 * Author: Francesco Montesano, MPE, Garching                               
 *                                                                          
 * Purpose: read the pk data to use in the mcmc chain
 *==========================================================================*/

#include <cstdlib>
#include <iostream>
#include <sstream>
#include <string>

#include <gsl/gls_linalg.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>

#include "parse_ini.h"
#include "gsl_funcs.h"

class Read_pk_files{
  private:
    ParseIni dataset;

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
    void alloc_gsl(){

    }

  public:
    /*==========================================================================
     * Constructor of the object
     * Parameter
     * ---------
     *  pk_file: string with the file name containing the name of the file to
     *  read in and the bins to consider
     *==========================================================================*/
    Read_pk_files(std::string pk_dataset);
    /*==========================================================================
     * Destructor of the object
     *==========================================================================*/
    ~Read_pk_files(std::string pk_dataset){
      gsl_vector_free(data);
      gsl_matrix_free(invcov);
      gsl_matrix_free(Wij);
      gsl_vector_free(kj);
      gsl_vector_free(W0j);
      gsl_vector_free(G2i);
    }

    /*==========================================================================
     * return the gsl vector containing the wavenumber where to evaluate the
     * model power spectrum
     * output
     * ------
     *  k: gsl vector
     *==========================================================================*/
    gsl_vector *get_k(){return{kj}};

};  
