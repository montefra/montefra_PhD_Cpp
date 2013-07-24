/*==========================================================================
 * Version: 5.00           date:xx/xx/xx                                    
 *                                                                          
 * Author: Francesco Montesano, MPE, Garching                               
 *                                                                          
 * Purpose: read the pk data to use in the mcmc chain
 *==========================================================================*/

#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>

#include "parse_ini.h"

class Read_pk_files{
  private:
    ParseIni dataset;

    struct nbins{
      int kitot, kimin, kimax, kiuse;  //bins in the measured power spectrum
      //bins in the kj file: used for the evaluation of the model power
      //spectrum before convolving with the window function
      int kjtot, kjmin, kjmax, kjuse;  
    } n_bins;
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


  public:
    //gsl vectors and matrices containing the content of files
    gsl_vector *data; // measured P(k)
    gsl_vector *kj, *G2, *W0j;     //P(k), kj, G^2(ki), W0j
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
    }

};  
