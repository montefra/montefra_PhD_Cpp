/*==========================================================================*/
/* Author: Francesco Montesano, montefra@mpe.mpg.de                         */
/* 06/02/2012,    version: 0.01                                             */
/*                                                                          */
/* compute the averaged power spectrum from a simulation box                */
/*==========================================================================*/
#ifndef FKPPOWER_SPECTRUM_H
#define FKPPOWER_SPECTRUM_H

#include <algorithm>
#include <fstream>
#include <string>
#include <vector>

//TCLAP command line interface package
#include <tclap/CmdLine.h> 
//class to compute the real to complex and complex to real FFT
#include "ps_r2c_c2r.h" 

#include "common_functions.h"

/*=======================================================================*/
/* get the description and the version of the code                       */
/*=======================================================================*/
std::string get_description();
std::string get_version();

/*=======================================================================
 * Print error, finalise mpi and exit
 * Parameters
 * ----------
 *  message: error message to print to std err
 *  code: error code
 *  myrank: rank of the processor
 *  root: root processor
 *  com: MPI communicator
 *=======================================================================*/
void on_error(std:string message, int error_code, int myrank, int, root, MPI_Comm
    com){
  if(myrank == root)
    std::cerr << message << std::endl;
  MPI_Barrier(com);
  MPI_Finalize();
  exit(error_code);
}

/*=======================================================================*/
/* check input files                                                     */
/*=======================================================================*/
/*==========================================================================
 * Check the number of files 
 * Parameters
 * ----------
 * files: vector of strings with the input files
 * winonly: if the window function required
 * cross: if the cross power spectrum or window function required
 * myrank: rank of the processor
 * root: root processor
 * com: MPI communicator
 *==========================================================================*/
void check_input_files(std::vector<std::string> &files, bool winonly, bool
    cross, int myrank, int root, MPI_Comm com);

/*==========================================================================
 * compute alpha, square of the inverse normalisation and 
 * shot noise amplitude given the sums over the weights
 *==========================================================================*/
/*==========================================================================
 * for the window function
 * Parameters
 * ----------
 *  sums: array containing sum(w), sum(n*w^2), sum(w^2)
 *  alpha: in this case 1 (for complementarity with the power spectrum case)
 *  N2: square of the inverse normalisation
 *  noise: shot noise amplitude
 *==========================================================================*/
void alpha_N_sh(double *sums, double *alpha, double *N2, double *noise){
  *alpha=1.;  //alpha is one in this case
  *N2 = 1./sums[1];
  *noise = sums[2]*N2;
}
/*==========================================================================
 * for the power spectrum
 * Parameters
 * ----------
 *  sumsdat: array containing sum(w), sum(n*w^2), sum(w^2) from the data
 *  sumsran: array containing sum(w), sum(n*w^2), sum(w^2) from the random
 *  alpha: in this case 1 (for complementarity with the power spectrum case)
 *  N2: square of the inverse normalisation
 *  noise: shot noise amplitude
 *==========================================================================*/
void alpha_N_sh(double *sumsdat, double *sumsran, double *alpha, double *N2,
    double *noise){
  *alpha = sumscat[0]/sumsran[0];
  *N2 = 1./(*alpha*sumsran[1]);
  *noise = (*alpha+1.) * *alpha*sumsran[2] * *N2;
}

/*==========================================================================
 * class created to read all the files. Provides a simple interface,
 * a constructor that allocates arrays used while reading the files
 *   and saves all the common information to avoid passing them each time
 * a function that read the files and save the content into the desired grid
 *==========================================================================*/
class readfiles{

  public:
    /*==========================================================================
     * constructor of the class
     * Parameters
     * ----------
     * cell_size: side of a cell
     * min_coordinate: minimum coordinate of the box used to convert the
     *   the position of each particle in grid units
     * pw: value of P(k) to be used to compute the FKP weights
     *==========================================================================*/
    readfiles(double cell_size, double min_coordinate, double pw){
      fl.pos = new double[3];  //array containing the position
      fl.w = new double[3];  //weight: w_fkp, w_fc+w_rf-1 and w_sys
      this->cell_size = cell_size;
      this->min_coordinate = min_coordinate;
      this->pw = pw;
    }
    ~readfiles(){
      delete [] fl.pos;   //clear velocity and position
      delete [] fl.w;
    }

    /*==========================================================================
     * Read the input ascii file and fill the grid
     * Parameters
     * ----------
     * ifile: input file name
     * grid: 'ps_r2c_c2r_mpi_inplace' object
     * zrange: array of two float: only objects within the redshift range considered
     * zrange: range of redshift
     * ignorew: vector<int> (-1,1,2): if !=-1 set to 1 fl.w[1] or fl.w[2]
     * repeatw: bool: if true, assign object fl.w[1] times with fl.w[1]=1
     * output
     * ------
     * sums: array containing sum(w), sum(n*w^2), sum(w^2)
     *==========================================================================*/
    double *read_file(std::string ifile, ps_r2c_c2r_mpi_inplace &grid,
        std::vector<double> zrange, std::vector<int> ignorew, bool repeatw);

  private: 

    struct fileline{
      double *pos, *w; // positions and weight from the input file
      double n, z;  // number density and redshift
    } fl;
    //size of the cell, minimum coordinate used to transform the position in
    //grid coordinated and FKP weight
    double cell_size, min_coordinate, pw; 

    /*==========================================================================
     * allocate and initialise to 0 the array that contains the sums
     * sum(w), sum(n*w^2), sum(w^2)
     *==========================================================================*/
    double *init_sums(){
      double *sums = new double[3]; 
      for(int i=0; i<3; ++i) sums[i] =0.; 
      return(sums);
    }

    /*==========================================================================
     * read in the file with the various options. 
     * For argument refer to the description of 'read_file'
     *==========================================================================*/
    void read_asitis(std::ifstream &inif, ps_r2c_c2r_mpi_inplace &grid, double* sums);
    void read_zrange(std::ifstream &inif, ps_r2c_c2r_mpi_inplace &grid, double* sums,
        std::vector<double> zrange);
    void read_ignorew(std::ifstream &inif, ps_r2c_c2r_mpi_inplace &grid, double* sums,
        std::vector<int> ignorew);
    void read_repeatw(std::ifstream &inif, ps_r2c_c2r_mpi_inplace &grid, double* sums);
    void read_zrange_ignorew(std::ifstream &inif, ps_r2c_c2r_mpi_inplace &grid,
        double* sums, std::vector<double> zrange, std::vector<int> ignorew);
    void read_zrange_repeatw(std::ifstream &inif, ps_r2c_c2r_mpi_inplace &grid,
        double* sums, std::vector<double> zrange);
    void read_ignorew_repeatw(std::ifstream &inif, ps_r2c_c2r_mpi_inplace &grid,
        double* sums, std::vector<int> ignorew);
    void read_zrange_ignorew_repeatw(std::ifstream &inif, ps_r2c_c2r_mpi_inplace
        &grid, double* sums, std::vector<double> zrange, std::vector<int> ignorew);

    /*==========================================================================
     * read a line from the file, store position, weights, number density and
     * redshift into 'fileline' struct. 
     * Convert position to grid units and compute the FKP weight
     * Parameters
     * ----------
     * inif: ifstream object with the file to read
     *==========================================================================*/
    void read_line(std::ifstream &inif){
      double dummy;  //dummy variable
      //read a line
      inif >> fl.pos[0] >> fl.pos[1] >> fl.pos[2] >> fl.w[0] >> fl.w[1] >>
        fl.w[2] >> fl.n >> fl.z >> dummy;
      to_gridunits();
      fl.w[0] /= (1.+pw*fl.n);   //multiply the intrinsic weight by the fkp weight
    }
    /*==========================================================================
     * convert the position in grid units
     *==========================================================================*/
    void to_gridunits(){
        for(int i=0; i<3; ++i)
          fl.pos[i] = (fl.pos[i]-min_coordinate)/cell_size;
    }
    /*==========================================================================
     * add to the sum of the weights the weights of the current object
     *==========================================================================*/
    void sums_weights(double *sums){
      sums[0] += fl.w[0]*fl.w[1]*fl.w[2];   //sum(w) 
      sums[1] += fl.w[0]*fl.w[0]*fl.w[1]*fl.w[1]*fl.w[2]*fl.w[2] * fl.n;  //sum( n*w^2 )
      //sums[1] += fl.w[0]*fl.w[0]*fl.w[1]*fl.w[2]*n;  //sum( n*w^2 )
      sums[2] += fl.w[0]*fl.w[0]*fl.w[1]*fl.w[1]*fl.w[2]*fl.w[2];    //sum(w^2)
    }
};

#endif
