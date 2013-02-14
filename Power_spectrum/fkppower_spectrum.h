/*==========================================================================*/
/* Author: Francesco Montesano, montefra@mpe.mpg.de                         */
/* 06/02/2012,    version: 0.01                                             */
/*                                                                          */
/* compute the averaged power spectrum from a simulation box                */
/*==========================================================================*/
#ifndef FKPPOWER_SPECTRUM_H
#define FKPPOWER_SPECTRUM_H

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

/*=======================================================================*/
/* check input files                                                     */
/*=======================================================================*/
/*==========================================================================
 * Do the checks when the cross power spectrum is required
 * if files.size > 4, files is resized to 4
 * Parameters
 * ----------
 * files: vector of strings with the input files
 * myrank: rank of the processor
 * root: root processor
 * com: MPI communicator
 *==========================================================================*/
void check_crosspk(std::vector<std::string> &files, int myrank, int root,
    MPI_Comm com);


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
     * ignorew: int (-1,1,2): if !=-1 set to 1 fl.w[1] or fl.w[2]
     * repeatw: bool: if true, assign object fl.w[1] times with fl.w[1]=1
     * output
     * ------
     * sums: array containing sum(w), sum(n*w^2), sum(w^2)
     *==========================================================================*/
    double *read_file(std::string ifile, ps_r2c_c2r_mpi_inplace &grid,
        std::vector<double> zrange, int ignorew, bool repeatw);

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
    void read_asitis(std::string ifile, ps_r2c_c2r_mpi_inplace &grid, double* sums);
    void read_zrange(std::string ifile, ps_r2c_c2r_mpi_inplace &grid, double* sums,
        std::vector<double> zrange);
    void read_ignorew(std::string ifile, ps_r2c_c2r_mpi_inplace &grid, double* sums,
        int ignorew);
    void read_repeatw(std::string ifile, ps_r2c_c2r_mpi_inplace &grid, double* sums);
    void read_zrange_ignorew(std::string ifile, ps_r2c_c2r_mpi_inplace &grid,
        double* sums, std::vector<double> zrange, int ignorew);
    void read_zrange_repeatw(std::string ifile, ps_r2c_c2r_mpi_inplace &grid,
        double* sums, std::vector<double> zrange);
    void read_ignorew_repeatw(std::string ifile, ps_r2c_c2r_mpi_inplace &grid,
        double* sums, int ignorew);
    void read_zrange_ignorew_repeatw(std::string ifile, ps_r2c_c2r_mpi_inplace
        &grid, double* sums, std::vector<double> zrange, int ignorew);

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
