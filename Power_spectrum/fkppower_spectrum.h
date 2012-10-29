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

/*==========================================================================
 * read a line from the file, store position and velocity and return the mass
 * Parameters
 * ----------
 * inif: ifstream object with the file to read
 * output:
 * weight: intrinsic and systematic weights contained in the input file
 * pos[3]: position of the particle
 * n: number density
 * z: redshift
 *==========================================================================*/
double* read_line(std::ifstream &inif, double *pos, double *n, double *z);


/*==========================================================================
 * Read the input ascii file and fill the grid
 * Parameters
 * ----------
 * ifile: input file name
 * grid: 'ps_r2c_c2r_mpi_inplace' object
 * cell_size: side of a cell
 * min_coordinate: minimum coordinate of the box used to convert the
 *   the position of each particle in grid units
 * pw: value of P(k) to be used to compute the FKP weights
 * output
 * ------
 * sums: array containing sum(w), sum(n*w^2), sum(w^2)
 *==========================================================================*/
double *read_file(std::string ifile, ps_r2c_c2r_mpi_inplace &grid, 
    double cell_size, double min_coordinate, double pw);

#endif
