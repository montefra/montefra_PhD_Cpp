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
 * read a line from the file, store position and velocity and return the mass
 * Parameters
 * ----------
 * inif: ifstream object with the file to read
 * output:
 * weight: w_fkp, w_fc+w_rf-1 and w_sys
 * pos[3]: position of the particle
 * n: number density
 * -------
 * return:
 * z: redshift
 *==========================================================================*/
double read_line(std::ifstream &inif, double *pos, double *w, double *n);

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
