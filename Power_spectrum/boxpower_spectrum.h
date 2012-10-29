/*==========================================================================*/
/* Author: Francesco Montesano, montefra@mpe.mpg.de                         */
/* 06/02/2012,    version: 0.01                                             */
/*                                                                          */
/* compute the averaged power spectrum from a simulation box                */
/*==========================================================================*/

#ifndef BOXPOWER_SPECTRUM_H
#define BOXPOWER_SPECTRUM_H

#include <fstream>
#include <string>
#include <vector>

//TCLAP command line interface package
#include <tclap/CmdLine.h> 

#include "common_functions.h"
#include "read_snapshot.h"
//class to compute the real to complex and complex to real FFT
#include "ps_r2c_c2r.h" 

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
 * mass: mass of the particle
 * pos[3], vel[3]: position and velocity of the particle
 *==========================================================================*/
double read_line(std::ifstream &inif, double *pos, double *vel);

/*==========================================================================
 * Read the input ascii file containing haloes and fill the grid
 * no mass limit, no redshift space distortions
 * Parameters
 * ----------
 * ifile: input file name
 * grid: 'ps_r2c_c2r_mpi_inplace' object
 * cell_size: side of a cell
 * output
 * ------
 * n_haloes: number of haloes assigned to the grid
 *==========================================================================*/
ptrdiff_t read_haloes(std::string ifile, ps_r2c_c2r_mpi_inplace &grid, double cell_size);
/*==========================================================================
 * Read the input ascii file containing haloes and fill the grid
 * lower mass limit, no redshift space distortions
 * Parameters
 * ----------
 * ifile: input file name
 * grid: 'ps_r2c_c2r_mpi_inplace' object
 * mmin: lower limit for the mass
 * cell_size: side of a cell
 * output
 * ------
 * n_haloes: number of haloes assigned to the grid
 *==========================================================================*/
ptrdiff_t read_haloes(std::string ifile, ps_r2c_c2r_mpi_inplace &grid, double mmin, 
    double cell_size);
/*==========================================================================
 * Read the input ascii file containing haloes and fill the grid
 * mass bin, no redshift space distortions
 * Parameters
 * ----------
 * ifile: input file name
 * grid: 'ps_r2c_c2r_mpi_inplace' object
 * mmin, mmax: lower and upper limit for the mass
 * cell_size: side of a cell
 * output
 * ------
 * n_haloes: number of haloes assigned to the grid
 *==========================================================================*/
ptrdiff_t read_haloes(std::string ifile, ps_r2c_c2r_mpi_inplace &grid, double mmin, 
    double mmax, double cell_size);

/*==========================================================================
 * Read the input ascii file containing haloes and fill the grid
 * no mass limit, redshift space distortions
 * Parameters
 * ----------
 * ifile: input file name
 * grid: 'ps_r2c_c2r_mpi_inplace' object
 * cell_size: side of a cell
 * zdist: 0: x axis; 1: y axis; 2: z axis
 * ncells: number of cells in the direction of the redshif space 
 *   distortions. Used to implement periodic boundary conditions
 * output
 * ------
 * n_haloes: number of haloes assigned to the grid
 *==========================================================================*/
ptrdiff_t read_haloes(std::string ifile, ps_r2c_c2r_mpi_inplace &grid, double cell_size,
    int zdist, ptrdiff_t ncells);
/*==========================================================================
 * Read the input ascii file containing haloes and fill the grid
 * lower mass limit, redshift space distortions
 * Parameters
 * ----------
 * ifile: input file name
 * grid: 'ps_r2c_c2r_mpi_inplace' object
 * mmin: lower limit for the mass
 * cell_size: side of a cell
 * zdist: 0: x axis; 1: y axis; 2: z axis
 * ncells: number of cells in the direction of the redshif space 
 *   distortions. Used to implement periodic boundary conditions
 * output
 * ------
 * n_haloes: number of haloes assigned to the grid
 *==========================================================================*/
ptrdiff_t read_haloes(std::string ifile, ps_r2c_c2r_mpi_inplace &grid, double mmin, 
    double cell_size, int zdist, ptrdiff_t ncells);
/*==========================================================================
 * Read the input ascii file containing haloes and fill the grid
 * mass bin, redshift space distortions
 * Parameters
 * ----------
 * ifile: input file name
 * grid: 'ps_r2c_c2r_mpi_inplace' object
 * mmin, mmax: lower and upper limit for the mass
 * cell_size: side of a cell
 * zdist: 0: x axis; 1: y axis; 2: z axis
 * ncells: number of cells in the direction of the redshif space 
 *   distortions. Used to implement periodic boundary conditions
 * output
 * ------
 * n_haloes: number of haloes assigned to the grid
 *==========================================================================*/
ptrdiff_t read_haloes(std::string ifile, ps_r2c_c2r_mpi_inplace &grid, double mmin, 
    double mmax, double cell_size, int zdist, ptrdiff_t ncells);

/*==========================================================================
 * Read the input ascii file containing haloes and fill the grids
 * cross power, no redshift space distortions
 * Parameters
 * ----------
 * ifile: input file name
 * grid1, grid2: 'ps_r2c_c2r_mpi_inplace' objects
 * mass: 4 element vector containing the mass limits
 * cell_size: side of a cell
 * output
 * ------
 * n_haloes: number of haloes assigned to the grid
 *==========================================================================*/
ptrdiff_t *read_haloes(std::string ifile, ps_r2c_c2r_mpi_inplace &grid1, ps_r2c_c2r_mpi_inplace &grid2, 
    std::vector<double> mass, double cell_size);
/*==========================================================================
 * Read the input ascii file containing haloes and fill the grids
 * cross power, no redshift space distortions
 * Parameters
 * ----------
 * ifile: input file name
 * grid1, grid2: 'ps_r2c_c2r_mpi_inplace' objects
 * mass: 4 element vector containing the mass limits
 * cell_size: side of a cell
 * zdist: 0: x axis; 1: y axis; 2: z axis
 * ncells: number of cells in the direction of the redshif space 
 *   distortions. Used to implement periodic boundary conditions
 * output
 * ------
 * n_haloes: number of haloes assigned to the grid
 *==========================================================================*/
ptrdiff_t *read_haloes(std::string ifile, ps_r2c_c2r_mpi_inplace &grid1, ps_r2c_c2r_mpi_inplace &grid2, 
    std::vector<double> mass, double cell_size, int zdist, ptrdiff_t ncells);

/*=======================================================================*/
/* read dark matter catalogues from binary file                          */
/*=======================================================================*/
/*==========================================================================
 * Read the input binary file containing dark matter and fill the grid
 * no redshift space distortions. Gadget 2 binary file
 * Parameters
 * ----------
 * ifile: input file root
 * n_files: number of files
 * grid: 'ps_r2c_c2r_mpi_inplace' object
 * cell_size: side of a cell
 * output
 * ------
 * n_dm: number of dark matter particles assigned to the grid
 *==========================================================================*/
ptrdiff_t read_dm(std::string ifile, int n_files, ps_r2c_c2r_mpi_inplace &grid, 
    double cell_size);
/*==========================================================================
 * Read the input binary file containing dark matter and fill the grid
 * redshift space distortions. Gadget 2 binary file
 * Parameters
 * ----------
 * ifile: input file root
 * n_files: number of files
 * grid: 'ps_r2c_c2r_mpi_inplace' object
 * cell_size: side of a cell
 * zdist: 0: x axis; 1: y axis; 2: z axis
 * ncells: number of cells in the direction of the redshif space 
 *   distortions. Used to implement periodic boundary conditions
 * output
 * ------
 * n_dm: number of dark matter particles assigned to the grid
 *==========================================================================*/
ptrdiff_t read_dm(std::string ifile, int n_files, ps_r2c_c2r_mpi_inplace &grid, 
    double cell_size, int zdist, ptrdiff_t ncells);

#endif
