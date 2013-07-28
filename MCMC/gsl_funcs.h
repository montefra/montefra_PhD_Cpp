/*==========================================================================
 * Version: 5.00           date:xx/xx/xx                                    
 *                                                                          
 * Author: Francesco Montesano, MPE, Garching                               
 *                                                                          
 * Purpose: container for functions using gsl
 *==========================================================================*/

#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_spline.h>

/*==========================================================================
 * read a vector of size dim from file 'file_name'
 * Parameter
 * ---------
 *  file_name: string
 *    file name of the vector
 *  dim: int
 *    dimention of the input vector
 * output
 * ------
 *  vector: gsl vector
 *==========================================================================*/
gsl_vector *read_gsl_vector(std::string file_name, int dim);
/*==========================================================================
 * read a vector of size dim from file 'file_name', cut 'cut_dim' using imin as
 * the minimum index
 * Parameter
 * ---------
 *  file_name: string
 *    file name of the vector
 *  dim: int
 *    dimention of the input vector
 *  cut_dim: int
 *    dimension of the cut vector
 *  imin: int
 *    minimum index in the cut vector
 * output
 * ------
 *  vector: gsl vector
 *==========================================================================*/
gsl_vector *read_cut_gsl_vector(std::string file_name, int dim, int cut_dim,
    int imin);

/*==========================================================================
 * read a matrix of size xdim*ydim from file 'file_name'
 * Parameter
 * ---------
 *  file_name: string
 *    file name of the matrix
 *  xdim, ydim: int
 *    dimentions of the input matrix
 * output
 * ------
 *  matrix: gsl matrix
 *==========================================================================*/
gsl_matrix *read_gsl_matrix(std::string file_name, int xdim, int ydim);
/*==========================================================================
 * read a matrix of size xdim*ydim from file 'file_name', cut to 
 * 'cut_xdim*cut_ydim' using xmin and ymin as the mimimum x and y indeces
 * Parameter
 * ---------
 *  file_name: string
 *    file name of the matrix
 *  xdim, ydim: int
 *    dimentions of the input matrix
 *  cut_xdim, cut_ydim: int
 *    dimension of the cut matrix
 *  xmin, ymin: int
 *    minimum x and y index in the cut matrix
 * output
 * ------
 *  matrix: gsl matrix
 *==========================================================================*/
gsl_matrix *read_cut_gsl_matrix(std::string file_name, int xdim, int ydim, int
    cut_xdim, int cut_ydim, int xmin, int ymin);


/*==========================================================================
 * read file fname (made of two columns) and save it into a spline object
 * Parameters
 * ----------
 *  fname: string
 *    name of the file to read
 * output
 * ------
 *  spl: gsl spline object
 *    containing the read file
 *==========================================================================*/
gsl_spline *read_2_spline(std::string fname);

