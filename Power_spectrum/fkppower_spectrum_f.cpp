/*==========================================================================*/
/* Author: Francesco Montesano, montefra@mpe.mpg.de                         */
/* 06/02/2012,    version: 0.01                                             */
/*                                                                          */
/* compute the averaged power spectrum from a simulation box                */
/* Definition of the functions                                              */
/*==========================================================================*/

#include "fkppower_spectrum.h"

/*=======================================================================*/
/* Code specific functions definitions                                   */
/*=======================================================================*/
std::string get_description()
/* description of the code   */
{
  std::string description = "This code computes the angular averaged power spectrum ";
  description += "from catalogues with a realistic geometry using the FKP estimator.\n";
  description += "The input file must have the following structure:\n";
  description += "\tx y z w_fkp (w_fc+w_rf-1) w_systematics n(z) z column\n";
  description += "(the last column is ignored)";
  return(description);
}
std::string get_version()
/*version of the code */
{ return("0.01"); }

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
double read_line(std::ifstream &inif, double *pos, double *w, double *n){
  double z, dummy;  //dummy variable
  inif >> pos[0] >> pos[1] >> pos[2] >> w[0] >> w[1] >> w[2] >> *n >> z >> dummy;   //read a line
  return(z);
}

/*=======================================================================*/
/* read the galaxy or random catalogues from an ascii file               */
/*=======================================================================*/
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
    double cell_size, double min_coordinate, double pw){
  std::ifstream in(ifile.c_str());  //open the input file 

  double *pos=new double[3];  //array containing the position
  double *w = new double[2];  //weight: w_fkp, w_fc+w_rf-1 and w_sys
  double n, z;   //number density and redhift
  double *sums = new double[3]; //sum(w), sum(n*w^2), sum(w^2)
  for(int i=0; i<3; ++i) sums[i] =0.; //initialise the sums

  for(;;){  //loop through the lines of the file
    z = read_line(in, pos, w, &n);
    if (in.eof() == true) break;    //breack if the end of the file reached
    for(int i=0; i<3; ++i) 
      pos[i] = (pos[i]-min_coordinate)/cell_size;   //convert the position in grid units

    w[0] /= (1.+pw*n);   //multiply the intrinsic weight by the fkp weight
    sums[0] += w[0]*w[1]*w[2];   //sum(w) 
    sums[1] += w[0]*w[0]*w[1]*w[1]*w[2]*w[2] * n;  //sum( n*w^2 )
    sums[2] += w[0]*w[0]*w[1]*w[1]*w[2]*w[2];    //sum(w^2)

    grid.assign_particle(pos[0], pos[1], pos[2], w[0]*w[1]*w[2]);
    //grid.assign_particle_periodic(pos[0], pos[1], pos[2], w[0]*w[1]);
  }
  in.close();   //close the input file
  delete [] pos;   //clear velocity and position
  delete [] w;
  return(sums);  //return the number of haloes assigned to the grid
}

