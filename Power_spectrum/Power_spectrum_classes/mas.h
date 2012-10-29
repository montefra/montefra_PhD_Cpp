/*==========================================================================*/
/* Author: Francesco Montesano, montefra@mpe.mpg.de                         */
/* 06/02/2012,    version: 0.01                                             */
/*                                                                          */
/* MAS for the different cases. Each mas is declared in a class             */
/*==========================================================================*/

#ifndef MAS_H
#define MAS_H

#include <cmath>

#include <fftw3.h>
#include <mpi.h>
#include <fftw3-mpi.h>

/*==========================================================================
 * common declarations and functions for (almost) all the classes
 *==========================================================================*/
class MAS_common
{
  protected:
    ptrdiff_t xcells, ycells, zcells;  //number of cells per dimension
    ptrdiff_t zcellso2_1;   //zcells/2+1: used for the r2c_c2r transform
    ptrdiff_t local_x_start, local_x_finish;   //MPI: start and finish index per processor in the x direction
  
    /*==========================================================================
     * Compute the weights and the cells to which this weights are to be 
     * applied for the CIC MAS
     * Parameters
     * ----------
     * c: position in cells units
     * icells: number of cells in the dimension considered
     * n: number of cells to used for the MAF (for CIC must be 2)
     * output
     * ------
     * cells: 'n' dimension array containing the cells over which the particle 
     *   is spread
     * w: 'n' dimensional array containing the weights per cell
     *==========================================================================*/
    void W_c_CIC(double c, ptrdiff_t ncells, int n, double *w, ptrdiff_t *cells);
    /*==========================================================================
     * Compute the weights and the cells to which this weights are to be 
     * applied for the TSC MAS
     * Parameters
     * ----------
     * c: position in cells units
     * icells: number of cells in the dimension considered
     * n: number of cells to used for the MAF (for TSC must be 3)
     * output
     * ------
     * cells: 'n' dimension array containing the cells over which the particle 
     *   is spread
     * w: 'n' dimensional array containing the weights per cell
     *==========================================================================*/
    void W_c_TSC(double c, ptrdiff_t ncells, int n, double *w, ptrdiff_t *cells);
};

/*==========================================================================
 * Virtual base class for all the *r2c_c2r_mpi
 *==========================================================================*/
class MAS_r2c_c2r_mpi: public MAS_common
{
  public:
    virtual void MAS(double *rgrid, double x, double y, double z, double w) =0;
    virtual void MAS(fftw_complex *cgrid, double x, double y, double z, int ind, double w) =0;
};

/*==========================================================================
 * NGP mass assigment scheme for r2c_c2r transform with MPI
 *==========================================================================*/
class NGP_r2c_c2r_mpi: public MAS_r2c_c2r_mpi
{

  public:
    /*==========================================================================
     * Class constructor
     * Parameters
     * ----------
     * xcell, ycell, zcell: number of cells in the x, y and z direction
     * local_n0: first cell in the x direction per processor
     * local_0_start: number of cells in the x direction per processor          
     *==========================================================================*/
    NGP_r2c_c2r_mpi(ptrdiff_t xcell, ptrdiff_t ycell, ptrdiff_t zcell, 
	ptrdiff_t local_n0, ptrdiff_t local_0_start){
      xcells = xcell;
      ycells = ycell;
      zcells = zcell;
      zcellso2_1 = zcells/2+1;
      local_x_start = local_0_start;
      local_x_finish = local_n0+local_0_start;
    }
    /*==========================================================================
     * NGP MAS for the real grid
     * Parameters
     * ----------
     * rgrid: double FFTW grid
     * x, y, z: particle position in grid units
     * w: weight of the particle
     *==========================================================================*/
    void MAS(double *rgrid, double x, double y, double z, double w)  
    {
      //take the floor of the positions
      ptrdiff_t fx = floor(x), fy = floor(y), fz = floor(z);
      if(fx >= local_x_start && fx < local_x_finish)
	rgrid[ fz + 2*zcellso2_1 * (fy + ycells * (fx-local_x_start)) ] += w;
    }
    /*==========================================================================
     * NGP MAS for the complex grid
     * Parameters
     * ----------
     * cgrid: complex FFTW grid
     * x, y, z: particle position in grid units
     * ind: 0 or 1 to assigne the particle to the real or immaginary part of the grid
     * w: weight of the particle
     *==========================================================================*/
    void MAS(fftw_complex *cgrid, double x, double y, double z, int ind, double w)
    {
      //take the floor of the positions
      ptrdiff_t fx = floor(x), fy = floor(y), fz = floor(z);
      if( (fx >= local_x_start && fx < local_x_finish) && fz < zcellso2_1 ) 
	cgrid[ fz + zcellso2_1 * (fy + ycells * (fx-local_x_start)) ][ind] += w;
    }
};

/*==========================================================================
 * CIC mass assigment scheme for r2c_c2r transform with MPI
 *==========================================================================*/
class CIC_r2c_c2r_mpi: public MAS_r2c_c2r_mpi
{

  public:
    /*==========================================================================
     * Class constructor
     * Parameters
     * ----------
     * xcell, ycell, zcell: number of cells in the x, y and z direction
     * local_n0: first cell in the x direction per processor
     * local_0_start: number of cells in the x direction per processor          
     *==========================================================================*/
    CIC_r2c_c2r_mpi(ptrdiff_t xcell, ptrdiff_t ycell, ptrdiff_t zcell, ptrdiff_t local_n0, ptrdiff_t local_0_start){
      xcells = xcell;
      ycells = ycell;
      zcells = zcell;
      zcellso2_1 = zcells/2+1;
      local_x_start = local_0_start;
      local_x_finish = local_n0+local_0_start;
    }
    /*==========================================================================
     * CIC MAS for the real grid
     * Parameters
     * ----------
     * rgrid: double FFTW grid
     * x, y, z: particle position in grid units
     * w: weight of the particle
     *==========================================================================*/
    void MAS(double *rgrid, double x, double y, double z, double w);
    /*==========================================================================
     * CIC MAS for the complex grid
     * Parameters
     * ----------
     * cgrid: complex FFTW grid
     * x, y, z: particle position in grid units
     * ind: 0 or 1 to assigne the particle to the real or immaginary part of the grid
     * w: weight of the particle
     *==========================================================================*/
    void MAS(fftw_complex *cgrid, double x, double y, double z, int ind, double w);
};

/*==========================================================================
 * TSC mass assigment scheme with MPI
 *==========================================================================*/
class TSC_r2c_c2r_mpi: public MAS_r2c_c2r_mpi
{
  public:
    /*==========================================================================
     * Class constructor
     * Parameters
     * ----------
     * xcell, ycell, zcell: number of cells in the x, y and z direction
     * local_n0: first cell in the x direction per processor
     * local_0_start: number of cells in the x direction per processor          
     *==========================================================================*/
    TSC_r2c_c2r_mpi(ptrdiff_t xcell, ptrdiff_t ycell, ptrdiff_t zcell, ptrdiff_t local_n0, ptrdiff_t local_0_start){
      xcells = xcell;
      ycells = ycell;
      zcells = zcell;
      zcellso2_1 = zcells/2+1;
      local_x_start = local_0_start;
      local_x_finish = local_n0+local_0_start;
    }
    /*==========================================================================
     * TSC MAS for the real grid
     * Parameters
     * ----------
     * rgrid: double FFTW grid
     * x, y, z: particle position in grid units
     * w: weight of the particle
     *==========================================================================*/
    void MAS(double *rgrid, double x, double y, double z, double w);
    /*==========================================================================
     * TSC MAS for the complex grid
     * Parameters
     * ----------
     * cgrid: complex FFTW grid
     * x, y, z: particle position in grid units
     * ind: 0 or 1 to assigne the particle to the real or immaginary part of the grid
     * w: weight of the particle
     *==========================================================================*/
    void MAS(fftw_complex *cgrid, double x, double y, double z, int ind, double w);
};
#endif
