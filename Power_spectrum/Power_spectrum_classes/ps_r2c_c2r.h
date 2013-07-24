/*==========================================================================*/
/* Author: Francesco Montesano, montefra@mpe.mpg.de                         */
/* 06/02/2012,    version: 0.01                                             */
/*                                                                          */
/* set of derived classes which compute the power spectrum using real to    */
/* complex and complex to real FFT transformation as implemente in FFTW     */
/* Although built with the power spectrum in mind, these class can be used  */
/* in different codes requiring the computation of the FFT                  */
/*==========================================================================*/

#ifndef PS_R2C_C2R_H
#define PS_R2C_C2R_H

#include "ps_base.h"

#include <vector>
//mas 
#include "mas.h"
//fits file in/out
#include <fitsio.h>
//common functions
#include "../common_functions.h"

/*==========================================================================*/
class ps_r2c_c2r: public ps_base
/* base class for the complex to real and real to complex FFT               */
{
  protected:
    ptrdiff_t zcellso2_1;   //contains zcells/2+1. Used to get the correct padding 
                          //of the real grid
  public:
    virtual ~ps_r2c_c2r(){};   //virtual destructor

    double *rgrid;        //real grid
    fftw_complex *cgrid;  //fftw complex grid 
};
/*==========================================================================*/


/*==========================================================================*/
class ps_r2c_c2r_mpi_inplace: public ps_r2c_c2r
/* This class implement the real to complex and complex to real in place    */
/* transform with MPI                                                       */
{

  private:
    MAS_r2c_c2r_mpi *mas;  //mass assigment scheme virtual base class. Pointer set up in 'set_MAS'

  public:

    /*==========================================================================
     * ps_r2c_c2r_mpi_inplace class constructor
     * Parameters
     * ----------
     * ncells: number of FFT cells assuming a square 
     * com: MPI communicator
     *
     * output
     * ------
     * local_n0: first cell in the x direction per processor
     * local_0_start: number of cells in the x direction per processor          
     *==========================================================================*/
    ps_r2c_c2r_mpi_inplace(ptrdiff_t ncells, MPI_Comm com,    //class constructor
	ptrdiff_t *local_n0, ptrdiff_t *local_0_start);
    /* class destructor                                                         */
    ~ps_r2c_c2r_mpi_inplace();  //destructor


    /*==========================================================================
     * Operators
     *==========================================================================*/
    /*==========================================================================
     * Compound assignment
     *==========================================================================*/
    /*==========================================================================
     * add LHS into RHS
     * Parameters
     * ----------
     *  rhs: ps_r2c_c2r_mpi_inplace instance to be added
     * output
     * ------
     *  ps_r2c_c2r_mpi_inplace with the object and the added grid
     *==========================================================================*/
    ps_r2c_c2r_mpi_inplace & operator+=(const ps_r2c_c2r_mpi_inplace &rhs){
      for(ptrdiff_t i=0; i<2*alloc_local; ++i)
        this->rgrid[i] += rhs.rgrid[i];
      return *this;
    }
    /*==========================================================================
     * multiply RHS by constant
     * Parameters
     * ----------
     *  rhs: ps_r2c_c2r_mpi_inplace instance to be added
     * output
     * ------
     *  ps_r2c_c2r_mpi_inplace with the object and the added grid
     *==========================================================================*/
    ps_r2c_c2r_mpi_inplace & operator*=(const double rhs){
      for(ptrdiff_t i=0; i<2*alloc_local; ++i)
        this->rgrid[i] *= rhs;
      return *this;
    }
    /*==========================================================================
     * divide RHS by constant
     * Parameters
     * ----------
     *  rhs: ps_r2c_c2r_mpi_inplace instance to be added
     * output
     * ------
     *  ps_r2c_c2r_mpi_inplace with the object and the added grid
     *==========================================================================*/
    ps_r2c_c2r_mpi_inplace & operator/=(const double rhs){
      for(ptrdiff_t i=0; i<2*alloc_local; ++i)
        this->rgrid[i] /= rhs;
      return *this;
    }

    /*==========================================================================
     * free the fftw grid
     *==========================================================================*/
    void free_grid(){ fftw_free(rgrid); }

    /*==========================================================================
     * Create the fftw plans                                                            
     * Parameters
     * ----------
     * sign: sign of the FFT transform. Accepted values:
     * 'FFTW_FORWARD'  or -1 for the forward (real to complex) transform 
     * 'FFTW_BACKWARD' or +1 for the backward (complex to real) transform
     * 'FFTW_BOTH' or 0 to initialise both transforms 
     *==========================================================================*/
    void create_plan(int sign);

    /*==========================================================================
     * Initialise the grid setting all the elements to 0                        
     *==========================================================================*/
    void initialise() {
      for(ptrdiff_t i=0; i<2*alloc_local; ++i) rgrid[i] = 0.; 
    }

    /*==========================================================================
     * Set the mass assigment scheme
     * Parameters
     * ----------
     * MAS: string containing the name of the wanted MAS
     *==========================================================================*/
    void set_MAS(std::string MAS);

    /*==========================================================================
     * Assign the particle to the grid according to the chosen MAS
     * first version for the real grid, second for the complex one
     * Parameters
     * ----------
     * x, y, z: position in grid units
     * ind (for the complex grid only): 0 or 1 to assigne the particle to the
     *   real or immaginary part of the grid
     * w: weight of the particle
     *==========================================================================*/
    void assign_particle(double x, double y, double z, double w){
      mas->MAS(rgrid, x, y, z, w);   //call the function to do the assigment
    }
    void assign_particle(double x, double y, double z, int ind, double w){
      mas->MAS((fftw_complex*)rgrid, x, y, z, ind, w);   //call the function to do the assigment
    }
    /*==========================================================================
     * Assign the particle to the grid according to the chosen MAS enforcing 
     * periodic boundary conditions
     * first version for the real grid, second for the complex one
     * Parameters
     * ----------
     * x, y, z: position in grid units
     * ind (for the complex grid only): 0 or 1 to assigne the particle to the
     *   real or immaginary part of the grid
     * w: weight of the particle
     *==========================================================================*/
    void assign_particle_periodic(double x, double y, double z, double w){
      if( x < 0 ) x += xcells;
      else if( x>= xcells ) x-= xcells;
      if( y < 0 ) y += ycells;
      else if( y>= ycells ) y-= ycells;
      if( z < 0 ) z += zcells;
      else if( z>= zcells ) z-= zcells;
      mas->MAS_periodic(rgrid, x, y, z, w);   //call the function to do the assigment
    }
    void assign_particle_periodic(double x, double y, double z, int ind, double w){
      if( x < 0 ) x += xcells;
      else if( x>= xcells ) x-= xcells;
      if( y < 0 ) y += ycells;
      else if( y>= ycells ) y-= ycells;
      if( z < 0 ) z += zcells;
      else if( z>= zcells ) z-= zcells;
      mas->MAS_periodic((fftw_complex*)rgrid, x, y, z, ind, w);   //call the function to do the assigment
    }
    
    /*==========================================================================
     * compute the mean of the particles in the grid and substitute the number 
     * of particle per cell by n/n_mean - 1
     * Parameters
     * ----------
     * N: number of particles assigned to the full grid
     * output
     * ------
     * mean: N/n_cells
     *==========================================================================*/
    double to_overdensity(ptrdiff_t N){
      double mean = (double)N / totncells;    //mean number density
      for(ptrdiff_t i=0; i<2*alloc_local; ++i) rgrid[i] = rgrid[i] / mean - 1.;
      return( mean );
    }

    /*==========================================================================
     * compute the F(r) density field from the as rgrid - alpha*randomgrid
     * Parameters
     * ----------
     * rangrid: FFTW grid of double with the same size of rgrid containin the 
     *   random points
     * alpha: sum(w_g)/sum(w_r)
     *==========================================================================*/
    void to_Fr(const ps_r2c_c2r_mpi_inplace &rangrid, double alpha){
      to_Fr(rangrid, alpha, 1.);
    }

    /*==========================================================================
     * compute the F(r) density field from the as (rgrid - alpha*randomgrid)*N
     * Parameters
     * ----------
     * rangrid: FFTW grid of double with the same size of rgrid containin the 
     *   random points
     * alpha: sum(w_g)/sum(w_r)
     * N: multiplicative mormalisation 1/sqrt(sum(n(z)*w^2))
     *==========================================================================*/
    void to_Fr(const ps_r2c_c2r_mpi_inplace &rangrid, double alpha, double N){
      for(ptrdiff_t i=0; i<2*alloc_local; ++i)
        this->rgrid[i] = (this->rgrid[i] - alpha * rangrid.rgrid[i])/N;
    }

    /*==========================================================================
     * Apply a log(1+delta) transform to the grid
     * log(1+delta), if delta>0, delta otherwise (Neyrinck, 2011)
     *==========================================================================*/
    void ln1delta(){
      for(ptrdiff_t i=0; i<2*alloc_local; ++i)
	rgrid[i] = (rgrid[i] > 0) ? log(1.+rgrid[i]) : rgrid[i];
    }
    /*==========================================================================
     * Apply a log(1+delta) transform to the grid
     * delta0*log(1+delta/delta0), with delta0=fabs(min(delta))+epsilon
     * (Seo et al. 2011)
     * Parameters
     * ----------
     * epsilon: small number 
     *==========================================================================*/
    void ln1delta(double epsilon);

    /*==========================================================================
     * compute the sum of all the values of the grid. Respect the padding of
     * the real grid
     * Parameters
     * ----------
     * root: processor number where to collect the total sum 
     *   to rebroadcast to all the processors
     * output
     * ------
     * sum: total sum
     *==========================================================================*/
    double sum_grid_real(int root);

    /*==========================================================================
     * Compute the mean of the values in the cells
     * Parameters
     * ----------
     * root: processor number where to collect the total sum 
     *   to rebroadcast to all the processors
     * output
     * ------
     * sum: total sum
     *==========================================================================*/
    double mean_cells(int root){
      double summodes = sum_grid_real(root);  //execute the sum over all the cells
      ptrdiff_t ncells = xcells * ycells * zcells;  //total number of cells
      return( summodes/ncells );   //return the mean of the modes 
    }

    /*==========================================================================
      - alpha * rangrid[i])*N;* Subtract from each cell the mean of all the cells
     * Parameters
     * ----------
     * root: processor number where to collect the total sum 
     *   to rebroadcast to all the processors
     *==========================================================================*/
    void rescale_mean(int root){
      double mu = mean_cells(root);  //get the mean over all the cells
      for(ptrdiff_t i=0; i<2*alloc_local; ++i)
        rgrid[i] -= mu;  //subtract the total mean from each cell
    }

    /*==========================================================================
     * correct each mode by W(k)
     *==========================================================================*/
    void correct_modes();

    /*==========================================================================
     * divide each element in the grid by the total number of cells. To be used
     * when renormalising the power spectrum after a double FFT transform
     *==========================================================================*/
    void normalize_double_transform(){
      for(ptrdiff_t i=0; i<2*alloc_local; ++i)
	rgrid[i] /= totncells;
    }

    /*==========================================================================
     * compute the sum of |delta(k)|^2 in spherical shells within the required 
     * bins, set by 'set_psk', and count the corresponding number of modes
     * this version normalise and substract the shot noise from every mode before
     * correcting
     * Parameters
     * ----------
     *  normalisation: multiplicative normalisation for |delta(k)^2| (not noise)
     *  noise: shot noise amplitude to subtract from normalisation*|delta(k)^2|
     *==========================================================================*/
    void sum_modes2_sph(double normalisation, double noise){
      sum_modes2_sph(*this, normalisation, noise);
    }
    /*==========================================================================
     * compute the sum of |delta(k)|^2 in spherical shells within the required 
     * bins, set by 'set_psk', and count the corresponding number of modes
     * no normalisation or shot noise applied
     *==========================================================================*/
    void sum_modes2_sph(){
      sum_modes2_sph(1., 0.);
    }
    /*==========================================================================
     * compute the sum of re[delta1(k)*delta2(k)] in spherical shells within 
     * the required bins, set by 'set_psk', and count the corresponding 
     * number of modes
     * Parameters
     * ----------
     *  other_grid: grid object containing delta2
     *  normalisation: multiplicative normalisation for |delta(k)^2| (not noise)
     *  noise: shot noise amplitude to subtract from normalisation*|delta(k)^2|
     *==========================================================================*/
    void sum_modes2_sph(const ps_r2c_c2r_mpi_inplace &other_grid, double
        normalisation, double noise);
    /*==========================================================================
     * compute the sum of re[delta1(k)*delta2(k)] in spherical shells within 
     * the required bins, set by 'set_psk', and count the corresponding 
     * number of modes
     * Parameters
     * ----------
     * other_grid: grid object containing delta2
     *==========================================================================*/
    void sum_modes2_sph(const ps_r2c_c2r_mpi_inplace &other_grid){
      sum_modes2_sph(other_grid, 1., 0.);
    }

    /*==========================================================================
     * save rgrid into a fits file. Each processor write a single file with a 
     * copy of the grid allocated locally. 
     * The number of cells in each direction, the half+1 of the last direction,
     * the start and number of elements per processory in the first direction
     * and the amount of space to be allocated per processor
     * Parameters
     * ----------
     * froot: root of the output fits file. The file is froot.n_proc.fit
     * n_proc: processor number
     * comment: string containing the commend to write in the fits file
     * output
     * ------
     * status: cfitsio output status. If the status is non 0 a status report
     * is written on stderr
     *==========================================================================*/
    int write_fits(std::string froot, int n_proc, 
	std::string comment);

};
/*==========================================================================*/

#endif

