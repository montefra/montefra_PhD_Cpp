/*==========================================================================*/
/* Author: Francesco Montesano, montefra@mpe.mpg.de                         */
/* 06/02/2012,    version: 0.01                                             */
/*                                                                          */
/* Base class for the computation of the power spectrum                     */
/* contains definitions and declarations common to different                */
/* implementations. The computation is based on the FFTW                    */
/* Although built with the power spectrum in mind, these class can be used  */
/* in different codes requiring the computation of the FFT                  */
/*==========================================================================*/

#ifndef PS_BASE_H
#define PS_BASE_H

#include <cmath>
#include <cstring>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include <mpi.h>
#include <fftw3.h>
#include <fftw3-mpi.h>

//mass corrections
#include "mas_correction.h"

//value used to initialised both the forward and backward FFTW plans
#define FFTW_BOTH 0
//define pi if not already defined
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

class ps_base
{

  protected:

    ptrdiff_t xcells, ycells, zcells;    //number of cells per dimension
    ptrdiff_t totncells;   //total number of cells
    MPI_Comm comm;  //MPI communicator
    //processor specific variables: first cells and their number in the x direction
    ptrdiff_t local_nx, local_x_start;   

    ptrdiff_t alloc_local;  //size of space to be allocated per processor

    fftw_plan fplan, bplan;       // forward and backward plan for the FFT
    int whichplan; //sign of the FFT transform. Accepted values: -1,0,1 (see create_plan)

    double kN, deltak;   //Nyquist wavenumber and cell size in k space
    double *kcells;    // values of k of the FFT cells

    //averaged power spectrum variables
    ptrdiff_t psnbins;  //number of bins in the input power spectrum
    double psdeltak, pskmax, pskmin;   //maximum and minimum value of k in the output power spectrum
    ptrdiff_t pskmax_ind, pskmin_ind;   //ceil and floor index of kmax and kmin in the FFT grid
    double *psk, *psPK;   //pointers to array containing k, sum_modes
    long *psn_modes;  //pointers to array containing nmodes

    std::string which_MAS;   //Save the name of the MAS used
    bool correction_is_set;  //correction for the MAS is set or not
    correction_base *corr;  //pointer to a virtual base class for the correction of the mas

    /*==========================================================================
     * destroy the plan(s) 
     *==========================================================================*/
    void destroy_plan()
    {
      if(whichplan == FFTW_FORWARD || whichplan == FFTW_BOTH) fftw_destroy_plan(fplan);
      if(whichplan == FFTW_BACKWARD || whichplan == FFTW_BOTH) fftw_destroy_plan(bplan);
    }
    /*==========================================================================
     * Set to 'NULL' all the pointers declared in this class
     * used in order to do save delete even if the array non allocated
     *==========================================================================*/
    void ps_base_set_null_pointer(){
      kcells = NULL;   //values of k of the cells
      psk = NULL;     //power spectrum: k
      psPK = NULL;     //power spectrum: sum_modes
      psn_modes = NULL;     //power spectrum: number of modes
    }

    /*==========================================================================
     * allocate the power spectrum arrays and initialise them setting them to 0
     *==========================================================================*/
    void alloc_initialise_psk();

  public:

    virtual ~ps_base(){};   //virtual destructor

    /*==========================================================================
     * read the FFTW widsom: single processor
     * Parameters
     * ----------
     * fname: file that should contain the wisdom
     * output
     * ------
     * 1 for success, 0 for falure
     *==========================================================================*/
    int read_wisdom(std::string fname){
      return( fftw_import_wisdom_from_filename(fname.c_str()) != 0); 
    }
    /*==========================================================================
     * read the FFTW widsom and broadcast it: MPI
     * Parameters
     * ----------
     * fname: file that should contain the wisdom
     * myrank: processor rank
     * comm: MPI communicator
     * output
     * ------
     * on myrank=0: 1 for success, 0 for falure; 1 for all the other processors
     *==========================================================================*/
    int read_wisdom(std::string fname, int myrank, MPI_Comm comm){
      int err = 1;  //error message from reading the wisdom
      if(myrank == 0) err = fftw_import_wisdom_from_filename(fname.c_str());
      fftw_mpi_broadcast_wisdom(comm);   //broadcast the wisdom is has been read
      return( err );
    }
    /*==========================================================================
     * save the FFTW widsom: single processor
     * Parameters
     * ----------
     * fname: file that should contain the wisdom
     * output
     * ------
     * 1 for success, 0 for falure
     *==========================================================================*/
    int write_wisdom(std::string fname){
      return( fftw_export_wisdom_to_filename(fname.c_str()) );
    }
    /*==========================================================================
     * gather and save the FFTW widsom and broadcast it: MPI
     * Parameters
     * ----------
     * fname: file that should contain the wisdom
     * myrank: processor rank
     * comm: MPI communicator
     * output
     * ------
     * on myrank=0: 1 for success, 0 for falure; 1 for all the other processors
     *==========================================================================*/
    int write_wisdom(std::string fname, int myrank, MPI_Comm comm){
      int err = 1;  //error message from reading the wisdom
      fftw_mpi_gather_wisdom(comm);
      if(myrank == 0) err = fftw_export_wisdom_to_filename(fname.c_str());
      return( err );
    }

    /*==========================================================================
     * Set the values of k for each point in the FFT grid 
     * Parameters
     * ----------
     * L: physical or comuving size of the box 
     * output
     * ------
     * cell: size of the FFTW cell in the input space (L/n_cells)
     *==========================================================================*/
    double set_k(double L);
    /*==========================================================================
     * Set the values of k the output power spectrum
     * Parameters
     * ----------
     * kmax, kmin: maximum and minimum k for the output power spectrum
     * nbins: signed number of bins for the output power spectrum 
     *   (positive: linear; negative: logarithmic
     * output
     * ------
     * none
     *==========================================================================*/
    void set_psk(double kmax, double kmin, ptrdiff_t nbins);

    /*==========================================================================
     * Set the correction for the MAS aliases
     * Parameters
     * ----------
     * cor: choise of the correction:
     *   0: no correction
     *   1: correct dividing by W^2(k)
     *   1: correct dividing by sum_n W^2(k+2*n*kN)
     *==========================================================================*/
    void set_MAS_correction(int cor);

    /*==========================================================================
     * Execute the FFTW plan
     * Works only for FFTW_FORWARD (-1) or FFTW_BACKWARD (+1)
     *==========================================================================*/
    void execute_plan();
    /*==========================================================================
     * Execute the FFTW plan. If only one initialised, 'sign' must correspond
     * Parameters
     * ----------
     * sign: sign of the FFT transform to execute
     *   'FFTW_FORWARD'  or -1 for the forward (real to complex) transform 
     *   'FFTW_BACKWARD' or +1 for the backward (complex to real) transform
     *==========================================================================*/
    void execute_plan(int sign);

    /*==========================================================================
    * single processor: compute the average of |delta|^2 per bin and 
    * save an ascii file with the spherical averaged power spectrum. 
    * The output file has the structure: k, P(k), P(k)+shot noise
    * Parameters
    * ----------
    * ofile: name of the output file
    * normalisation: normalisation to apply to the power spectrum
    * noise: shot noise amplitude
    * header: (multiline) string to add at the beginning of file. 
    *   Default: "none", no header
    *==========================================================================*/
    void savePK(std::string ofile, double normalisation, double noise,
        std::string header="none");
    /*==========================================================================
     * MPI: broadcast the sum of the deltas to the root processor, compute the
     * average of |delta|^2 per bin and save an ascii file with the spherical
     * averaged power spectrum. The output file has the structure:
     * k, P(k), P(k)+shot noise
     * Parameters
     * ----------
     * ofile: name of the output file
     * normalisation: normalisation to apply to the power spectrum
     * noise: shot noise amplitude
     * myrank: rank of the processor
     * root: root processor
     * com: MPI communicator
     *==========================================================================*/
    void savePK(std::string ofile, double normalisation, double noise, 
	int myrank, int root, MPI_Comm com);
    /*==========================================================================
     * MPI: broadcast the sum of the deltas to the root processor, compute the
     * average of |delta|^2 per bin and save an ascii file with the spherical
     * averaged power spectrum. The output file has the structure:
     * k, P(k), P(k)+shot noise
     * Parameters
     * ----------
     * ofile: name of the output file
     * normalisation: normalisation to apply to the power spectrum
     * noise: shot noise amplitude
     * header: (multiline) string to add at the beginning of file. 
     * myrank: rank of the processor
     * root: root processor
     * com: MPI communicator
     *==========================================================================*/
    void savePK(std::string ofile, double normalisation, double noise,
        std::string header, int myrank, int root, MPI_Comm com);

};

/*==========================================================================
 * Convert a vector of strings in a char**
 * Parameters
 * ----------
 * vec: vector of strings
 * output
 * ------
 * chars: array of char arrays
 *==========================================================================*/
char** vector2char(std::vector<std::string> vec);

#endif
