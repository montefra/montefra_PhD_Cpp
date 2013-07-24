/*==========================================================================*/
/* Author: Francesco Montesano, montefra@mpe.mpg.de                         */
/* 06/02/2012,    version: 0.01                                             */
/*                                                                          */
/* set of derived classes which compute the power spectrum using real to    */
/* complex and complex to real FFT transformation as implemente in FFTW     */
/* Although built with the power spectrum in mind, these class can be used  */
/* in different codes requiring the computation of the FFT                  */
/*==========================================================================*/

#include "ps_r2c_c2r.h"

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
ps_r2c_c2r_mpi_inplace::ps_r2c_c2r_mpi_inplace(ptrdiff_t ncells, 
    MPI_Comm com, ptrdiff_t *local_n0, ptrdiff_t *local_0_start) 
{
  fftw_mpi_init();   //initialize mpi
  xcells = ycells = zcells = ncells;    //copy the size of the cell
  totncells = xcells * ycells * zcells;   //compute the totatl number of cells
  zcellso2_1 = zcells/2+1;   //save zcells/2+1 for the correct padding
  comm = com;   //save the MPI communnicator
  //compute the number of cells that each processor contains
  alloc_local = fftw_mpi_local_size_3d(xcells, zcells, zcellso2_1, 
    comm, local_n0, local_0_start);
  local_nx = *local_n0;               //each processor has different number of xcells
  local_x_start = *local_0_start;   //starting point of the x dimension in each processor
  rgrid = fftw_alloc_real(2 * alloc_local);   //allocate the real and complex array
  correction_is_set = false;  //correction for the MAS is not set
  ps_base_set_null_pointer();   //set to zero all the pointers from ps_base that might be used
  whichplan = -999;  //initialise this value: if plan is created it will be changed, otherwise not
}
/* ps_r2c_c2r_mpi_inplace destructor                                        */
ps_r2c_c2r_mpi_inplace::~ps_r2c_c2r_mpi_inplace()
{
  //free_grid();  //deallocate the fftw grid
  if(whichplan != -999) destroy_plan();   //destroy the plan, if ever created
  fftw_forget_wisdom();   //forget the wisdom
  fftw_mpi_cleanup();  //mpi cleanup
  delete[] kcells;   //dealloc the array containing the values of k of the cells
  delete[] psk;     //power spectrum: k
  delete[] psPK;     //power spectrum: sum_modes
  delete[] psn_modes;     //power spectrum: number of modes
}

/*==========================================================================
 * Create the fftw plans                                                            
 * Parameters
 * ----------
 * sign: sign of the FFT transform. Accepted values:
 * 'FFTW_FORWARD'  or -1 for the forward (real to complex) transform 
 * 'FFTW_BACKWARD' or +1 for the backward (complex to real) transform
 * 'FFTW_BOTH' or 0 to initialise both transforms 
 *==========================================================================*/
void ps_r2c_c2r_mpi_inplace::create_plan(int sign)
{
  if(abs(sign)>1){
    std::cerr << "Values accepted for variable 'sign':" << std::endl;
    std::cerr << "'FFTW_FORWARD'  or -1 for the forward (real to complex) transform" << std::endl;
    std::cerr << "'FFTW_BACKWARD' or +1 for the backward (complex to real) transform" << std::endl;
    std::cerr << "'FFTW_BOTH' or 0 to initialise both transforms" << std::endl; 
    exit(10);
    //MPI_Abort(comm, 10);
  }
  whichplan = sign;   //save locally which plans are created
  if(sign == FFTW_FORWARD)   //create the forward plan
    fplan = fftw_mpi_plan_dft_r2c_3d(xcells, ycells, zcells, rgrid, 
        (fftw_complex*)rgrid, comm, FFTW_MEASURE);
  else if(sign == FFTW_BACKWARD)    //create the forward plan
    bplan = fftw_mpi_plan_dft_c2r_3d(xcells, ycells, zcells, (fftw_complex*)rgrid, 
	rgrid, comm, FFTW_MEASURE); 
  else{   //create both
    fplan = fftw_mpi_plan_dft_r2c_3d(xcells, ycells, zcells, rgrid, 
        (fftw_complex*)rgrid, comm, FFTW_MEASURE);
    bplan = fftw_mpi_plan_dft_c2r_3d(xcells, ycells, zcells, (fftw_complex*)rgrid, 
	rgrid, comm, FFTW_MEASURE); 
  }
}


/*==========================================================================
 * Set the mass assigment scheme
 * Parameters
 * ----------
 * MAS: string containing the name of the wanted MAS
 *==========================================================================*/
void ps_r2c_c2r_mpi_inplace::set_MAS(std::string MAS){
  which_MAS=MAS;  //save the name of the MAS
  if(MAS.compare("NGP") == 0)   //NGP
    mas = new NGP_r2c_c2r_mpi(xcells, ycells, zcells, local_nx, local_x_start);
  else if(MAS.compare("CIC") == 0)   //CIC
    mas = new CIC_r2c_c2r_mpi(xcells, ycells, zcells, local_nx, local_x_start);
  else //TSC
    mas = new TSC_r2c_c2r_mpi(xcells, ycells, zcells, local_nx, local_x_start);
}


/*==========================================================================
 * Apply a log(1+delta) transform to the grid
 * delta0*log(1+delta/delta0), with delta0=fabs(min(delta))+epsilon
 * (Seo et al. 2011)
 * Parameters
 * ----------
 * epsilon: small number 
 *==========================================================================*/
void ps_r2c_c2r_mpi_inplace::ln1delta(double epsilon){
  //find first the minimum
  double min = 1.; //e5;  //minimum
  ptrdiff_t i;  //loop integer
  if( epsilon > 0){   //if epsilon positive, find the minimum
  /*  for(i=0; i<2*alloc_local; ++i)
      min = (rgrid[i] < min) ? rgrid[i] : min;   //find the minimum
      */
    min = fabs(min)+epsilon;  //increase the absolute value of the minimum by a small amount 
  }
  else min = 1.;   //if epsilon is exactly 0 do a simple log transform
  for(i=0; i<2*alloc_local; ++i)
    rgrid[i] = min * log(1+rgrid[i]/min);   //as in Seo et al. 2011
}

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
double ps_r2c_c2r_mpi_inplace::sum_grid_real(int root){
  double sum=0, locsum=0;  //sums
  for(ptrdiff_t i=0; i<local_nx; ++i)       //loop over the cells in the x direction
    for(ptrdiff_t j=0; j<ycells; ++j)       //loop over the cells in the y direction
      for(ptrdiff_t l=0; l<zcells; ++l) //loop over the cells in the z direction
	sum += rgrid[l + 2*zcellso2_1 * (j + ycells * i)]; 
  //collect all the sums to 'root' summing them
  MPI_Reduce(&sum, &locsum, 1, MPI_DOUBLE, MPI_SUM, root, comm);
  sum = locsum;   //copy the local sum to the sum
  //broadcast back the total sum to all the processors
  MPI_Bcast(&sum, 1, MPI_DOUBLE, root, comm);   //broadcast sum from root to all

  return(sum);  //return the total sum
}

/*==========================================================================
 * correct each mode by W(k)
 *==========================================================================*/
void ps_r2c_c2r_mpi_inplace::correct_modes(){
  for(ptrdiff_t i=local_nx-1; i>=0; --i){       //loop over the cells in the x direction
    for(ptrdiff_t j=ycells-1; j>=0; --j){       //loop over the cells in the y direction
      for(ptrdiff_t l=zcellso2_1-1; l>=0; --l){ //loop over the cells in the z direction

	//correction of the maf
	double correction = sqrt( corr->corr( kcells[i+local_x_start], kcells[j], kcells[l] ) );   
	rgrid[ 2*(l + zcellso2_1 * (j + ycells*i)) ] /= correction;    //real part of delta(k)
	rgrid[ 2*(l + zcellso2_1 * (j + ycells*i)) + 1 ] /= correction;   //imaginary part of delta(k)
      }
    }
  }
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
void ps_r2c_c2r_mpi_inplace::sum_modes2_sph(const ps_r2c_c2r_mpi_inplace &other_grid, 
    double normalisation, double noise){
  ptrdiff_t kind;   //index of |k| in the power spectrum array
  double modk, temp_modk, d2;   //|k| of the cell and temporary value, |delta(k)|^2 (corrected if requested)
  ptrdiff_t i, j, l;  //loop integers

  for(i=local_nx-1; i>=0; --i){       //loop over the cells in the x direction
    for(j=ycells-1; j>=0; --j){       //loop over the cells in the y direction
      for(l=zcellso2_1-1; l>=0; --l){ //loop over the cells in the z direction
        modk = sqrt( kcells[i+local_x_start]*kcells[i+local_x_start] +
	    kcells[j]*kcells[j] + kcells[l]*kcells[l]);   //|k| of the cell
	if( psnbins>0) temp_modk = modk - pskmin;   //subtract the minimum value of k for the output power spectrum
	else temp_modk = log10(modk/pskmin);  //do the same but in log space
	kind = floor(temp_modk/psdeltak);   //convert the |k| in power spectrum bin units 

	if( kind<0 || kind>=abs(psnbins) ) continue;   //if the value of k is smaller or larger than the range, go to the next loop

	double deltak_re1 = this->rgrid[2*(l + zcellso2_1 * (j + ycells*i))];   //grid1: real part of delta(k) 
	double deltak_re2 = other_grid.rgrid[2*(l + zcellso2_1 * (j + ycells*i))];   //grid2: real part of delta(k)
	double deltak_im1 = this->rgrid[2*(l + zcellso2_1 * (j + ycells*i)) + 1];    //grid1: imaginary part of delta(k)
	double deltak_im2 = other_grid.rgrid[2*(l + zcellso2_1 * (j + ycells*i)) + 1];    //grid2: imaginary part of delta(k)
	d2 = deltak_re1*deltak_re2 + deltak_im1*deltak_im2;    //|delta(k)|^2
        d2 = normalisation*d2 - noise;   //normalise the amplitude and subtract the shot noise
	d2 /= corr->corr(kcells[i+local_x_start], kcells[j], kcells[l]);   //correct for the maf

	psPK[kind] += 2.*d2;
	psn_modes[kind] += 2;
      }
    }
  }
  //if the mode with k=0 is included in the computation, take it out of the sum
  //must be counted only once, while all the others count two for simmetry of
  //the real to complex transform
  if(i == 0 && modk == 0){
    psPK[kind] -= d2;
    psn_modes[kind] -= 1;
  }
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
int ps_r2c_c2r_mpi_inplace::write_fits(std::string froot, int n_proc, 
    std::string comment){
  std::string sfname = froot + "." + to_string(n_proc) + ".fit";  //full file name
  char *cfname = new char[sfname.size()+1];   //char with the fits file name
  std::strcpy( cfname, sfname.c_str() );   //copy the file name to char

  fitsfile *fptr;   //pointer to a fits file
  int status = 0;   //cfitsio status

  fits_create_file(&fptr, cfname, &status);   //create the output fits file

  //create the FITS binary table
  std::vector<std::string> sname;   //vector of strings 
  sname.push_back("test");   //name of the column
  char **name = vector2char(sname);  //convert vector of strings into char**
  sname.clear();  //clean the vector
  sname.push_back("1D");  //add the type
  char **type = vector2char(sname);  //convert vector of strings into char**;
  fits_create_tbl(fptr, BINARY_TBL, 0, 1, name, type, NULL, "rgrid", &status);  //create the table
  //write the part of rgrid in the processor as the unique column of the file
  fits_write_col(fptr, TDOUBLE, 1, 1, 1, 2*alloc_local, rgrid, &status);

  //write the comment of the fits
  char *ccomment = new char[comment.size()+1];   //pointer to the char containing the comment
  std::strcpy(ccomment, comment.c_str());  //transform the comment to char
  fits_write_comment(fptr, ccomment, &status);    //comment
   
  fits_write_date(fptr, &status);    //write the date

  //save the number of cells per dimension
  fits_write_key(fptr, TLONG, "xcells", &xcells, "number of cells in the x dimension", &status);
  fits_write_key(fptr, TLONG, "ycells", &ycells, "number of cells in the y dimension", &status);
  fits_write_key(fptr, TLONG, "zcells", &zcells, "number of cells in the z dimension", &status);
  //save zcells/2+1: only for r2c, c2r transform
  fits_write_key(fptr, TLONG, "zcellso2_1", &zcellso2_1, "zcells/2+1", &status);
  //processor specific variables
  fits_write_key(fptr, TLONG, "local_n0", &local_nx, 
      "number of indeces in the x direction in the processor", &status);
  fits_write_key(fptr, TLONG, "local_0_start", &local_x_start,
      "first index in the x direction in the processor", &status);
  fits_write_key(fptr, TLONG, "alloc_local", &alloc_local, 
      "Half of number of double to be allocated per processor", &status);

  fits_close_file(fptr, &status);    //close the fits file

  if(status != 0) fits_report_error(stdout, status) ;
  return(status);
}

