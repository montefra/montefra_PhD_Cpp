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

#include "ps_base.h"

/*==========================================================================
 * Set the values of k for each point in the FFT grid 
 * Parameters
 * ----------
 * L: physical or comuving size of the box 
 * output
 * ------
 * cellsize: size of the FFTW cell in the input space (L/n_cells)
 *==========================================================================*/
double ps_base::set_k(double L){
  double cellsize = L/xcells;   //size of each cell

  kN = M_PI/cellsize;   //Nyquist wavesnumber
  deltak = 2.*M_PI/L;   //cell size in kspace 

  kcells = new double[xcells];   //alloc the k corresponding to the cells
  kcells[0] = 0.;   //set to 0 the first element of the k array
  ptrdiff_t ii;  //loop integer
  //fill almost all the array with k with kcells[i]=-kcells[n-i] for 1<=i<xyzcells[0]/2
  for(ii=1; ii<xcells/2; ++ii){
    kcells[ii] = deltak*ii;
    kcells[xcells-ii] = -kcells[ii]; 
  }
  if(xcells%2 == 0) //if the number of cells is even, the xcells/2 element of k is kN 
    kcells[xcells/2] = kN;   
  else{ //otherwise kcells[xcells/2] = kcells[xcells/2] and they are both < kN
    kcells[xcells/2] = deltak*xcells/2;     
    kcells[xcells/2+1] = -kcells[xcells/2];
  }
  return(cellsize);
}
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
void ps_base::set_psk(double kmax, double kmin, ptrdiff_t nbins){
  psnbins = nbins;  //number of bins in the input power spectrum
  if( kmax < 0 ) pskmax = kN;   //maximum value of k in the output power spectrum
  else pskmax = kmax;
  pskmin = kmin;   //minimum value of k in the output power spectrum

  pskmax_ind = ceil(pskmax/deltak);  //save the indeces to which kmax and kmin belongs 
  pskmin_ind = floor(pskmin/deltak); //in the FFT grid

  //allocate the three power spectrum array
  alloc_initialise_psk();

  if(psnbins < 0){   //logarithimic binning
    kmax = log10(kmax);
    kmin = log10(kmin);
  }
  ptrdiff_t anbins = abs(nbins);  //temporary save the absolute value of the number of bins
  psdeltak = (pskmax - pskmin)/anbins;  //k binning
  ptrdiff_t i; //loop integer
  for(i=0; i< anbins; ++i) psk[i] = psdeltak*(i+1./2.) + pskmin;   //meadian of the bin
  if(psnbins < 0)   //logarithimic binning
    for(i=0; i< anbins; ++i) psk[i] = pow(10., psk[i]);   //convert k to 10^k
}
/*==========================================================================
 * allocate the power spectrum arrays and initialise them setting them to 0
 *==========================================================================*/
void ps_base::alloc_initialise_psk(){
  psk = new double[abs(psnbins)];   //alloc the k of the power spectrum
  psPK = new double[abs(psnbins)];   //alloc the P(k)
  psn_modes = new long[abs(psnbins)];   //alloc the number of modes per bin
  for(ptrdiff_t i=0; i<abs(psnbins); ++i){   //initialise the power spectrum and nmodes to 0
    psPK[i] = 0.;
    psn_modes[i] = 0;
  }
}

/*==========================================================================
 * Set or reset the correction for the MAS aliases
 * Parameters
 * ----------
 * cor: choise of the correction:
 *   0: no correction
 *   1: correct dividing by W^2(k)
 *   2: correct dividing by sum_n W^2(k+2*n*kN)
 *==========================================================================*/
void ps_base::set_MAS_correction(int cor){
  if( correction_is_set == true )   //if the correction had already been set, 
    delete corr;    //deallocate the pointer
  else
    correction_is_set = true;   //if not, now will be set

  if( which_MAS.compare("NGP") == 0 ){   //NGP
    if( cor == 1 )   //correction 1
      corr = new W_k_correction(kN, 1); 
    else   //no correction or correction 2
      corr = new no_correction; 
  }
  else if( which_MAS.compare("CIC") == 0 ){   //CIC
    if( cor == 0 )  //no correction
      corr = new no_correction; 
    else if( cor == 1 )  //correction 1
      corr = new W_k_correction(kN, 2);
    else   //correction 2
      corr = new W_k_2nkN_correction_CIC(kN);
  }
  else if( which_MAS.compare("TSC") == 0 ){   //TSC
    if( cor == 0 )  //no correction
      corr = new no_correction; 
    else if( cor == 1 )  //correction 1
      corr = new W_k_correction(kN, 3);
    else   //correction 2
      corr = new W_k_2nkN_correction_TSC(kN);
  }
}

/*==========================================================================
 * Execute the FFTW plan
 * Works only for FFTW_FORWARD (-1) or FFTW_BACKWARD (+1)
 *==========================================================================*/
void ps_base::execute_plan(){
  if(whichplan == FFTW_BOTH){
    std::cerr << "The plan has been created with 'FFTW_BOTH' (0) and now ";
    std::cerr << "I don't now what to do. Use 'void execute_plan(int sign)'";
    std::cerr << "instead to be able to chose among the two options" << std::endl;
    exit(1);
    //MPI_Abort(comm, 1);
  }
  else if(whichplan == FFTW_FORWARD) fftw_execute(fplan);
  else fftw_execute(bplan);
}
/*==========================================================================
 * Execute the FFTW plan. If only one initialised, 'sign' must correspond
 * Parameters
 * ----------
 * sign: sign of the FFT transform to execute
 *   'FFTW_FORWARD'  or -1 for the forward (real to complex) transform 
 *   'FFTW_BACKWARD' or +1 for the backward (complex to real) transform
 *==========================================================================*/
void ps_base::execute_plan(int sign){
  if( sign == 0 ){
    std::cerr << "Please decide which plan you want to execute: " << std::endl;
    std::cerr << "'FFTW_FORWARD'  or -1 for the forward (real to complex) transform"
	<< std::endl;
    std::cerr << "'FFTW_BACKWARD' or +1 for the backward (complex to real) transform"
	<< std::endl;
    exit(2);
    //MPI_Abort(comm, 2);
  }
  if(whichplan != 0){   //if only one plan created
    if(sign == whichplan) execute_plan();  //if the created and requested plan corresponds
    else{
      std::cerr << "The plan created and the one requested for execution do not correspond"
	  << std::endl;
      exit(3);
      //MPI_Abort(comm, 3);
    }
  }
  else{
    if(sign == FFTW_FORWARD) fftw_execute(fplan);
    else fftw_execute(bplan);
  }
}


/*==========================================================================
 * single processor: compute the average of |delta|^2 per bin and 
 * save an ascii file with the spherical averaged power spectrum. 
 * The output file has the structure: k, P(k), P(k)+shot noise
 * Parameters
 * ----------
 * ofile: name of the output file
 * normalisation: normalisation to apply to the power spectrum
 * noise: shot noise amplitude
 *==========================================================================*/
void ps_base::savePK(std::string ofile, double normalisation, double noise)
{
  std::ofstream out(ofile.c_str());   //open the output file
  out << "#	k	P(k)	P(k)+noise" << std::endl;   //write an header
  //set the precision for the output 
  out.setf(std::ios_base::scientific);
  out.precision(6);
  out.width(9);

  for(ptrdiff_t i=0; i<abs(psnbins); ++i){
    if(psn_modes[i] == 0)   //if there are no modes in the bins, print 0
      out << psk[i] << "\t" << 0. << "\t" << 0. << std::endl;
    else{   // if there are modes
      psPK[i] = normalisation * psPK[i] / psn_modes[i];  //do the averange
      out << psk[i] << "\t" << psPK[i] - noise << "\t" << psPK[i] << std::endl;
    }
  }
  out.close();  //close the file
}
/*==========================================================================
 * MPI: broadcast the sum of the deltas to the root processor, compute the
 * average of |delta|^2 per bin and save an ascii file with the spherical
 * averaged power spectrum. The output file has the structure:
 * k, P(k), P(k)+shot noise
 * Parameters
 * ----------
 * ofile: name of the output file
 * normalisation: multiplicative normalisation to apply to the power spectrum
 * noise: shot noise amplitude
 * myrank: rank of the processor
 * root: root processor
 * com: MPI communicator
 *==========================================================================*/
void ps_base::savePK(std::string ofile, double normalisation, double noise, 
    int myrank, int root, MPI_Comm com)
{
  double *psPK_root;   //arrays in the root to which sum |delta|^2 is broadcasted
  long *psn_modes_root;    //arrays in the root to which number of modes is broadcasted
  if(myrank==root){   //allocate these array in the root
    psPK_root = new double[abs(psnbins)];   //alloc the P(k)
    psn_modes_root = new long[abs(psnbins)];   //alloc the number of modes per bin
  }
  //sum all the delta^2 and number of modes in the root processor
  MPI_Reduce(psPK, psPK_root, abs(psnbins), MPI_DOUBLE, MPI_SUM, root, com);
  MPI_Reduce(psn_modes, psn_modes_root, abs(psnbins), MPI_LONG, MPI_SUM, root, com);

  if(myrank==root){   //only in the root

    std::ofstream out(ofile.c_str());   //open the output file
    out << "#	k	P(k)	P(k)+noise" << std::endl;   //write an header
    //set the precision for the output 
    out.setf(std::ios_base::scientific);
    out.precision(6);
    out.width(9);

    for(ptrdiff_t i=0; i<abs(psnbins); ++i){
      if(psn_modes_root[i] == 0)   //if there are no modes in the bins, print 0
        out << psk[i] << "\t" << 0. << "\t" << 0. << std::endl;
      else{   // if there are modes
        psPK_root[i] = normalisation * psPK_root[i] / psn_modes_root[i];  //do the averange
        out << psk[i] << "\t" << psPK_root[i] - noise << "\t" << psPK_root[i] << std::endl;
      }
    }
    out.close();  //close the file
    delete []psPK_root;    //dealloc the local array
    delete [] psn_modes_root;
  }
}

/*==========================================================================
 * Convert a vector of strings in a char**
 * Parameters
 * ----------
 * vec: vector of strings
 * output
 * ------
 * chars: array of char arrays
 *==========================================================================*/
char** vector2char(std::vector<std::string> vec){
  size_t svec = vec.size();   //get the number of elements in the vector
  char **chars = new char*[svec];   //allocate the array of char pointers
  for(size_t i=0; i<svec; ++i){     //loop over the vector elements
    chars[i] = new char[vec[i].size()+1];   //allocate enought space to copy each string
    std::strcpy( chars[i], vec[i].c_str() );  //copy the string to the char
  }
  return(chars);  //returns the array of chars
}
