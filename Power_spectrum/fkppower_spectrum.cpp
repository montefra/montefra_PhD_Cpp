/*==========================================================================*/
/* Author: Francesco Montesano, montefra@mpe.mpg.de                         */
/* 06/02/2012,    version: 0.01                                             */
/*                                                                          */
/* compute the averaged power spectrum from a simulation box                */
/*==========================================================================*/

//define that this code compute the power spectrum from a box used in 'power_spectrum.h' if 
//some for some code specific declaration
#define MAIN_CPP   

#include "fkppower_spectrum.h"

int main(int argc, char* argv[])
{

  /*==========================================================================*/
  /* Initialise MPI                                                           */
  /*==========================================================================*/
  int myrank, root=0, n_proc;    //rank of the processor, root processor 
  // and number of processes in the group of the communicator
  MPI_Comm com=MPI_COMM_WORLD;   //MPI communicator

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(com, &myrank);
  MPI_Comm_size(com, &n_proc);

  try{    // Catch ecceptions
    /*==========================================================================*/
    /* command line parsing                                                     */
    /*==========================================================================*/
    // Define the command line object, and insert a message
    // that describes the program. The second argument is the 
    // delimiter (usually space) and the last one is the version number. 
    TCLAP::CmdLine cmd(get_description(), ' ', get_version());
    //very very ugly way to include a file containing the command line options
#include "parse_options.h"
    cmd.add( outfile );  //output file name. Add to the command line
    cmd.add( infiles );   //input file names. Add to the command line
    cmd.add( verbose );  //verbose mode
    cmd.add( pw );   //pw for the FKP weight
    cmd.add( kmax );   //maximum k
    cmd.add( kmin );   //minimum k
    cmd.add( nbins); //number of bins. Add to the command line
    cmd.add( ncells ); //number of cells. Add to the command line
    cmd.add( bmax );  //maximum coordinate of the fft box
    cmd.add( bmin );  //maximum coordinate of the fft boxsize of the fft box
    cmd.add( corr );   //MAS correction
    cmd.add( mas );   //mass assigment scheme
    cmd.add( epsilon );  //add the switch to ln(1+delta) transform
    cmd.parse(argc, argv);  //parse the command line

    /*==========================================================================*/
    /* Check the input and output files                                         */
    /* Check the number of output files and decided what is wanted              */
    /* 1: window function; 2: power spectrum; 4: cross power spectru            */
    /*==========================================================================*/
    std::vector<std::string> vinfiles = infiles.getValue();
    size_t ninfiles = vinfiles.size();   //number of input files

    if( ninfiles != 1 && ninfiles != 2 && ninfiles != 4 ){
      if(myrank == root){ 
        std::cerr << "The number of input files is incorrect. ";
	std::cerr << "1, 2 or 4 files are accepted to compute, ";
	std::cerr << "respectivelly, the window function, the ";
	std::cerr << "power spectrum or the cross power spectrum" << std::endl;
      }
      MPI_Barrier(com);
      MPI_Finalize();
      exit(3);
    }

    for(size_t i=0; i< ninfiles; ++i)  //check the input files
      if(fileexists(vinfiles[i]) == false){
        if(myrank == root)
          std::cerr << "Input file " +vinfiles[i]+ 
            " does not exists. Cleanup and exit" << std::endl;
        MPI_Barrier(com);
        MPI_Finalize();
        exit(4);
      } 
    if(fileexists(outfile.getValue()) == true){ //check the output file
      if(myrank == root)
        std::cerr << "Output file " +outfile.getValue()+ " does already exists. " +
          "Change or rename it. Cleanup and exit" << std::endl;
      MPI_Barrier(com);
      MPI_Finalize();
      exit(5);
    }   


    /*==========================================================================*/
    /* initialise the class for the FFT and read and/or write the wisdom        */
    /*==========================================================================*/
    ptrdiff_t local_n0, local_0_start; //first cell and number of cells in the
                                         //x direction per processor
    if(verbose.getValue() && myrank == root) std::cout << "Initialising the grid" << std::endl;

    ps_r2c_c2r_mpi_inplace grid2(ncells.getValue(), com, 
	&local_n0, &local_0_start); //constructor of the second grid that will be FFT
    if( ninfiles < 4 ) //if the window function or the auto power spectrum required 
      grid2.free_grid();  //free the space of grid2

    ps_r2c_c2r_mpi_inplace grid1(ncells.getValue(), com, 
	&local_n0, &local_0_start); //constructor of the first grid that will be FFT
    if( ninfiles == 1 ) //if the window function required 
      grid1.free_grid();  //free the space of grid1

    ps_r2c_c2r_mpi_inplace gridran(ncells.getValue(), com, 
	&local_n0, &local_0_start); //constructor for the grid containing the random

    //FFTW wisdom file name
    std::string wisdom_fname = "r2c_c2r_inplace_mpi."+to_string(ncells.getValue())+".wis";
    //try to read the fftw wisdom
    int err = gridran.read_wisdom(wisdom_fname, myrank, com);
    if(verbose.getValue() && myrank == root){ 
      if(err == 1) std::cout << "Wisdom read successfully" << std::endl;
      else std::cout << "Wisdom not read or not found" << std::endl;
    }

    gridran.create_plan(FFTW_FORWARD);    
    if( ninfiles != 1 ) //if the power spectra wanted
      grid1.create_plan(FFTW_FORWARD);    //create the real to complex plan
    if( ninfiles == 4 )   //if the cross power spectrum wanted
      grid1.create_plan(FFTW_FORWARD);    //create the real to complex plan
    if(verbose.getValue() && myrank == root)
      std::cout << "Plan created" << std::endl;

    //try to write the fftw wisdom
    err = gridran.write_wisdom(wisdom_fname, myrank, com);
    if(verbose.getValue() && myrank == root){ 
      if(err == 1) std::cout << "Wisdom written successfully" << std::endl;
      else std::cout << "Wisdom not written" << std::endl;
    }
    //set to 0 all the elements of the grid
    gridran.initialise();
    if( ninfiles != 1 ) //if the power spectra wanted
      grid1.initialise();
    if( ninfiles == 4 )   //if the cross power spectrum wanted
      grid2.initialise();
    
    /*==========================================================================*/
    /* set the values of k associated to the grid and set mas                   */
    /*==========================================================================*/
    //set the values of k for the FFT grid and get the size of the cell
    double cell_size = gridran.set_k( bmax.getValue() - bmin.getValue() );
    if( ninfiles != 1 ) //if the power spectra wanted
      cell_size = grid1.set_k( bmax.getValue() - bmin.getValue() );
    if( ninfiles == 4 )   //if the cross power spectrum wanted
      cell_size = grid2.set_k( bmax.getValue() - bmin.getValue() );
    
    //set the mas
    gridran.set_MAS(mas.getValue());
    if( ninfiles != 1 ) //if the power spectra wanted
      grid1.set_MAS(mas.getValue());
    if( ninfiles == 4 )   //if the cross power spectrum wanted
      grid2.set_MAS(mas.getValue());

    /*==========================================================================*/
    /* read the input files                                                     */
    /*==========================================================================*/
    //read the file name and fill the grid(s)
    double *sumsran, *sumscat1, *sumscat2;  //arrays containing the sums done when reading the files
    double alpha1, alpha2;   //arrays containing the values of alpha=sum(w_g)/sum(w_r)
    if( ninfiles == 1 ){   //if the window function is required
      sumsran = read_file(vinfiles[0], gridran, cell_size, bmin.getValue(), pw.getValue());
      
      alpha1 = 1.;  //no alpha in this case
    }
    else{  //if the power spectra are required
      sumscat1 = read_file(vinfiles[0], grid1, cell_size, bmin.getValue(), pw.getValue());
      sumsran = read_file(vinfiles[1], gridran, cell_size, bmin.getValue(), pw.getValue());
      //convert the galaxy and random density fields into F(r)
      alpha1 = grid1.to_Fr(gridran.rgrid, sumscat1[0]/sumsran[0]);  

      if( ninfiles == 4 ){  // if hte cross power spectrum wanted
	gridran.initialise();  //reset the random grid to zero
	sumscat2 = read_file(vinfiles[2], grid2, cell_size, bmin.getValue(), pw.getValue());
	sumsran = read_file(vinfiles[3], gridran, cell_size, bmin.getValue(), pw.getValue());
	//convert the galaxy and random density fields into F(r)
	alpha2 = grid2.to_Fr(gridran.rgrid, sumscat2[0]/sumsran[0]);  
      }
      //free the random grid: non used anymore if the power spectra are required
      gridran.free_grid();
    }

    if(verbose.getValue() && myrank == root){
      std::cout << "Catalogue(s) read and grid(s) filled. ";
      if(ninfiles != 1) std::cout << "F(r) density field(s) created." ;
      std::cout << std::endl;
    }

    /*==========================================================================*/
    /* apply epsilon transform if required                                      */
    /*==========================================================================*/

    //if required do a logarithmic transform
    if(epsilon.isSet() == true){
      if( ninfiles == 1 ){   //if the window function is wanted
	if(epsilon.getValue() < 0) gridran.ln1delta();
	else gridran.ln1delta(epsilon.getValue());
      }
      else{ //if the power spectra are wanted
	if(epsilon.getValue() < 0) grid1.ln1delta();
	else grid1.ln1delta(epsilon.getValue());
	if( ninfiles == 4){ //if cross power spectrum
	  if(epsilon.getValue() < 0) grid2.ln1delta();
	  else grid2.ln1delta(epsilon.getValue());
	}
      }
    }

    /*==========================================================================*/
    /* execute the FFT transform                                                */
    /*==========================================================================*/
    if( ninfiles == 1 )   //if the window function is wanted
      gridran.execute_plan();
    else{ //if the power spectra are wanted
      grid1.execute_plan();
      if( ninfiles == 4) //if cross power spectrum
	grid2.execute_plan();
    }
    if(verbose.getValue() && myrank == root)
      std::cout << "FFT done" << std::endl;

    /*==========================================================================*/
    /* set up all the things needed for the computation of the power spectrum   */
    /* done only on grid1, as grid 2 will be 'hosted' if the cross power is     */
    /* required                                                                 */
    /*==========================================================================*/
    if( ninfiles == 1 ){   //if the window function is wanted
      gridran.set_MAS_correction(corr.getValue());   //set the correction for the MAS
      //set the k bins for the output power spectrum
      gridran.set_psk(kmax.getValue(), kmin.getValue(), nbins.getValue());
    }
    else{ //if the power spectra are wanted
      grid1.set_MAS_correction(corr.getValue());   //set the correction for the MAS
      //set the k bins for the output power spectrum
      grid1.set_psk(kmax.getValue(), kmin.getValue(), nbins.getValue());
    }
    
    /*==========================================================================*/
    /* compute normalisation and shot noise                                     */
    /*==========================================================================*/

    double N2, noise;  //inverse normalisation sqare and shot noise
    if( ninfiles == 1 ){  //window function
      N2 = 1./(alpha1*sumsran[1]);
      noise = sumsran[2]*N2;
    }
    else if( ninfiles == 2 ){
      N2 = 1./(alpha1*sumsran[1]);
      //noise = (sumscat1[2] + alpha1*alpha1*sumsran[2]) * N2;
      noise = (alpha1+1.) * alpha1*sumsran[2] * N2;
    }
    else{
      std::cerr << "Implement normalisation and shot noise for the cross power spectra" 
        << std::endl;
      MPI_Abort(com, 6);
    }

    /*==========================================================================*/
    /* sum the amplitude and number of modes in spherical shells                */
    /* set in 'set_pks'                                                         */
    /*==========================================================================*/
    if( ninfiles == 1 )   //if the window function is wanted
      gridran.sum_modes2_sph();   //window function
    else if( ninfiles == 2 )   //if the power spectra are wanted
      grid1.sum_modes2_sph();   //auto power spectrum
    else
      grid1.sum_modes2_sph(grid2.rgrid);   //cross power spectrum
//  apply normalisation and noise before correcting
//    if( ninfiles == 1 )   //if the window function is wanted
//      gridran.sum_modes2_sph(N2, noise);   //window function
//    else if( ninfiles == 2 )   //if the power spectra are wanted
//      grid1.sum_modes2_sph(N2, noise);   //auto power spectrum
//    else
//      grid1.sum_modes2_sph(grid2.rgrid);   //cross power spectrum

    /*==========================================================================*/
    /* print the output file                                                    */
    /*==========================================================================*/
    if(verbose.getValue() && myrank == root)
      std::cout << "Saving the power spectrum to file: " << outfile.getValue() << std::endl;
    if( ninfiles == 1 )   //if the window function is wanted
      gridran.savePK(outfile.getValue(), N2, noise, myrank, root, com);
    else
      grid1.savePK(outfile.getValue(), N2, noise, myrank, root, com);
    //the shot noise and normalisation already applied
    //grid1.savePK(outfile.getValue(), 1., 0., myrank, root, com);

  } 
  catch (TCLAP::ArgException &e){  // catch any exceptions
    if(myrank == root) 
      std::cerr << "Error: '" << e.error() << "' for arg '" << e.argId() << "'." <<  std::endl; 
    MPI_Finalize();
    exit(10);
  }

  /*==========================================================================*/
  /* Finalize MPI                                                             */
  /*==========================================================================*/
  MPI_Finalize();

  exit(0);
}

