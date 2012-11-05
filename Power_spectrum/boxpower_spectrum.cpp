/*==========================================================================*/
/* Author: Francesco Montesano, montefra@mpe.mpg.de                         */
/* 06/02/2012,    version: 0.01                                             */
/*                                                                          */
/* compute the averaged power spectrum from a simulation box                */
/*==========================================================================*/

//define that this code compute the power spectrum from a box used in 'power_spectrum.h' if 
//some for some code specific declaration
#define MAIN_CPP

#include "boxpower_spectrum.h"

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
    cmd.add( infile );   //input file name. Add to the command line
    cmd.add( outfile );  //output file name. Add to the command line
    cmd.add( verbose );  //verbose mode
    cmd.add( kmax );   //maximum k
    cmd.add( kmin );   //minimum k
    cmd.add( nbins); //number of bins. Add to the command line
    cmd.add( ncells ); //number of cells. Add to the command line
    cmd.add( bsize );  //size of the fft box
    cmd.add( corr );   //MAS correction
    cmd.add( mas );   //mass assigment scheme
    cmd.add( dm );   //add the possibility to add dark matter
    cmd.add( zdist );   //add redshift space distortion switch
    cmd.add( epsilon );  //add the switch to ln(1+delta) transform
    cmd.add( mlim );   //add the mass limits
    cmd.parse(argc, argv);  //parse the command line

    /*==========================================================================*/
    /* Check the input and output files                                         */
    /*==========================================================================*/

    if( dm.isSet() == false){  //if the input file is a halo catalogue
      if(fileexists(infile.getValue()) == false){
	std::cerr << "Input file " +infile.getValue()+ " does not exists" << std::endl;
	MPI_Abort( com, 2);
      }
    }
    else{
      for( int i=0; i<dm.getValue(); ++i ){
	if(fileexists(infile.getValue()+"."+to_string(i)) == false){
	  std::cerr << "Input file " +infile.getValue()+"."+to_string(i)+ 
	    " does not exists" << std::endl;
	  MPI_Abort( com, 3);
	}
      }
    }
    if(fileexists(outfile.getValue()) == true){
      std::cerr << "Output file " +outfile.getValue()+ " does already exists. " +
	"Change or rename it" << std::endl;
      MPI_Abort( com, 4);
    } 

    /*==========================================================================*/
    /* Check the mass bins                                                      */
    /*==========================================================================*/
    std::vector<double> mlims;
    size_t n_mlim=0;  //size of the mass limits
    if( mlim.isSet() == true){   //if the option is called
      mlims = mlim.getValue();   //save the mass limits in a local variable
      n_mlim = mlims.size();    //save the size of the mass limits
      if( n_mlim != 1 && n_mlim != 2 && n_mlim != 4 ){   //if the number of bins is not 1,2 or 4
        if(myrank == root){
	  std::cerr << "The option 'mass-limits' can be called 1,2 or 4 times. ";
	  std::cerr << "Mass limits ignored" << std::endl;
	}
	n_mlim = 0;
      }
    }

    /*==========================================================================*/
    /* initialise the class for the FFT and read and/or read the wisdom         */
    /*==========================================================================*/
    ptrdiff_t local_n0, local_0_start; //first cell and number of cells in the
                                         //x direction per processor
    if(verbose.getValue() && myrank == root) std::cout << "Initialising the grid" << std::endl;
    //grid 2 initilized first.
    ps_r2c_c2r_mpi_inplace grid2(ncells.getValue(), com, 
	&local_n0, &local_0_start); //constructor
    if(n_mlim != 4)        //if the cross power spectrum is not required free the space in memory
      grid2.free_grid();   //Done to be able to compile as it is not possible to initialise something in the if
    ps_r2c_c2r_mpi_inplace grid1(ncells.getValue(), com, 
	&local_n0, &local_0_start); //constructor

    //FFTW wisdom file name
    std::string wisdom_fname = "r2c_c2r_inplace_mpi."+to_string(ncells.getValue())+".wis";
    //try to read the fftw wisdom
    int err = grid1.read_wisdom(wisdom_fname, myrank, com);
    if(verbose.getValue() && myrank == root){ 
      if(err == 1) std::cout << "Wisdom read successfully" << std::endl;
      else std::cout << "Wisdom not read or not found" << std::endl;
    }

    if(epsilon.isSet() == false)
      grid1.create_plan(FFTW_FORWARD);    //create the real to complex plan
    else
      grid1.create_plan(FFTW_BOTH);    //create the real to complex and complex to real plan
    if(n_mlim == 4){
      if(epsilon.isSet() == false)
	grid2.create_plan(FFTW_FORWARD);    //create the real to complex plan
      else
	grid2.create_plan(FFTW_BOTH);    //create the real to complex and complex to real plan
    }
    if(verbose.getValue() && myrank == root)
      std::cout << "Plan created" << std::endl;

    //try to write the fftw wisdom
    err = grid1.write_wisdom(wisdom_fname, myrank, com);
    if(verbose.getValue() && myrank == root){ 
      if(err == 1) std::cout << "Wisdom written successfully" << std::endl;
      else std::cout << "Wisdom not written" << std::endl;
    }
    //set to 0 all the elements of the grid
    grid1.initialise();
    if(n_mlim == 4) grid2.initialise();
    
    /*==========================================================================*/
    /* set the values of k associated to the grid and set mas                   */
    /*==========================================================================*/
    //set the values of k for the FFT grid and get the size of the cell
    double cell_size = grid1.set_k( bsize.getValue() );
    if(n_mlim == 4 ) cell_size = grid2.set_k( bsize.getValue() );
    
    //set the mas
    grid1.set_MAS(mas.getValue());
    if(n_mlim == 4 ) grid2.set_MAS(mas.getValue());

    /*==========================================================================*/
    /* set up all the things needed for the computation of the power spectrum   */
    /* done only on grid1, as grid 2 will be 'hosted' if the cross power is     */
    /* required                                                                 */
    /*==========================================================================*/
    grid1.set_MAS_correction(corr.getValue());   //set the correction for the MAS
    //set the k bins for the output power spectrum
    grid1.set_psk(kmax.getValue(), kmin.getValue(), nbins.getValue());
    

    /*==========================================================================*/
    /* read the input files                                                     */
    /*==========================================================================*/
    //read the file name and fill the grid(s)
    ptrdiff_t *n_part = new ptrdiff_t[3];  //number of particles fed to the grid. All three 
      //elements used only when the cross power spectrum is wanted
    if( zdist.isSet() == false ){   //no redshift space distortions
      if( dm.isSet() == true )    //dark matter
	n_part[0] = read_dm(infile.getValue(), dm.getValue(), grid1, cell_size);
      else if( n_mlim == 0 )   //no mass limits
	n_part[0] = read_haloes(infile.getValue(), grid1, cell_size);
      else if( n_mlim == 1 )   //mass lower limin
	n_part[0] = read_haloes(infile.getValue(), grid1, mlims[0], cell_size);
      else if( n_mlim == 2 )    //mass bin
	n_part[0] = read_haloes(infile.getValue(), grid1, mlims[0], mlims[1], cell_size);
      else    //two mass bins for the cross power spectrum
        n_part = read_haloes(infile.getValue(), grid1, grid2, mlims, cell_size);
    }
    else{   //redshift space distortions
      if( dm.isSet() == true )   //dark matter
	n_part[0] = read_dm(infile.getValue(), dm.getValue(), grid1, cell_size, 
	    zdist.getValue());
      else if( n_mlim == 0 )    //no mass limits
	n_part[0] = read_haloes(infile.getValue(), grid1, cell_size, 
	    zdist.getValue());
      else if( n_mlim == 1 )   //mass lower limit
	n_part[0] = read_haloes(infile.getValue(), grid1, mlims[0], cell_size,
	    zdist.getValue());
      else if( n_mlim == 2 )   //mass bin
	n_part[0] = read_haloes(infile.getValue(), grid1, mlims[0], mlims[1], cell_size,
	    zdist.getValue());
      else    //two mass bins for the cross power spectrum
	n_part = read_haloes(infile.getValue(), grid1, grid2, mlims, cell_size, 
	    zdist.getValue());
    }
    if(verbose.getValue() && myrank == root)
      std::cout << "Catalogue(s) read and grid(s) filled" << std::endl;

    /*==========================================================================*/
    /* transform numbers of object to overdensity of objects and do ln transform*/
    /*==========================================================================*/
    //transform the number per object to overdensity and save the mean density
    double dmean[2];   //mean density for the (two) grid(s)
    dmean[0] = grid1.to_overdensity(n_part[0]) / pow(cell_size, 3.); 
    if(n_mlim == 4 ) dmean[1] = grid2.to_overdensity(n_part[1]) / pow(cell_size, 3.);

    /*==========================================================================*/
    /* execute the FFT transform                                                */
    /*==========================================================================*/
    grid1.execute_plan(FFTW_FORWARD);
    if(n_mlim == 4 ) grid2.execute_plan(FFTW_FORWARD);

    /*==========================================================================*/
    /* If the logarithmic transform required correct for the MAS, FT back, do   */
    /* the tranfrom, re-execute the forward plan.                               */
    /*==========================================================================*/
    if(epsilon.isSet() == true){
      //correct for the effects on the mas assigment scheme
      grid1.correct_modes();
      if(n_mlim == 4 ) grid2.correct_modes();

      //go back to real space
      grid1.execute_plan(FFTW_BACKWARD);
      if(n_mlim == 4 ) grid2.execute_plan(FFTW_BACKWARD);

      //rescale the amplitude of the modes after the double transform
      grid1.normalize_double_transform();
      if(n_mlim == 4 ) grid2.normalize_double_transform();

      //do the log tranform
      if(epsilon.getValue() < 0){
	grid1.ln1delta();
	if(n_mlim == 4 ) grid2.ln1delta();
      }
      else{
	grid1.ln1delta(epsilon.getValue());
	if(n_mlim == 4 ) grid2.ln1delta(epsilon.getValue());
      }

      //re-execute the forward transform
      grid1.execute_plan(FFTW_FORWARD);
      if(n_mlim == 4 ) grid2.execute_plan(FFTW_FORWARD);

      //The MAS has already been corrected. Set to no correction
      grid1.set_MAS_correction( 0 );   
    }

    if(verbose.getValue() && myrank == root)
      std::cout << "FFT done" << std::endl;

    /*==========================================================================*/
    /* sum the amplitude and number of modes in spherical shells                */
    /* set in 'set_pks'                                                         */
    /*==========================================================================*/
    if(n_mlim != 4 ) grid1.sum_modes2_sph();   //autopower spectrum
    else grid1.sum_modes2_sph(grid2.rgrid);   //cross power spectrum

    /*==========================================================================*/
    /* compute normalisation and shot noise                                     */
    /*==========================================================================*/
    // normalization is V/N^2 = v_{cell}/N_{cell}
    double normalisation = pow( cell_size/ncells.getValue(), 3.);
    //shot noise 
    // 1/n_mean, if power spectrum
    // n_{common}/n1*n2 = VN_{common}/N1_{particle}*N2_{particle} if cross power spectrum
    /* dimensional n_mean */
    double noise = 1. / dmean[0];
    if( n_mlim == 4 ) noise = n_part[2] * pow(bsize.getValue(), 3.) / (n_part[0]*n_part[1]);


    /*==========================================================================*/
    /* print the output file                                                    */
    /*==========================================================================*/
    if(verbose.getValue() && myrank == root)
      std::cout << "Saving the power spectrum to file: " << outfile.getValue() << std::endl;
    grid1.savePK(outfile.getValue(), normalisation, noise, myrank, root, com);

  } 
  catch (TCLAP::ArgException &e){  // catch any exceptions
    if(myrank == root) 
      std::cerr << "Error: '" << e.error() << "' for arg '" << e.argId() << "'." <<  std::endl; 
    MPI_Finalize();
    exit(3);
  }

  /*==========================================================================*/
  /* Finalize MPI                                                             */
  /*==========================================================================*/

  MPI_Finalize();

  exit(0);
}

