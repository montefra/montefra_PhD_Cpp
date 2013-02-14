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
    cmd.add(outfile);  //output file name. Add to the command line
    cmd.add(infiles);   //input file names. Add to the command line
    cmd.add(winonly);   //the input files are for the window function only
    cmd.add(crosspk);   //cross power spectrum
    cmd.add(verbose);  //verbose mode
    cmd.add(pw);   //pw for the FKP weight
    cmd.add(kmax);   //maximum k
    cmd.add(kmin);   //minimum k
    cmd.add(nbins); //number of bins. Add to the command line
    cmd.add(ncells); //number of cells. Add to the command line
    cmd.add(bmax);  //maximum coordinate of the fft box
    cmd.add(bmin);  //maximum coordinate of the fft boxsize of the fft box
    cmd.add(corr);   //MAS correction
    cmd.add(mas);   //mass assigment scheme
    cmd.add(whichnoise);   //noise computation method
    cmd.add(ignorew);   //ignore some of the weights from the input
    cmd.add(zrange);  //set zrange
    cmd.add(repeatw);  //repeat particle w times instead of assigning weigh w
    cmd.add(epsilon);  //add the switch to ln(1+delta) transform
    cmd.parse(argc, argv);  //parse the command line

    /*==========================================================================*/
    /* Check the input and output files                                         */
    /*==========================================================================*/
    std::vector<std::string> vinfiles = infiles.getValue();
    size_t ninfiles = vinfiles.size();   //number of input files

    // if the cross power spectrum wanted. PARTIALLY IMPLEMENTED
    // CHECK NORMALISATION, SHOT NOISE AND WINDOWFUNCTION
    if(crosspk.getValue() == true){
      check_crosspk(vinfiles, myrank, root, com);
    }

    // if the standard power spectrum required check that the number of files is even
    if(winonly.getValue() == false && ninfiles%2 != 0){
      if(myrank == root)
        std::cerr << "An even number of files is needed" << std::endl;
      MPI_Barrier(com);
      MPI_Finalize();
      exit(2);
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
    /* check the redshift range                                                 */
    /*==========================================================================*/
    std::vector<double> vzrange(2, -1); //default valued for the redshift range
    if(zrange.isSet() == true){
      size_t nzrange = zrange.getValue().size(); //get the number of given values
      if(nzrange<2){ //
        if(myrank == root)
          std::cerr << "Two values required by '--zrange' option" << std::endl;
        MPI_Barrier(com);
        MPI_Finalize();
        exit(6);
      }
      else if(nzrange>2){
        if(myrank == root)
          std::cerr << "Skipping values after the second in option '--zrange'"
            << std::endl;
        for(int i=0; i<2; ++i) vzrange[i] = zrange.getValue()[i];
      }
      else vzrange = zrange.getValue();
      if(vzrange[1]<=0){
        if(myrank == root)
          std::cerr << "A negative upper value for the redshift range " <<
            "does not make sense" << std::endl;
        MPI_Barrier(com);
        MPI_Finalize();
        exit(7);
      }
    }

    /*==========================================================================*/
    /* initialise the class for the FFT and read and/or write the wisdom        */
    /*==========================================================================*/
    ptrdiff_t local_n0, local_0_start; //first cell and number of cells in the
                                         //x direction per processor
    if(verbose.getValue() && myrank == root) std::cout << "Initialising the grid" << std::endl;

    //initialise the FFT objects
    ps_r2c_c2r_mpi_inplace *grid1;
    ps_r2c_c2r_mpi_inplace *grid2;    
    ps_r2c_c2r_mpi_inplace *gridran;

    // the one grid is used in any case
    grid1 = new ps_r2c_c2r_mpi_inplace(ncells.getValue(), com, &local_n0,
        &local_0_start); 
    if(winonly.getValue() == false) //power spectrum
      gridran = new ps_r2c_c2r_mpi_inplace(ncells.getValue(), com, &local_n0,
          &local_0_start); 
    if(crosspk.getValue() == true) //cross power spectrum
      grid2 = new ps_r2c_c2r_mpi_inplace(ncells.getValue(), com, &local_n0,
          &local_0_start); 

    //FFTW wisdom file name
    std::string wisdom_fname = "r2c_c2r_inplace_mpi."+to_string(ncells.getValue())+".wis";
    //try to read the fftw wisdom
    int err = grid1->read_wisdom(wisdom_fname, myrank, com);
    if(verbose.getValue() && myrank == root){ 
      if(err == 1) std::cout << "Wisdom read successfully" << std::endl;
      else std::cout << "Wisdom not read or not found" << std::endl;
    }

    //create the real to complex plan
      grid1->create_plan(FFTW_FORWARD);
    if(winonly.getValue() == false) //power spectrum
      gridran->create_plan(FFTW_FORWARD);    
    if(crosspk.getValue() == true) //cross power spectrum
      grid2->create_plan(FFTW_FORWARD);
    if(verbose.getValue() && myrank == root)
      std::cout << "Plan created" << std::endl;

    //try to write the fftw wisdom
    err = grid1->write_wisdom(wisdom_fname, myrank, com);
    if(verbose.getValue() && myrank == root){ 
      if(err == 1) std::cout << "Wisdom written successfully" << std::endl;
      else std::cout << "Wisdom not written" << std::endl;
    }
    //set to 0 all the elements of the grid
    grid1->initialise();
    if(winonly.getValue() == false) //power spectrum
      gridran->initialise();
    if(crosspk.getValue() == true) //cross power spectrum
      grid2->initialise();
    
    /*==========================================================================*/
    /* set the values of k associated to the grid and set mas                   */
    /*==========================================================================*/
    //set the values of k for the FFT grid and get the size of the cell
    double cell_size = grid1->set_k( bmax.getValue() - bmin.getValue() );
    if(winonly.getValue() == false) //power spectrum
      cell_size = gridran->set_k( bmax.getValue() - bmin.getValue() );
    if(crosspk.getValue() == true) //cross power spectrum
      cell_size = grid2->set_k( bmax.getValue() - bmin.getValue() );
    
    //set the mas
    grid1->set_MAS(mas.getValue());
    if(winonly.getValue() == false) //power spectrum
      gridran->set_MAS(mas.getValue());
    if(crosspk.getValue() == true) //cross power spectrum
      grid2->set_MAS(mas.getValue());

    /*==========================================================================*/
    /* read the input files, compute normalisation and shot noise and create    */
    /* density fields to fourier transform                                      */
    /*==========================================================================*/

    // arrays containing sum(w), sum(n*w^2), sum(w^2)
    size_t dimsums = 3;  //dimensions of the sums
    double *sumstemp=new double[dimsums];  //temporary container for the sums
    double *sumscat1=new double[dimsums];   //sums from the catalogues
    double *sumscat2=new double[dimsums];   //sums from the catalogues
    double *sumsran=new double[dimsums];   //sums from the random catalogues
    for(size_t i=0; i<dimsums; ++i){ //initialize to 0
      sumscat1[i] = 0.;
      sumscat2[i] = 0.;
      sumsran[i] = 0.;
    }
    
    // values of alpha=sum(w_cat)/sum(w_ran) or 1 for window functions
    double alpha1=1., alpha2=1.;
    //inverse normalisation square and shot noise 
    double N2, noise; 

    //initilise the class that reads the file(s)
    readfiles rf(cell_size, bmin.getValue(), pw.getValue());

    //window function only
    if(winonly.getValue() == true){
      for(size_t i=0; i<ninfiles; ++i){
        sumstemp = rf.read_file(vinfiles[i], *grid1, vzrange, -1, false);
        for(size_t j=0; j<dimsums; ++j) sumscat1[j] += sumstemp[j];
      }
      N2 = 1./sumscat1[1];
      noise = sumscat1[2]*N2;
    }
    else{
      if(crosspk.getValue()==false){  // power spectrum
        for(size_t i=0; i<ninfiles/2; ++i){
          sumstemp = rf.read_file(vinfiles[2*i], *grid1, vzrange,  //read the catalogue
              ignorew.getValue(), repeatw.getValue());
          for(size_t j=0; j<dimsums; ++j) sumscat1[j] += sumstemp[j];
          sumstemp = rf.read_file(vinfiles[2*i+1], *gridran, vzrange, -1,    //read the random
              false);
          for(size_t j=0; j<dimsums; ++j) sumsran[j] += sumstemp[j];
        }
        alpha1 = sumscat1[0]/sumsran[0];
        N2 = 1./(alpha1*sumsran[1]);
        //N2 = 1./sumscat1[1];
        if(whichnoise.getValue() ==  1)
          noise = (alpha1+1.) * alpha1*sumsran[2] * N2;
        else
          noise = (sumscat1[2] + alpha1*alpha1*sumsran[2]) * N2;
        //convert the galaxy and random density fields into F(r)
        grid1->to_Fr(gridran->rgrid, alpha1);
        //grid1->to_Fr(gridran->rgrid, alpha1, sqrt(N2));
      }
      else{  // cross power spectrum
        std::cerr << "NOT IMPLEMENTED" << std::endl;
      }
      gridran->free_grid();
    }
    delete [] sumstemp; //not used anymore

    if(verbose.getValue() && myrank == root){
      std::cout << "Catalogue(s) read and grid(s) filled. ";
      if(winonly.getValue() == false) std::cout << "F(r) density field(s) created." ;
      std::cout << std::endl;
    }

    /*==========================================================================*/
    /* apply epsilon transform if required                                      */
    /*==========================================================================*/

    //if required do a logarithmic transform
    if(epsilon.isSet() == true){
      if(epsilon.getValue() < 0) grid1->ln1delta();
      else grid1->ln1delta(epsilon.getValue());

      if(crosspk.getValue() == true){ // cross power spectrum
        std::cerr << "NOT IMPLEMENTED" << std::endl;
      }
    }

    /*==========================================================================*/
    /* execute the FFT transform                                                */
    /*==========================================================================*/

    grid1->execute_plan();
    if(crosspk.getValue() == true){ // cross power spectrum
      std::cerr << "NOT IMPLEMENTED" << std::endl;
    }
    if(verbose.getValue() && myrank == root)
      std::cout << "FFT done" << std::endl;

    /*==========================================================================*/
    /* set up all the things needed for the computation of the power spectrum   */
    /* done only on grid1, as grid 2 will be 'hosted' if the cross power is     */
    /* required                                                                 */
    /*==========================================================================*/
    //set the correction for the MAS
    grid1->set_MAS_correction(corr.getValue());   
    //set the k bins for the output power spectrum
    grid1->set_psk(kmax.getValue(), kmin.getValue(), nbins.getValue());
    

    /*==========================================================================*/
    /* sum the amplitude and number of modes in spherical shells                */
    /* set in 'set_pks'                                                         */
    /*==========================================================================*/
    grid1->sum_modes2_sph();   //auto power spectrum
    if(crosspk.getValue() == true){ // cross power spectrum
      std::cerr << "NOT IMPLEMENTED" << std::endl;
      grid1->sum_modes2_sph(grid2->rgrid);   //cross power spectrum
      grid2->free_grid();  //cleanup the grid
    }

    /*==========================================================================*/
    /* print the output file                                                    */
    /*==========================================================================*/
    if(verbose.getValue() && myrank == root)
      std::cout << "Saving the power spectrum to file: " << outfile.getValue() << std::endl;
    // create the header to the output file
    std::string header;
    header = "# \t sum(w) \t sum(w^2n(z)) \t sum(w^2)\n";
    header += "#data";
    for(size_t i=0; i<dimsums; ++i) header += "\t"+to_string(sumscat1[i], 4);
    header += "\n";
    if(crosspk.getValue() == true){
      std::cerr << "NOT IMPLEMENTED" << std::endl;
      header += "#data2";
      for(size_t i=0; i<dimsums; ++i) header += "\t"+to_string(sumscat1[i], 4);
      header += "\n";
    }
    header += "#random";
    for(size_t i=0; i<dimsums; ++i) header += "\t"+to_string(sumsran[i], 4);

    grid1->savePK(outfile.getValue(), N2, noise, header, myrank, root, com);
    //the shot noise and normalisation already applied
    //grid1->savePK(outfile.getValue(), 1., noise, myrank, root, com);

    delete [] sumscat1; //not used anymore
    delete [] sumscat2; //not used anymore
    delete [] sumsran; //not used anymore
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

