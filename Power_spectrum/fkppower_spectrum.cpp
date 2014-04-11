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
    cmd.add(cross);   //cross power spectrum
    cmd.add(twogrids); //use two grids to compute the density field
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
    //cmd.add(whichnoise);   //noise computation method
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

    // Check the number of the input files
    check_input_files(vinfiles, winonly.getValue(), cross.getValue(),
        twogrids.getValue(), myrank, root, com);

    for(size_t i=0; i< ninfiles; ++i)  //check the input files
      if(fileexists(vinfiles[i]) == false)
        on_error(std::string("Input file '")+vinfiles[i]+std::string("' does not exists"),
            5, myrank, root, com);
    if(fileexists(outfile.getValue()) == true) //check the output file
        on_error("Output file '"+outfile.getValue()+
            "' already exists. Change or rename it", 6, myrank, root, com);

    /*==========================================================================*/
    /* check the redshift range                                                 */
    /*==========================================================================*/
    std::vector<double> vzrange(2, -1); //default valued for the redshift range
    if(zrange.isSet() == true){
      size_t nzrange = zrange.getValue().size(); //get the number of given values
      if(nzrange<2)
        on_error("Two values required by '--zrange' option", 7, myrank, root,
            com);
      else if(nzrange>2){
        if(myrank == root)
          std::cerr << "Skipping values after the second in option '--zrange'"
            << std::endl;
        for(int i=0; i<2; ++i) vzrange[i] = zrange.getValue()[i];
      }
      else vzrange = zrange.getValue();
      if(vzrange[1]<=0)
        on_error("A negative upper value for the redshift range does not make sense", 
            7, myrank, root, com);
    }
    /*==========================================================================*/
    /* check ignorew                                                            */
    /* if ignorew is not given create a vector with one element -1              */
    /* else check how many numbers are given and avoid repetitions              */
    /*==========================================================================*/
    std::vector<int> vignorew; //default value for ignorew
    if(ignorew.isSet() == true){
      std::vector<int> temp = ignorew.getValue();
      for(int i=1; i<3; ++i)
        if(std::find(temp.begin(), temp.end(), i) != temp.end())
          vignorew.push_back(i);
    }
    else vignorew.push_back(-1);
    std::vector<int> vignore_randoms(1, -1); //randoms don't have weights to ignore

    /*==========================================================================*/
    /* initialise the class for the FFT and read and/or write the wisdom        */
    /*==========================================================================*/
    ptrdiff_t local_n0, local_0_start; //first cell and number of cells in the
                                         //x direction per processor
    if(verbose.getValue() && myrank == root) 
      std::cout << "Initialising the grid" << std::endl;

    //initialise the FFT objects
    ps_r2c_c2r_mpi_inplace *grid1;
    ps_r2c_c2r_mpi_inplace *grid2;    
    ps_r2c_c2r_mpi_inplace *gridran;

    // the one grid is used in any case
    grid1 = new ps_r2c_c2r_mpi_inplace(ncells.getValue(), com, &local_n0,
        &local_0_start); 
    if(winonly.getValue()==false)
      gridran = new ps_r2c_c2r_mpi_inplace(ncells.getValue(), com, &local_n0,
          &local_0_start); 
    if(twogrids.getValue() || cross.getValue())
      grid2 = new ps_r2c_c2r_mpi_inplace(ncells.getValue(), com, &local_n0,
          &local_0_start); 

    //FFTW wisdom file name
    std::string wisdom_fname =
      "r2c_c2r_inplace_mpi."+to_string(ncells.getValue())+".wis";
    //try to read the fftw wisdom
    int err = grid1->read_wisdom(wisdom_fname, myrank, com);
    if(verbose.getValue() && myrank == root){ 
      if(err == 1) std::cout << "Wisdom read successfully" << std::endl;
      else std::cout << "Wisdom not read or not found" << std::endl;
    }

    //create the real to complex plan
    grid1->create_plan(FFTW_FORWARD);
    if(winonly.getValue()==false)
      gridran->create_plan(FFTW_FORWARD);    
    if(twogrids.getValue() || cross.getValue())
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
    if(winonly.getValue()==false)
      gridran->initialise();
    if(twogrids.getValue() || cross.getValue())
      grid2->initialise();
    
    /*==========================================================================*/
    /* set the values of k associated to the grid and set mas                   */
    /*==========================================================================*/
    //set the values of k for the FFT grid and get the size of the cell
    //set the mas
    double cell_size = grid1->set_k( bmax.getValue() - bmin.getValue() );
    grid1->set_MAS(mas.getValue());
    if(winonly.getValue()==false){
      cell_size = gridran->set_k( bmax.getValue() - bmin.getValue() );
      gridran->set_MAS(mas.getValue());
    }
    if(twogrids.getValue() || cross.getValue()){
      cell_size = grid2->set_k( bmax.getValue() - bmin.getValue() );
      grid2->set_MAS(mas.getValue());
    }
    
    /*==========================================================================*/
    /* read the input files, compute normalisation and shot noise and create    */
    /* density fields to fourier transform                                      */
    /*==========================================================================*/

    // arrays containing sum(w), sum(n*w^2), sum(w^2)
    size_t dimsums = 3;  //dimensions of the sums
    double *sumscat1=new double[dimsums];   //sums from the catalogues
    double *sumscat2=new double[dimsums];   //sums from the catalogues
    double *sumsran1=new double[dimsums];   //sums from the random catalogues
    double *sumsran2=new double[dimsums];   //sums from the random catalogues
    for(size_t i=0; i<dimsums; ++i){ //initialize to 0
      sumscat1[i] = 0.;
      sumscat2[i] = 0.;
      sumsran1[i] = 0.;
      sumsran2[i] = 0.;
    }
    
    // values of alpha=sum(w_cat)/sum(w_ran) or 1 for window functions
    double alpha1=1., alpha2=1.;
    //normalisation square and shot noise for the two samples
    double N2_1, noise1, N2_2, noise2; 

    //initilise the class that reads the file(s)
    readfiles rf(cell_size, bmin.getValue(), pw.getValue());

    //window function only
    if(winonly.getValue() == true){
      //read anyway the first file
      sumscat1 = rf.read_file(vinfiles[0], *grid1, vzrange, vignore_randoms, false);
      if(ninfiles == 2){  //read the second
        if(twogrids.getValue() || cross.getValue())  //if two grids are wanted 
          sumscat2 = rf.read_file(vinfiles[1], *grid2, vzrange, vignore_randoms, false);
        else{  //just add the second file on the same grid and then add the weights
          sumscat2 = rf.read_file(vinfiles[1], *grid1, vzrange, vignore_randoms, false);
          for(size_t i=0; i<dimsums; ++i) sumscat1[i] += sumscat2[i];
        }
      }
      //compute alpha,  normalisation and shot noise (this is done anyway)
      alpha_N_sh(sumscat1, &alpha1, &N2_1, &noise1);
      *grid1 /= sqrt(N2_1);  //multiply by the normalisation
      if(twogrids.getValue() || cross.getValue()){  //if two grids are wanted 
      // cross window function or from two catalogues
        alpha_N_sh(sumscat2, &alpha2, &N2_2, &noise2);
        *grid2 /= sqrt(N2_2);  //multiply by the normalisation
      }
    }
    else{
      // read anyway the first catalogue and random
      sumscat1 = rf.read_file(vinfiles[0], *grid1, vzrange,
          vignorew, repeatw.getValue());
      sumsran1 = rf.read_file(vinfiles[1], *gridran, vzrange,
          vignore_randoms, false);
      if(ninfiles == 2){
          //compute alpha inverse normalisation and shot noise (this is done anyway)
          alpha_N_sh(sumscat1, sumsran1, &alpha1, &N2_1, &noise1);
          //convert the galaxy and random density fields into F(r)
          grid1->to_Fr(*gridran, alpha1, sqrt(N2_1));
      }
      else if(ninfiles == 4){  //read the second catalogue and random
        if(twogrids.getValue() || cross.getValue()){  //if two grids are wanted 
          //compute alpha, inverse normalisation and shot noise
          alpha_N_sh(sumscat1, sumsran1, &alpha1, &N2_1, &noise1);
          //convert the galaxy and random density fields into F(r)
          grid1->to_Fr(*gridran, alpha1, sqrt(N2_1));
          // cross power spectrum or power spectrum on two grids
          gridran->initialise();  //reset to 0 the random grid
          // read the second catalogue and random
          sumscat2 = rf.read_file(vinfiles[2], *grid2, vzrange,
              vignorew, repeatw.getValue());
          sumsran2 = rf.read_file(vinfiles[3], *gridran, vzrange,
              vignore_randoms, false);
          alpha_N_sh(sumscat2, sumsran2, &alpha2, &N2_2, &noise2);
          //convert the galaxy and random density fields into F(r)
          grid2->to_Fr(*gridran, alpha2, sqrt(N2_2));
        }
        else{  //just add the second file and random on the same grids and then add the weights
          sumscat2 = rf.read_file(vinfiles[2], *grid1, vzrange,
              vignorew, repeatw.getValue());
          sumsran2 = rf.read_file(vinfiles[3], *gridran, vzrange,
              vignore_randoms, false);
          for(size_t i=0; i<dimsums; ++i){
            sumscat1[i] += sumscat2[i];
            sumsran1[i] += sumsran2[i];
          }
          //compute alpha inverse normalisation and shot noise (this is done anyway)
          alpha_N_sh(sumscat1, sumsran1, &alpha1, &N2_1, &noise1);
          //convert the galaxy and random density fields into F(r)
          grid1->to_Fr(*gridran, alpha1, sqrt(N2_1));
        }
      }
      // random grid not used anymore
      gridran->free_grid();
    }
    
    // if two (four) files are provided to compute the window function (power
    // spectrum) but no cross wf (ps) required, sum the first and second grid
    if(cross.getValue() == true || twogrids.getValue() == true){
      double N_tot = N2_1+N2_2;  //square normalisation for the full sample
      noise1 *= N2_1/N_tot;  //rescale also the shot noise
      noise2 *= N2_2/N_tot;
      N_tot = sqrt(N_tot);   //normalisation for the full sample
      *grid1 *= sqrt(N2_1)/N_tot;  //rescale the two grids to match the full sample amplitude
      *grid2 *= sqrt(N2_2)/N_tot;
      if(twogrids.getValue()){
        *grid1 += *grid2;  //sum the two grids if you don't want the cross 
        grid2->free_grid(); // grid 2 not needed anymore
      }
    }

    if(verbose.getValue() && myrank == root){
      std::cout << "Catalogue(s) read and grid(s) filled. ";
      std::cout << "F(r) density field(s) created." ;
      std::cout << std::endl;
    }

    /*==========================================================================*/
    /* apply epsilon transform if required                                      */
    /*==========================================================================*/

    //if required do a logarithmic transform
    if(epsilon.isSet() == true){
      if(epsilon.getValue() < 0) grid1->ln1delta();
      else grid1->ln1delta(epsilon.getValue());
      if(cross.getValue()){
        if(epsilon.getValue() < 0) grid2->ln1delta();
        else grid2->ln1delta(epsilon.getValue());
      }
    }

    /*==========================================================================*/
    /* execute the FFT transform                                                */
    /*==========================================================================*/

    grid1->execute_plan();
    if(cross.getValue() == true) // cross power spectrum
      grid2->execute_plan();
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
    if(cross.getValue() == true){ // cross power spectrum
      grid1->sum_modes2_sph(*grid2);   //cross power spectrum
      grid2->free_grid();  //cleanup the grid not used anymore
    }
    else grid1->sum_modes2_sph();   //auto power spectrum

    /*==========================================================================*/
    /* print the output file                                                    */
    /*==========================================================================*/
    if(verbose.getValue() && myrank == root)
      std::cout << "Saving the power spectrum to file: " << outfile.getValue() << std::endl;
    // create the header to the output file
    std::string header = create_header(sumscat1, sumscat2, sumsran1, sumsran2,
        dimsums, twogrids.getValue(), cross.getValue(), winonly.getValue());

    double noise;  //shot noise to substract from the spherical averaged power spectrum
    if(cross.getValue()==true) noise=0.;  //no shot noise in the cross power spectrum
    else if(twogrids.getValue() == true)
      noise=noise1+noise2;  //if catalogues normalised and then summed, shot noise is sum of individual shot noises
    else noise=noise1;

    grid1->savePK(outfile.getValue(), 1, noise, header, myrank, root, com);

    delete [] sumscat1; //not used anymore
    delete [] sumscat2; //not used anymore
    delete [] sumsran1; //not used anymore
    delete [] sumsran2; //not used anymore
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

