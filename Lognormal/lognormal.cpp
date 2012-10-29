/*==========================================================================*/
/* version: 1.00, 05/09/2011                                                */       
/* author: Francesco Montesano (montefra@mpe.mpg.de)                        */
/*                                                                          */
/* Purpose: create lognormal mocks in cubic grids                           */
/*==========================================================================*/

#define MAIN_CPP
#include "lognormal.h"

int main(int argc, char* argv[])
{
  size_t i, j;    //loop integers
  ptrdiff_t ii, jj, ll;   // loop integers
  int errors;   //get output signal from functions in order to check the error

  /*==========================================================================*/
  /* Initialise MPI                                                           */
  /*==========================================================================*/
  int myrank, root=0, n_proc;    //rank of the processor, root processor (I use 0) and number of processes in the group of the communicator
  MPI_Comm comm = MPI_COMM_WORLD;   //MPI communicator used throught the code

  MPI_Init(&argc,&argv);
  MPI_Comm_rank(comm, &myrank);
  MPI_Comm_size(comm, &n_proc);

  //int debugging = 1;
  //while(debugging);
  /*==========================================================================*/
  /* Command line parsing interface                                           */
  /*==========================================================================*/
  string description;  //contains the description of the code. Then used to create the output file names
  description = "The code computes lognormal mock catalogues with positions and velocities according to ";
  description += "the given power spectrum and save them in ASCII files. ";
  description += "An initial power spectrum is needed and is transformed to a gaussian one, ";
  description += "which is then used to create the gaussian random fields\n ";
  description += "P(k) -> xi(r) -> xi_G(r)=ln(1+xi(r)) -> P_G(k)\n ";
  description += "Any of the above functions can be given as input. Which one is selected by an option. ";
  description += "P_G(k) can be saved for later use if required. ";
  description += "If the input power spectrum or correlation function are given in units of h, ";
  description += "positions and velocities are in units of Mpc/h";


  try{    // Catch tclap ecceptions
    // Define the command line object, and insert a message
    // that describes the program. The second argument is the 
    // delimiter (usually space) and the last one is the version number. 
    TCLAP::CmdLine cmd(description, ' ', "0.01");

    //product of H and a. If given, the velocity is printed in units of km/(h*s)
    TCLAP::ValueArg<double> H_a("", "H_times_a", "H(a)*a. If it is given the velocity is printed out in units of km/(s*h)",
	                            false, 1., "double", cmd);

    //path to a temporary directory where to save the output files from each processor
    TCLAP::ValueArg<string> tempdir("", "temp-dir", 
	                            string("Write the files for each processor into a temporary directory")+
				    string("If this option is set 'outroot' is split in path+file and the temporary ")+
				    string("files are saved as 'temp-dir'/'outrootfile_#processor'. Are then moved, ")+
				    string("merged und deleted. 'temp-dir' must not have the trailing '/'."),
				    false, "", "string", cmd);

    //number density or number of objects used to populate the box
    TCLAP::ValueArg<size_t> totalnumber("t", "total-number", 
	                                "Number of object to be used in each box. If this option is set '-d' ignored",
	                                false, 0, "unsigned int", cmd);
    TCLAP::ValueArg<double> ndensity("d", "number-density", "Number density of objects in [h/Mpc]^3. [Default: 1.e-4]",
	                            false, 1.e-4, "double", cmd);

    //number of output boxes
    TCLAP::ValueArg<size_t> nboxes("n", "number-boxes", "Number of lognormal catalogues to create. [Default: 100]",
	                         false, 100, "unsigned int", cmd);
    //offset used in the file name creation
    TCLAP::ValueArg<size_t> offset("o", "offset", 
	                           string("Offset to add to the number of the catalogue in the output file name.")+
				   string("Usefull for different calls to the code when creating a large number of catalogues"),
				   false, 0, "unsigned int", cmd);

    //number of cells assuming a square box
    TCLAP::ValueArg<ptrdiff_t> ncells("g", "cell-number", 
	                                "Number of cells per dimension of the box. The box is assumed cubic. [Default: 512]",
	                                false, 512, "unsigned int", cmd);
    //size of the output box. Assumed squared
    TCLAP::ValueArg<double> boxsize("l", "box-size", 
	                            "Size of the output box in Mpc/h. The box is assumed cubic. [Default: 2000]", 
	                            false, 2000, "double", cmd);

    //save P_G(k) for later use
    TCLAP::ValueArg<string> outP_Gk("s", "save-pgk", "Save P_G(k) for later use. If file already exists, file not written", 
				    false, "", "string", cmd);

    //maximum and minimum scales and binning scheme of the P(k)->xi(r) or inverse transform
    TCLAP::SwitchArg ftbin("", "ftbin", 
	                   "If the option is called, linear binning in the P(k)<->xi(r) FT used instead of logarithmic",
			   cmd, false);
    TCLAP::ValueArg<double> kmin("", "kmin", 
				 "Minimum value of the wavenumber in h/Mpc when xi(r)->P(k). [Default: " + to_string(dkmin) + "]",
				 false, dkmin, "double", cmd);
    TCLAP::ValueArg<double> kmax("", "kmax", 
				 "Maximum value of the wavenumber in h/Mpc when xi(r)->P(k). [Default: " + to_string(dkmax) + "]",
				 false, dkmax, "double", cmd);
    TCLAP::ValueArg<double> rmin("", "rmin", 
				 "Minimum value of the scale in Mpc/h when P(k)->xi(r). [Default: " + to_string(drmin) + "]",
				 false, drmin, "double", cmd);
    TCLAP::ValueArg<double> rmax("", "rmax", 
				 "Maximum value of the scale in Mpc/h when P(k)->xi(r). [Default: " + to_string(drmax) + "]",
				 false, drmax, "double", cmd);

    //binning of the input file, used to compute 'dk'
    vector<string> infunc;  //vector containing the function that can be used as input file

    infunc.push_back("lin");   //set the possible options: linear binning
    infunc.push_back("log");   //set the possible options: logarithmic binning
    infunc.push_back("irr");   //set the possible options: irregular binning
    TCLAP::ValuesConstraint<string> binTypes(infunc);
    infunc.clear();   //clear the vector for later use
    TCLAP::ValueArg<string> infileBinning("b", "binning", 
					  "Binning of the input file. Used only if P(k) or P_G(k) are read from file", 
					  false, "log", &binTypes, cmd);

    //number of columns in the input file
    TCLAP::ValueArg<int> inputcols("c", "columns-input", 
	                           "Number of columns of the input file. Must be >=2. [Default: 2]",
	                           false, 2, "int", cmd);

    //chose between P(k), xi(r), xi_G(r) and P_G(k). The first is the default
    infunc.push_back("Pk");
    infunc.push_back("xir");
    infunc.push_back("xi_Gr");
    infunc.push_back("P_Gk");
    TCLAP::ValuesConstraint<string> infileTypes(infunc);
    TCLAP::ValueArg<string> infileType("i", "inType", 
	                               "Input file can be any of these: P(k), xi(r), xi_G(r) and P_G(k). [Default: P(ka)]", 
				       false, "Pk", &infileTypes, cmd);
    
    //create input file name argument and add to the command line parser
    TCLAP::UnlabeledValueArg<string> infile("infile", "Input file name", true, "", "string", cmd); 

    //create the argument that contain f=dln(D)/dln(a), with D the linear growth factor and a the scale factor and add to the command line parser
    TCLAP::UnlabeledValueArg<double> ddot("f", "f=dln(D)/dln(a), with 'D' the linear growth factor and 'a' the scale factor", 
	                                  true, 0., "double", cmd); 

    //create output file name root argument and add to the command line parser
    TCLAP::UnlabeledValueArg<string> outfileroot("outroot", string("Output file root. File name: 'outroot.'+c+'.dat'. ")+
	                                         string("'c' is a counter with an appropriate number of padding 0."), 
	                                         true, "", "string", cmd); 

    // Parse the argv array.
    cmd.parse(argc, argv);

    /*==========================================================================*/
    /* check the output file names and throw an exeption if any of them exists  */
    /* for MPI creates a file per process and then merge them                   */
    /*==========================================================================*/
    int n_digits;  // number of digits on the input number of catalogues to be produced
    size_t numboxes = nboxes.getValue();  //number of boxes to be generated

    if( numboxes < 1 ){
      cerr << "At least one mock must be created" << endl;
      MPI_Abort(comm, 66);
    }
    n_digits = getNumberOfDigits(numboxes + offset.getValue());  //calculates the number of digits
    if(myrank == root){
      for(i=0; i<numboxes; ++i){   //loop through the required boxes
	description = outfileroot.getValue() + "." + padwithzero(i+1+offset.getValue(), n_digits) + string(".dat");  //get the file root
	if(fileexists(description) == true) throw outputexception(description);
      }
    }

    /*==========================================================================*/
    /* save the output root, file name root and the temporary root              */
    /*==========================================================================*/
    string outrootpath, outroottemppath, outrootfname;   //strings contining the output path, the temporary output path 
							 //and the output root file name
    outrootpath = extractPath(outfileroot.getValue());       //get the output path
    outrootfname = extractFilename(outfileroot.getValue());  //get the output file name root
    if(tempdir.isSet() == 0) outroottemppath = outrootpath;  //if a temporary path is not given used the output one
    else{
      outroottemppath = tempdir.getValue();               //otherwise save it
      errors = mkdir(outroottemppath.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);   //try to create the output temporary dir
      if(errors != 0 && errors != -1){   //if the directory is not created and doesn't already exists
	cerr << "I wasn't able to create the temporary directory. Used the output one to store the temporary files." << endl;
	outroottemppath = outrootpath; 
      }
    }

    /*==========================================================================*/
    /* get the value of f=dln(D)/dln(a) from the command line and save in a     */
    /* variable. If H and a are also provided they are multiplied to f          */
    /* and the velocity is printed out in km/(h*s) instead of Mpc/h             */
    /*==========================================================================*/

    double f;   //f=dln(D)/dln(a)  

    f = ddot.getValue();  //get f
    f *= H_a.getValue();  //multiply by H*a. if not given the default is 1, so there is not difference

    /*==========================================================================*/
    /* read the input file and do any of the steps required to get P_G(k)       */
    /* according to the requirements from the option '-i'                       */
    /* P(k) -> xi(r) -> xi_G(r)=ln(1+xi(r)) -> P_G(k)                           */
    /* the array containing dk is allocated and filled                          */
    /*==========================================================================*/
    double * r, * xi, * k, * pk;   //arrays containing correlation and power spectra
    size_t pkxisize;  //size of the input arrays

    if( infileType.getValue().compare(infunc[0]) == 0){    //if the input file is P(k)
      pkxisize = readandconvert(infile.getValue(), inputcols.getValue(), k, pk);   //read the input file and returns arrays
      allocarray(xi, pkxisize, __FILE__, __LINE__);    //allocate xi
      allocarray(r, pkxisize, __FILE__, __LINE__);    //allocate r
      fillarray(r, pkxisize, rmin.getValue(), rmax.getValue(), ftbin.getValue());
      pk2corfunc(k, pk, pkxisize, r, xi, pkxisize);  //convert the power spectrum into the correlation function and get the values of r
      logxi(xi, pkxisize);   //xi_G = ln(1+xi)
      corfunc2pk(r, xi, pkxisize, k, pk, pkxisize);   //convert the xi_G(r) into P_G(k)
      delete[] r;   //deallocate 'r' not used anymore
      delete[] xi;   //deallocate 'xi' not used anymore
    }
    else if( infileType.getValue().compare(infunc[1]) == 0){    //if the input file is xi(r)
      pkxisize = readandconvert(infile.getValue(), inputcols.getValue(), r, xi);   //read the input file and returns arrays
      logxi(xi, pkxisize);   //xi_G = ln(1+xi)
      allocarray(pk, pkxisize, __FILE__, __LINE__);    //allocate xi
      allocarray(k, pkxisize, __FILE__, __LINE__);    //allocate r
      fillarray(k, pkxisize, kmin.getValue(), kmax.getValue(), ftbin.getValue());
      corfunc2pk(r, xi, pkxisize, k, pk, pkxisize);   //convert the xi_G(r) into P_G(k)
      delete[] r;   //deallocate 'r' not used anymore
      delete[] xi;   //deallocate 'xi' not used anymore
    }
    else if( infileType.getValue().compare(infunc[2]) == 0){    //if the input file is xi_G(r)
      pkxisize = readandconvert(infile.getValue(), inputcols.getValue(), r, xi);   //read the input file and returns arrays
      allocarray(pk, pkxisize, __FILE__, __LINE__);    //allocate xi
      allocarray(k, pkxisize, __FILE__, __LINE__);    //allocate r
      fillarray(k, pkxisize, kmin.getValue(), kmax.getValue(), ftbin.getValue());
      corfunc2pk(r, xi, pkxisize, k, pk, pkxisize);   //convert the xi_G(r) into P_G(k)
      delete[] r;   //deallocate 'r' not used anymore
      delete[] xi;   //deallocate 'xi' not used anymore
    }
    else    //if the input file is P_G(k)
      pkxisize = readandconvert(infile.getValue(), inputcols.getValue(), k, pk);   //read the input file and returns arrays
    /*==========================================================================*/
    /* if required save k, P_G(k) in a file for later use                       */
    /*==========================================================================*/
    if(myrank == root){   //only the root write the output file
      if( outP_Gk.isSet() ){
	if( writetwoarrays(outP_Gk.getValue(), k, pk, pkxisize) == 1)
	  cerr << "The file " << outP_Gk.getValue() << " already exists. P_G(k) won't be saved."  << endl;
      }
    }
    MPI_Barrier(comm);   //wait for the root to write the output file name
    if(myrank == root) cout << "Input file read and processed." << endl;


    /*==========================================================================*/
    /* Compute sigma: is the variance of the gaussian power spectrum            */
    /* field used to fill the grid                                              */
    /* also add k=0 and p(k=0)=0 at the beginning of the k and pk arrays        */
    /*==========================================================================*/
    
    double b, h;    //box and cell size (temporarly volume of the box) in Mpc/h

    b = boxsize.getValue();   // size of the box. Modify for non cubic grid
    h = b*b*b;   //save the volume of the box temporarly

    allocarray(r, pkxisize, __FILE__, __LINE__);    //allocate 'r' for temporary use
    for(i=0; i<pkxisize; ++i) r[i] = pk[i];  //copy pk
    delete[] pk;   //deallocate pk
    allocarray(pk, pkxisize+1, __FILE__, __LINE__);    //reallocate pk with one more element for themporary storing
    pk[0] = 0.;  //set to zero the first element
    for(i=0; i<pkxisize; ++i) pk[i+1] =  sqrt( r[i]/h );   //fill pk with sigma(k)=sqrt(P(k)/V)
    for(i=0; i<pkxisize; ++i) r[i] = k[i];  //copy k in r
    delete[] k;   //deallocate k
    allocarray(k, pkxisize+1, __FILE__, __LINE__);    //reallocate k with one more element for themporary storing
    k[0] = 0.;  //set to zero the first element
    for(i=0; i<pkxisize; ++i) k[i+1] = r[i];   //copy back k

    delete[] r;    //delete the array 'r' not used anymore

    /*==========================================================================*/
    /* initialise the spline of sigma^2(k)                                      */
    /*==========================================================================*/
    gsl_spline *spl;   // declaration of the interpolation scheme
    gsl_interp_accel *acc;   // declaration of the accelarator for the interpolation;

    spl = gsl_spline_alloc(gsl_interp_akima, pkxisize);   //alloc a gsl_spline
    if(spl == 0){
      cerr << __FILE__ << ": " << __LINE__-2 << ". ERROR: the allocation of the spline object failed" << endl;
      exit(10);
    }
    errors= gsl_spline_init(spl, k, pk, pkxisize);
    if( errors != 0 ){
      cerr << __FILE__ << ": " << __LINE__-1 << ". ";
      if( errors == 4 ) cerr << "ERROR: x values must be monotonically increasing" << endl;
      else cerr << "ERROR: initialisation of the gsl spline failed" << endl;
      exit(20);
    }
    acc = gsl_interp_accel_alloc();   //allocate the accelerator for the interpolation
    if(spl == 0){
      cerr << __FILE__ << ": " << __LINE__-2 << ". ERROR: the allocation of the interpolation accelerator failed" << endl;
      exit(11);
    }

    delete[] k;   //delete the array k and pk. They content is saved in the spline
    delete[] pk;

    if(myrank == root) cout << "Variance computed and saved into an interpolation object." << endl;

    /*==========================================================================*/
    /* compute the values of k associated with the grid in Fourier space and    */
    /*==========================================================================*/

    double kN;   //Nyquist frequency
    ptrdiff_t* xyzcells;  //save the number of cells from the command line options
    ptrdiff_t totcells; //total number of cells
    ptrdiff_t* xyzhalf;  //half of the dimensions

    allocarray(xyzcells, 3, __FILE__, __LINE__);    //allocate xyzcells with dimension 3
    xyzcells[0] = xyzcells[1] = xyzcells[2] = ncells.getValue();  //save the number of cells from the command line options. Modify for non cubic grid
    totcells = xyzcells[0] * xyzcells[1] * xyzcells[2];  //total number of cells
    allocarray(xyzhalf, 3, __FILE__, __LINE__);    //allocate xyzcells with dimension 3
    for(i=0; i<3; ++i) xyzhalf[i] = xyzcells[i]/2;     //save half of the number of cells, rounded to down
    h = b/double(xyzcells[0]);       //size of a cell in Mpc/h: L_{box}/N_{cells}. Modify for non cubic grid
    kN = M_PI / h;       //set the Nyqist frequency

    allocarray(k, xyzcells[0], __FILE__, __LINE__);    //allocate k to get the values of k on the grid
    
    k[0] = 0.;   //set to 0 the first element of the k array
    for(ii=1; ii<xyzhalf[0]; ++ii){
      k[ii] = 2.*M_PI/boxsize.getValue()*ii;    //fill almost all the array with k with k[i]=k[n-i] for 1<=i<xyzcells[0]/2
      k[xyzcells[0]-ii] = -k[ii];    //fill almost all the array with k with k[i]=k[n-i] for 1<=i<xyzcells[0]/2
    }
    if(xyzcells[0]%2 == 0) k[xyzhalf[0]] = kN;   //if the number of cells is even, the xyzhalf[0] element of k is kN 
    else k[xyzhalf[0]] = k[xyzhalf[0]+1] = 2.*M_PI/boxsize.getValue()*xyzhalf[0];     //otherwise k[xyzhalf[0]] = -k[xyzhalf[0]+1] and they are both < kN

    /*==========================================================================*/
    /* set the random generator                                                 */
    /*==========================================================================*/

    gsl_rng* rnd;   //random generator object

    //rnd = gsl_rng_alloc(gsl_rng_ranlxs1);   //allocate the random generator object
    rnd = gsl_rng_alloc(gsl_rng_taus2);   //allocate the random generator object
    MPI_Barrier(comm);   //wait for all the processors before setting the seed to the c random generator
    srand(time(0)+myrank);   //seed the c++ random generator rand with a different seed for each processor

    /*==========================================================================*/
    /* Initialize the grid for the FFTW                                         */
    /*==========================================================================*/
    ptrdiff_t local_nx, local_x_start;  //number of cells, number of cells in the x direction and 
    ptrdiff_t vlocal_nx, vlocal_x_start;  //number of cells, number of cells in the x direction and 
                                                   //first cell in the x direction per processor
    fftw_r2c_c2r lngrid(xyzcells[0], comm, &local_nx, &local_x_start, true);  //constructor for mpi with inplace transform
    //lngrid.read_wisdom(root, myrank, comm);  //MPI: wisdom read by in the root processor and then broadcasted to 'comm'
    lngrid.create_plan(FFTW_BOTH);  //create the forward and backwark plans for FFTW.
    //lngrid.write_wisdom(root, myrank, comm);            //MPI: write wisdom   

    fftw_r2c_c2r velgrid(xyzcells[0], comm, &vlocal_nx, &vlocal_x_start, true);  //constructor for mpi with inplace transform for the fft of the velocity field
    velgrid.create_plan(FFTW_BACKWARD);  //create both the backward plan for FFTW for the velocity field.

    if( local_nx != local_nx && local_x_start != vlocal_x_start ){
      cerr << "The lognormal and the velocity field grid are not divided in the same way among processors" << endl;
      MPI_Abort(comm, 6);
    }

    if(myrank == root) cout << "Grids for the fft initialised." << endl;
 
    /*==========================================================================*/
    /* set the number of particles per processor                                */
    /*==========================================================================*/
    size_t totparticles;    //total number of particles

    if(totalnumber.isSet() != 0) totparticles = totalnumber.getValue();    //if the total number of particles given, use it
    else    //otherwise compute it from the density
      totparticles = ndensity.getValue() * b*b*b;

    /*==========================================================================*/
    /* Loop through the number of required realizations                         */
    /*==========================================================================*/

    double modd, phi;   //modulus and phase of the complex density field
    double modk, sigmak;   //|k|=sqrt(kx^2 + ky^2 + kz^2 ) and the corresponding value of the variance
    double maxd;   //contains the sum (or the max) of delta(x)
    double maxdroot;   //maxd on the root used to compute the absolute sum (or maximum) of dr
    ofstream out;  //outputstream with the file name
    double x, y, z;   //coordinates of the output particle
    size_t localnpart;  //local number of particles for the random positions way
    vector<ptrdiff_t> saved;   //list the cells to which the particles saved belongs. Used to save their velocities

    MPI_Barrier(comm);   //wait before beginning all the loops
    if(myrank == root) cout << "Now beginning with the production of the LN catalogues" << endl;

    for(i=0; i<numboxes; ++i){   //loop through the required boxes
      gsl_rng_set(rnd, rand());   //in each loop seed the generator with a different random number from the c++ random generator
      maxd = 0.;   //reinitialised to 0 at each step
      //lngrid.initialise();   //initialize the grid with all zeros

      /*==========================================================================*/
      /* fill the complex grid with gaussian random numbers                       */
      /*==========================================================================*/

      for(ii=0; ii<local_nx; ++ii) for(jj=0; jj<xyzcells[1]; ++jj)   //fill the complex grid for with random gaussian points
	for(ll=0; ll<xyzhalf[2]+1; ++ll){   //optimised for MPI
	  //loop through the cells of the fftw array
	  modk = sqrt( k[ii+local_x_start]*k[ii+local_x_start] + k[jj]*k[jj] + k[ll]*k[ll] );   //modulus of the wave number
	  sigmak = gsl_spline_eval(spl, modk, acc);   //get the value of sigma at modk

	  //modd = gsl_ran_gaussian(rnd, sigmak);   //get modulus and the phase of the complex density field
	  modd = gsl_ran_gaussian_ziggurat(rnd, sigmak);   //get modulus and the phase of the complex density field
	  phi = gsl_rng_uniform(rnd)*2.*M_PI;
	  //assign the complex amplitued delta = modd( cos(phi) + i*sin(phi) )
	  lngrid.set_value_complex(ii, jj, ll, 0, modd*cos(phi));
	  lngrid.set_value_complex(ii, jj, ll, 1, modd*sin(phi));
        }

      lngrid.make_hermitian();   //make the complex grid hermitian
      /*==========================================================================*/
      /* do the fourier transform                                                 */
      /*==========================================================================*/
      lngrid.execute_fft(FFTW_BACKWARD);

      /*==========================================================================*/
      /* compute the log normal density field delta_LN(x) = exp(delta(x))         */
      /* and extract beta= Npart/sum(dr) (or the maximum of dr)                   */
      /*==========================================================================*/
      for(ii=0; ii<local_nx; ++ii) for(jj=0; jj<xyzcells[1]; ++jj) for(ll=0; ll<xyzcells[2]; ++ll){   
	modd = exp( lngrid.get_value_real(ii, jj, ll) );   //get the value of the real part of the grid and take the exponent of it
	lngrid.set_value_real(ii, jj, ll, modd);   //reassigne the value 
	maxd += modd;
	//maxd = (modd>maxd) ? modd : maxd;   //if new 'dr' larger than max, 
      }
      //find the absolute maximum among the processors
      //MPI_Reduce(&maxd, &maxdroot, 1, MPI_DOUBLE, MPI_MAX, root, comm);   //get all the value and check for the absolute maximum
      MPI_Reduce(&maxd, &maxdroot, 1, MPI_DOUBLE, MPI_SUM, root, comm);   //get all the value and check for the absolute maximum
      if(myrank == root) maxd = maxdroot;   // on the root copy the absolute maximum into the given variable
      MPI_Bcast(&maxd, 1, MPI_DOUBLE, root, comm);   //broadcast the absolute maximum to all the processors
      maxd = totparticles/maxd;   //save beta in 'maxd'

      /*==========================================================================*/
      /* Create the catalogue with x, y and z position and print to file          */
      /*==========================================================================*/

      description = outroottemppath + "/" + outrootfname + "." + padwithzero(i+1+offset.getValue(), n_digits) + 
	            string("_position_") + to_string(myrank) + string(".dat");  //set the output file name for each processor
      out.open(description.c_str());  //open the output file
      out.setf(ios_base::scientific);   //set precision and style
      out.precision(6);
      out.width(9);

//      for(ii=0; ii<local_nx; ++ii) for(jj=0; jj<xyzcells[1]; ++jj) for(ll=0; ll<xyzcells[2]; ++ll){   //loop through the cells
//	if( maxd*lngrid.get_value_real(ii, jj, ll) >= gsl_rng_uniform(rnd) ){
//	  x = gsl_rng_uniform(rnd);   //add a random displacement in the cell
//	  y = gsl_rng_uniform(rnd);
//	  z = gsl_rng_uniform(rnd);
//	  out << double(ii+local_x_start+x)*h << "\t" << double(jj+y)*h << "\t" << double(ll+z)*h << endl;  //write the valid coordinates
//	  saved.push_back( lngrid.ijl_to_ind_real(ii, jj, ll) );   //save the position of the particle in the array. 
//	                                                           //This allow to use one index instead of the three
//	}
//      }

//      if(myrank == root) localnpart = gsl_ran_poisson(rnd, totparticles);  //the root get a random poisson number
//      MPI_Bcast(&localnpart, 1, MPI_UNSIGNED_LONG, root, comm);   //broadcast the absolute maximum to all the processors
//      localnpart = (double)local_nx* localnpart / (double)xyzcells[0];   //poisson distibution of particles with mean 'totparticles' then distributed among the processors
//      //localnpart = (double)local_nx*totparticles / (double)xyzcells[0];   //distribute the number of particles among the processors
//      j=0;   //initialise the counter
//      while(j<localnpart){
//	//get random coordinates in cells units
//	x = gsl_rng_uniform(rnd)*local_nx;
//	y = gsl_rng_uniform(rnd)*xyzcells[1];
//	z = gsl_rng_uniform(rnd)*xyzcells[2];
//	if( maxd*gsl_rng_uniform(rnd) <= lngrid.get_value_real(ptrdiff_t(x), ptrdiff_t(y), ptrdiff_t(z)) ){
//	  out << (x+double(local_x_start))*h << "\t" << y*h << "\t" << z*h << endl;  //write the valid coordinates
//	  saved.push_back( lngrid.ijl_to_ind_real(ptrdiff_t(x), ptrdiff_t(y), ptrdiff_t(z)) );   //save the position of the particle in the array. 
//	  j++;
//	}
//      }

      for(ii=0; ii<local_nx; ++ii) for(jj=0; jj<xyzcells[1]; ++jj) for(ll=0; ll<xyzcells[2]; ++ll){   //loop through the cells
	localnpart = gsl_ran_poisson(rnd, maxd*lngrid.get_value_real(ii, jj, ll));   //get a poisson random integer with mean 'maxd*delta(x,y,z)
	for(j=0; j<localnpart; ++j){
	  x = 0.; // gsl_rng_uniform(rnd);   //add a random displacement in the cell
	  y = 0.; // gsl_rng_uniform(rnd);
	  z = 0.; // gsl_rng_uniform(rnd);
	  out << double(ii+local_x_start+x)*h << "\t" << double(jj+y)*h << "\t" << double(ll+z)*h << endl;  //write the valid coordinates
	  saved.push_back( lngrid.ijl_to_ind_real(ii, jj, ll) );   //save the position of the particle in the array. 
	}
      }

      out.close();  //close the output file

      /*==========================================================================*/
      /* go back to k space and work on the velocities                            */
      /*==========================================================================*/

      lngrid.execute_fft(FFTW_FORWARD);   //transform the delta(x) -> delta(k)

      localnpart = saved.size(); //number of particles saved by each processor

      int l;  //loop through the three directions of in order to compute the velocity for x,y,z
      for(l=0; l<3; ++l){

	if(l == 0) _return_index = &_return_index_ii;        //assign the pointer to the right function 
	else if(l == 1) _return_index = &_return_index_jj;   //according to the direction of the velocity
	else _return_index = &_return_index_ll;              //wanted

	for(ii=0; ii<local_nx; ++ii) for(jj=0; jj<xyzcells[1]; ++jj) for(ll=0; ll<xyzhalf[2]+1; ++ll){   //optimised for MPI

	  //square of the modulus of the wave number and multiply by the total number of cells 
	  //in order to obtain the right amplitude of the fluctuations
	  modk = (k[ii+local_x_start]*k[ii+local_x_start] + k[jj]*k[jj] + k[ll]*k[ll]) * totcells;
	  kN = k[return_index(ii+local_x_start, jj, ll)];    //get the correct value of k according to the direction
	  // u_n(k) propto i k_n/k^2 delta(k) with n=x,y,z
	  velgrid.set_value_complex(ii, jj, ll, 0, -kN*lngrid.get_value_complex(ii, jj, ll, 1)/modk );  //fill the grid with the
	  velgrid.set_value_complex(ii, jj, ll, 1, kN*lngrid.get_value_complex(ii, jj, ll, 0)/modk );    //u_n(k)
	}
	if(local_x_start == 0){
	  velgrid.set_value_complex(0, 0, 0, 0, 0.);   //set to zero the fluctuations in for |k|=0
	  velgrid.set_value_complex(0, 0, 0, 1, 0.);
	}

	//velgrid.make_hermitian();   //make the complex grid hermitian

	velgrid.execute_fft();  //execute the fourier transform of the velocity field

	/*==========================================================================*/
	/* Create a file with the velocities of the particles previously saved      */
	/*==========================================================================*/
	description = outroottemppath + "/" + outrootfname + "." + padwithzero(i+1+offset.getValue(), n_digits) + 
		      string("_vel_") + to_string(l) + string("_") + to_string(myrank) + string(".dat");  //set the output file name for each processor
	out.open(description.c_str());  //open the output file
	out.setf(ios_base::scientific);   //set precision and style
	out.precision(6);
	
	//print out the file with the velocities in the l-th direction corresponding to the particles printed
	for(j=0; j<localnpart; ++j) out << f * velgrid.get_value_real_ijl(saved[j]) << endl;

	out.close();  //close the output file
      }

      saved.clear(); //the vector of the saved particles is cleared before making the next catalogue

      //paste the position and velocity files together in a temporary file
      description = outroottemppath + "/" + outrootfname + "." + padwithzero(i+1+offset.getValue(), n_digits) + 
	            string("_*_") + to_string(myrank) + string(".dat");  //set the output file name for each processor
      string temp = outroottemppath + "/" + outrootfname + "." + padwithzero(i+1+offset.getValue(), n_digits) + 
		    string("_") + to_string(myrank) + string(".dat");  //set the output file name for each processor
      system( ("paste " + description + " > " + temp).c_str() );
      system( ("rm " + description).c_str() );  //delete the file with positions and velocities

      //if a temporary directory used move all the files to the required one
      if(outrootpath.compare(outroottemppath) != 0) system( ("mv " + temp + " " + outrootpath + "/").c_str() );

      MPI_Barrier(comm);   //wait for all the processor to paste their position and velocity files and merge them together
      if(myrank == root){   //the root merge all the files into a unique one with cat and delete the temporary ones
	description = "cat ";  //cat temp_files > final_file
	description += outrootpath + "/" + outrootfname + "." + padwithzero(i+1+offset.getValue(), n_digits) + "_*"; 
	description += " > ";   
	description += outrootpath + "/" + outrootfname + "." + padwithzero(i+1+offset.getValue(), n_digits) + ".dat";
	system(description.c_str());   //execute the merging

	description = "rm ";   //rm temp_files
	description += outrootpath + "/" + outrootfname + "." + padwithzero(i+1+offset.getValue(), n_digits) + "_*"; 
	system(description.c_str());   //execute the deletion

      }
      MPI_Barrier(comm);   //wait for the root to merge the files and delete the temporary ones

      if( (myrank == root) && ((i+1)%10 == 0) ) cout << "Catalogue # " << i+1 << " created" << endl;
    }

    delete[] k;   //deallocate the array containing k, not used anymore
    gsl_interp_accel_free(acc);    //free the gsl spline accelerator
    gsl_spline_free(spl);   //free the gsl spline object
    gsl_rng_free(rnd);   //free the random generator
  } 
  catch (TCLAP::ArgException &e){  // catch any exceptions
    cerr << "Error: " << e.error() << " for arg " << e.argId() << std::endl; 
    MPI_Abort(comm, 3);
  }
  catch (inputexception &e){
    cerr << "Error reading input file " << e.what() << "." << endl;
    MPI_Abort(comm, 4);
  }
  catch (outputexception &e){
    cerr << "The output file " << e.what() << " already exists. ";
    cerr << "Delete it, move it, rename it or change the output root." << endl;
    MPI_Abort(comm, 5);
  }

  /*==========================================================================*/
  /* Finalize MPI                                                             */
  /*==========================================================================*/
  MPI_Finalize();

  exit(0);
}
