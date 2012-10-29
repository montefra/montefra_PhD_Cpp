/*==========================================================================*/
/* Version: 0.10           date:19/01/12                                    */
/*                                                                          */
/* Author: Francesco Montesano, MPE, Garching                               */
/*                                                                          */
/* Purpose: calculate P22(k) as in equation 165                             */
/* in Bernaredau et al. 2002, Phys. Rep. 367, 1                             */
/*                                                                          */
/* WARNING:distances expressed in unit of h                                 */
/*==========================================================================*/

#define p22_CPP
#include "p22.h"


int main(int argc, char* argv[])
{
  try{    // Catch tclap ecceptions
    // Define the command line object, and insert a message
    // that describes the program. The second argument is the 
    // delimiter (usually space) and the last one is the version number. 

    //description of the code
    std::string description = "Compute 'P22(k)' as in equation 165 in Bernaredau et al. 2002, Phys. Rep. 367, 1. "; 
    description += "The integrations are done usind the 'gauleg'\n";

    TCLAP::CmdLine cmd(description, ' ', "0.10");
    
    //create input file name argument for the linear power spectrum and add to the command line parser
    TCLAP::UnlabeledValueArg<std::string> ifile("Plin(k)", "Input linear power spectrum file name (file format k, P(k))",
	                                  true, "", "Plin(k) string", cmd); 

    //create output file name argument for the P22 power spectrum and add to the command line parser
    TCLAP::UnlabeledValueArg<std::string> ofile("P22(k)", "Output P22 power spectrum file name (file format k, P(k))",
	                                  true, "", "P22(k) string", cmd); 

    //max and min value of k for the output P22(k)
    double ddef=0.5;   //set the default for kmax
    TCLAP::ValueArg<double> okmax("", "kmax", "Maximum wavenumber of P22(k). [Default: "+to_string(ddef)+" h/Mpc]",
				 false, ddef, "double", cmd);
    ddef=1e-3;   //set the default for kmin
    TCLAP::ValueArg<double> okmin("", "kmin", "Minimum wavenumber of P22(k). [Default: "+to_string(ddef)+" h/Mpc]",
				 false, ddef, "double", cmd);
    //number of bins in P22(k)
    int idef=300;   //set the default for the number of bins
    TCLAP::ValueArg<int> obins("n", "n-bins", "Number of bins in P22(k). Positive: linear binning, "+
	                       std::string("negative: logarithmic binning. [Default: "+to_string(idef)+"]"),
			       false, idef, "int", cmd);

    //max and min values of the wave number and number of bins used for the gauled integration
    TCLAP::ValueArg<double> rkmax("", "r-kmax", "Maximum wavenumber in the radial integration of P22(k). "+
	                          std::string("If not given the maximum k in Plin used"),
				  false, -1, "double", cmd);
    TCLAP::ValueArg<double> rkmin("", "r-kmin", "Minimum wavenumber in the radial integration of P22(k). "+
	                          std::string("If not given the minimum k in Plin used"),
				  false, -1, "double", cmd);
    idef=-1000;   //set the default for the number of bins of the radial integration
    TCLAP::ValueArg<int> rbins("", "r-bins", "Number of bins in the radial integraion of P22(k). Positive: linear binning, "+
	                       std::string("negative: logarithmic binning. [Default: "+to_string(idef)+"]"),
			       false, idef, "int", cmd);
    idef=300;   //set the default for the number of bins of the angular integration
    TCLAP::ValueArg<size_t> abins("", "a-bins", "Number of bins in the angular integraion of P22(k). "+
	                       std::string("[Default: "+to_string(idef)+"]"),
			       false, idef, "size_t", cmd);

    ddef=1e10;
    TCLAP::ValueArg<double> kstar("", "kstar", "Dumping of the linear power spectrum. If not given, no dumping applied",
	                    false, ddef, "double", cmd);


    // Parse the argv array.
    cmd.parse(argc, argv);

    /*==========================================================================*/
    /* check the files                                                          */
    /*==========================================================================*/
    //input file name
    if( fileexists(ifile.getValue()) == false) throw inputexception(ifile.getValue());
    //output file names
    if( fileexists(ofile.getValue()) == true ) throw outputexception(ofile.getValue());

    /*==========================================================================*/
    /* read the input file and save into a gsl interpolation object             */
    /* the variables are already declared in 'mcmc_incl_dec.h'                  */
    /*==========================================================================*/
    spllin = th2spl(ifile.getValue());
    acclin = gsl_interp_accel_alloc();    //allocate the accelerator for the linear power spectrum

    /*==========================================================================*/
    /* initialise the gauleg integration                                        */
    /*==========================================================================*/
    double kmax, kmin;  //save the values of kmax and kmin 
    if(rkmin.isSet() == true) kmin = rkmin.getValue();  //get the minimum k for the radial integration
    else kmin = spllin->x[0];  //if not set, use the min(k) of the input power spectrum
    if(rkmax.isSet() == true) kmax = rkmax.getValue();  //get the maximum k for the radial integration
    else kmax = spllin->x[spllin->size-1];  //if not set, use the max(k) of the input power spectrum

    p22 pMC( spllin, acclin );  //initilise the class computing the p22

    pMC.init_gauleg(kmax, kmin, rbins.getValue(), abins.getValue());

    /*==========================================================================*/
    /* compute the P22(k) and save in the outfile on the fly                    */
    /*==========================================================================*/
    /* open output file */
    ofstream out(ofile.getValue().c_str());
    out.setf(ios_base::scientific);
    out.precision(6);
    out.width(9);

    int i; //loop integer
    double k, p22;  //step by step
    int bins = obins.getValue(); //save the number of output bins
    kmax = okmax.getValue(); //save the values of kmax and kmin 
    kmin = okmin.getValue();
    if(bins<0){  //if logarithmic binning required
      kmax = log(kmax);
      kmin = log(kmin);
    }
    double binsize = (kmax - kmin)/abs(bins-1);   //size of k-bin

    if(kstar.isSet() == false){
      for(i=0; i<abs(bins); ++i){  //loop through the requred output bins
	k = binsize*i+kmin;   
	if(bins < 0) k = exp(k);   //if logarithmic binning

	p22 = pMC.computep22(k);    //compute p22 in k and return it
	out << k << "\t" << p22 << endl;  //print out  
      }
    }else{
      for(i=0; i<abs(bins); ++i){  //loop through the requred output bins
	k = binsize*i+kmin;   
	if(bins < 0) k = exp(k);   //if logarithmic binning

	p22 = pMC.computep22(k, kstar.getValue());    //compute p22 in k if the damping is required and return it
	out << k << "\t" << p22 << endl;  //print out  
      }
    }

    out.close();  //close the outputfile

    /*==========================================================================*/
    /* deallocations                                                            */
    /*==========================================================================*/
    gsl_spline_free(spllin);
    gsl_interp_accel_free(acclin);

  } 
  catch (TCLAP::ArgException &e){  // catch any exceptions
    cerr << "Error: " << e.error() << " for arg " << e.argId() << std::endl; 
    exit( 3);
  }
  catch (inputexception &e){
    cerr << "Error reading input file " << e.what() << "." << endl;
    exit( 4);
  }
  catch (outputexception &e){
    cerr << "The output file " << e.what() << " already exists. ";
    cerr << "Delete it, move it, rename it or change the output root." << endl;
    exit( 5);
  }
  catch (bad_alloc&)
  {
    cout << "Error allocating memory with 'new'." << endl;
    exit(6);
  }

  exit(0);
}

