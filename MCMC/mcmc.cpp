/*==========================================================================*/
/* Version: 4.00           date:17/01/12                                    */
/*                                                                          */
/* Author: Francesco Montesano, MPE, Garching                               */
/*                                                                          */
/* Purpose: implement the Markov Chain Monte Carlo calculating the          */
/* likelihood function with the RPT prescription                            */
/*==========================================================================*/

#define MCMC_CPP
#include "mcmc.h"

/*==========================================================================*/
/*                                 MAIN                                     */
/*==========================================================================*/

int main(int argc, char* argv[])
{
  int i,j;   //loop integer

  try{    // Catch ecceptions
    /*==========================================================================*/
    /* Command line parsing interface description                               */
    /*==========================================================================*/
    std::string description;  //contains the description of the code. Then used to create the output file names
    description = "Run MonteCarlo Markov Chains applying to the power spectrum ";
    description += "given the model in function 'double p_obj(double k, double *par). ";
#ifndef WINMAT    
    description += "The root of the output files, file names of the measured power spectrum and of the covariance matrix ";
#else
    description += "The root of the output files, file names of the measured power spectrum, of the ";
    description += "covariance matrix,of  the window matrix, of the values of k where to ";
    description += "evaluate the model before convolving with the window matrix, ";
    description += "of the line of the window matrix evaluated at k=0 and ";
    description += "of the window function evaluated in k=0 and in the k of the input power spectrum " ;
#endif
    description += "are mandatory.\n";

    // Define the command line object, and insert a message
    // that describes the program. The second argument is the 
    // delimiter (usually space) and the last one is the version number. 
    TCLAP::CmdLine cmd(description, ' ', "4.00");

    //create input file name argument for the power spectrum and add to the command line parser
    TCLAP::UnlabeledValueArg<std::string> Pkfn("P(k)", "Power spectrum file name (file format k, P(k), variance)",
	true, "", "P(k) string", cmd); 
    //create input file name argument for the covariance matrix and add to the command line parser
    TCLAP::UnlabeledValueArg<std::string> covfn("covariance", std::string("Covariance matrix file name")+ 
	std::string(" (file format i,j,ki,kj, covariance)"), true, "", "Covariance (string)", cmd); 
#ifdef WINMAT    
    //create input file name argument for the window matrix and add to the command line parser
    TCLAP::UnlabeledValueArg<std::string> wijfn("W_ij", "Window matrix file name (matrix file)", true, "", "W_ij (string)", cmd); 
    //create input file name argument for kj and add to the command line parser
    TCLAP::UnlabeledValueArg<std::string> kjfn("k_j", std::string("File name of the file containing the values where to ")+
	std::string("evaluate the model before convolving with W_ij (one column file)"), 
	true, "", "k_j (string)", cmd); 
    //create input file name argument for the power spectrum and add to the command line parser
    TCLAP::UnlabeledValueArg<std::string> w0jfn("W_0j", "W(ki=0,kj), used for the integral contraint (one column file)", 
	true, "", "W_0j (string)", cmd); 
    //create input file name argument for the covariance matrix and add to the command line parser
    TCLAP::UnlabeledValueArg<std::string> G2kifn("G^2(k_i)", "Window function estimated in the k of P(k)", 
	true, "", "G^2(k_i) (string)", cmd); 
#endif
    //inifile
    TCLAP::UnlabeledValueArg<std::string> inifile("inifile", std::string("File containing the parameters initial guess, range ") + 
	std::string("and width of the step in the MCMC chain"), true, "", "inifile (string)", cmd);

    //create output file name root argument and add to the command line parser
    TCLAP::UnlabeledValueArg<std::string> outroot("outroot", std::string("Output file root. Two files created: ")+
	std::string("'outroot.dat': MCMC chain (format: weight, likelihood, parameters); 'outroot.paramnames': ")+ 
	std::string("parameter name file (coded in the function 'void paramnames(string ofile, int n_params)')."), 
	true, "", "outroot (string)", cmd); 

    //linear power spectrum
    std::string dfname = "/data01/montefra/BOSS/Theory/lasdamas_matterpower_z0.55.dat";  //default file name of the linear ps
    TCLAP::ValueArg<std::string> linPk("l", "lin-pk", "Name of the linear power spectrum. [Default: " + dfname + "]", 
	false, dfname, "string", cmd);
    //1loop power spectrum
    dfname = "/data01/montefra/BOSS/Theory/lasdamas_1loop_z0.55.dat";  //default file name of the 1 loop ps
    TCLAP::ValueArg<std::string> mcPk("1", "1loop-pk", std::string("Name of the 1loop mode ") + 
	std::string("coupling power spectrum. [Default: ") + dfname + "]", false, dfname, "string", cmd);

    //enable the dumping of the 1 loop power spectrum
    TCLAP::SwitchArg dumpMC("d", "dump-1loop", "Turn on the dumping of the 1loop power spectrum", cmd);

    //include redshift space distortions 
    TCLAP::ValueArg<double> f("f", "log_der_D", std::string("f=dln(D(a))/dln(a). If f>0 model of redshift space distortions ")+
	std::string("used; otherwise not used."), false, -1., "double", cmd);

    // Unbias the inverse covariance
    //Hartlap et al. 2007: "the inverse of the maximum-likelihood
    //estimator of the covariance is biased".
    TCLAP::ValueArg<int> unbias("u", "unbias-cov", 
	"Number of realisations used to estimate the covariance matrix. If positive unbias the inverse covariace matrix as in Hartlap et al. 2007.",
        false, -1, "int", cmd);

#ifdef WINMAT    
    //max and min value of kj where to evaluate the model before convolving with W_ij
    TCLAP::ValueArg<double> ckjmax("", "kjmax", std::string( "Maximum value of the wavenumber in h/Mpc ")+ 
	std::string("where to evaluate the model. [Default: maximum of k_j]"), false, -1., "double", cmd);
    TCLAP::ValueArg<double> ckjmin("", "kjmin", std::string("Minimum value of the wavenumber in h/Mpc ")+ 
	std::string("where to evaluate the model. [Default: minimum of k_j]"), false, -1., "double", cmd);
#endif
    //max and min value of k where to compare model and observed P(k)
    TCLAP::ValueArg<double> ckmax("", "kmax", 
	"Maximum value of the wavenumber in h/Mpc used in the fit. [Default: maximum k in P(k)]", false, -1., "double", cmd);
    TCLAP::ValueArg<double> ckmin("", "kmin", 
	"Minimum value of the wavenumber in h/Mpc used in the fit. [Default: minimum k in P(k)]", false, -1., "double", cmd);

    int dmax_step=300000;  //default number of steps
    TCLAP::ValueArg<int> max_step_cmd("m", "max-step", "Number of steps of the MCMC chain. [Default: " + to_string(dmax_step) + "]", 
	false, dmax_step, "int", cmd);

    //verbose mode
    TCLAP::SwitchArg verbose("v", "verbose", "Verbose mode", cmd);

    // Parse the argv array.
    cmd.parse(argc, argv);

    /*==========================================================================*/
    /* check input and output files                                             */
    /*==========================================================================*/

    if( fileexists(Pkfn.getValue()) == false)   //measured power spectrum
      throw inputexception(Pkfn.getValue());    
    if( fileexists(covfn.getValue()) == false)  //covariance
      throw inputexception(covfn.getValue());
#ifdef WINMAT    
    if( fileexists(wijfn.getValue()) == false)  //window function files
      throw inputexception(wijfn.getValue());
    if( fileexists(kjfn.getValue()) == false)
      throw inputexception(kjfn.getValue());
    if( fileexists(w0jfn.getValue()) == false)
      throw inputexception(w0jfn.getValue());
    if( fileexists(G2kifn.getValue()) == false)
      throw inputexception(G2kifn.getValue());
#endif
    if( fileexists(linPk.getValue()) == false)  //linear and 1loop power spectra
      throw inputexception(linPk.getValue());
    if( fileexists(mcPk.getValue()) == false)  //linear and 1loop power spectra
      throw inputexception(mcPk.getValue());

    //check the inifile
    if( fileexists(inifile.getValue()) == false)   //measured power spectrum
      throw inputexception(inifile.getValue());    

    //output file names
    std::string ofiles[2];   // output file with the mcmc chain and the file with the parameter names
    ofiles[0] = outroot.getValue()+".dat";
    ofiles[1] = outroot.getValue()+".paramnames";
    for(i=0; i<2; i++) if( fileexists(ofiles[i]) == true )
	throw outputexception(ofiles[i]);

    /*==========================================================================*/
    /* allocate the parameters array and set the default values                 */
    /*==========================================================================*/
    int n_params;    // number of parameters. WARNING hard coded in 'int paramnames(string ofile, vector<string> &pars)'
    std::vector<std::string> parsname;   //contain the short name of the parameters. WARNING hard coded in 'int paramnames(string ofile, vector<string> &pars)'
    double *params, *params1, *upper, *lower, *sigma;   //array containing the parameter values, the upper and lower boundaries and the sigma
    n_params = paramnames(ofiles[1], parsname);
    params = new double[n_params];
    params1 = new double[n_params];
    upper = new double[n_params];
    lower = new double[n_params];
    sigma = new double[n_params];
    //if the ini file is given substitute the default values with the custom ones
    readini(inifile.getValue(), parsname, params, upper, lower, sigma);
    if( verbose.getValue() == true)
      std::cout << "Inifile read" << std::endl;

    /*==========================================================================*/
    /* check cosmology related parameters                                       */
    /*==========================================================================*/
    cosmo.f = f.getValue();   //get the value of f=dln(D(a))/dln(a)

    if( cosmo.f > 0 && params[4] <= 0 && sigma[4] == 0){
      std::cerr << "To model redshift space distortions the bias must be a free parameter" << std::endl;
      exit(7);
    }

    /*==========================================================================*/
    /* copy the maximum and mimimum values of k to be used                      */
    /*==========================================================================*/
    double kimax, kimin; //maximum and minimum k to be used in the measured power spectrum
#ifdef WINMAT    
    double kjmax, kjmin; //maximum and mimimum k to be used in the model power spectrum before convolving with W_ij
#endif
    kimax = ckmax.getValue();
    kimin = ckmin.getValue();
#ifdef WINMAT    
    kjmax = ckjmax.getValue();
    kjmin = ckjmin.getValue();
#endif

    /*==========================================================================*/
    /* read and save the input files                                            */
    /*==========================================================================*/
    std::vector<in_data> temp;    // temporary vector used to read the data
    int ddpk;   // dimention of the input data power spectrum. Then used as temporary storage
    int dkj;   // dimention of the kj(when used) and temporary storage
    /* data power spectrum */
    ddpk = input_data(Pkfn.getValue(), temp, 1);  //read the input data ps 
#ifndef WINMAT    
    kj = vec2gslvec(temp, &kimin, kimax, 0, false);  //and convert k into a gsl_vector
#endif
    data = vec2gslvec(temp, &kimin, kimax, 1, true);  //and convert P(k) into a gsl_vector
    temp.clear();   //clear temp for the covariance matrix 
    /* covariance matrix */
    dkj = input_data(covfn.getValue(), temp, 2);   //read the input covariance matrix
    // tranform the vector to a gsl matrix and invert it
    invcov = invert(temp, data->size, int(kimin), unbias.getValue(), false);
    temp.clear();   //clear temp for the kj
#ifdef WINMAT    
    /* kj */
    dkj = input_data(kjfn.getValue(), temp, 3);   //read the input kj
    kj = vec2gslvec(temp, &kjmin, kjmax, 1, true);     //and convert kj into a gsl_vector
    temp.clear();   //clear temp
    /* window matrix */
    window = winmat(wijfn.getValue(), int(kimin), int(kjmin), data->size, kj->size, ddpk, dkj);   //read and and save into a gsl_matrix
    /* W0j */
    dkj = input_data(w0jfn.getValue(), temp, 3);   //read the input W0j and save the dimention in dkj (used as temporary from here on)
    W0j = gsl_vector_alloc(kj->size);     //allocate the vector for W0j (has the same dimention as kj)
    for(i=0; i<int(kj->size); ++i) gsl_vector_set(W0j, i, temp[int(kjmin)+i].y);    //copy the interesting part of the W0j into the gsl vector
    temp.clear();   //clear temp
    /* G^2([0,ki]) */
    dkj = input_data(G2kifn.getValue(), temp, 3);   //read the input G^2([0,ki]) and save the dimention in dkj (used as temporary from here on)
    G20 = temp[0].y;   //save G^2(0) from the first place in G2kifn
    G2 = gsl_vector_alloc(int(data->size));     //allocate the vector for G2 (has the same dimention as the data)
    for(i=0; i<int(data->size); ++i) gsl_vector_set(G2, i, temp[int(kimin)+i+1].y);    //copy the interesting part of the G^2(ki) into the gsl vector
    temp.clear();   //clear temp
#endif

    /* read the linear and 1 loop power spectrum and save them into a spline */
    spllin = th2spl(linPk.getValue());
    acclin = gsl_interp_accel_alloc();    //allocate the accelerator for the linear power spectrum

    pmc *p22;   //pointer to the virtual base class used to compute the mode couple power spectrum
    if( dumpMC.getValue() == false )  //if the P_MC does not have to be dumped
      p22 = new pmc_interp( th2spl(mcPk.getValue()) );   //simple interpolation
    else  //otherwise 
      p22 = new pmc_dump( th2spl(mcPk.getValue()) );  //interpolation times dumping

    if( verbose.getValue() )
      std::cout << "Input files read" << std::endl;

    /*==========================================================================*/
    /* allocate the gsl vectors and matricies used to compute the likelihood    */
    /*==========================================================================*/
    theory = gsl_vector_alloc(kj->size);
    wtheory = gsl_vector_alloc(data->size);
    tempv = gsl_vector_alloc(data->size);
    tempm = gsl_matrix_alloc(data->size,data->size);
    tempm2 = gsl_matrix_alloc(data->size,data->size);

    /*==========================================================================*/
    /* set random number generators                                             */
    /*==========================================================================*/
    gsl_rng *r, *rg;   // gsl random number variables
    r = gsl_rng_alloc(gsl_rng_taus2);
    rg = gsl_rng_alloc(gsl_rng_taus2);
    gsl_rng_set(r, time(NULL));
    gsl_rng_set(rg, time(NULL)+1);

    /*==========================================================================*/
    /* do the MCMC                                                              */
    /*==========================================================================*/
    int max_step = max_step_cmd.getValue();  //number of step of the chain
    double likely_i, likely_i1;   // value of the likelihood in the first and in the second step
    double weight = 1.;  //  weight of each value for the parameters
    double *tempparams, *tempparams1;  //temporary store the values of the parameters

    /* open output file: parameters of the MCMC are written on fly */
    ofstream out(ofiles[0].c_str());
    out.setf(ios_base::scientific);
    out.precision(6);
    out.width(9);

    if( verbose.getValue() )
      std::cout << "Beginning the mcmc chains" << std::endl;

    /* loop with the Metropolis-Hasting algorithm*/
    for(i=0; i<max_step; ++i){
      for(j=0; j<n_params; ++j) for(;;){   //random walk from the old parameters and check if it is in the boundaries
	params1[j] = params[j] + gsl_ran_gaussian(rg, sigma[j]);
	if(params1[j] >= lower[j] && params1[j] <= upper[j]) break;
      }

      if(i == 0){
        p22->set_kstar( params[0] );  //set the value of kstar
	if( params[5] == -1 && sigma[5] == 0 && upper[5] == -1 && lower[5] == -1 )  //if sigma_v = -1
	  tempparams = substitute_value( params, n_params, 5, 1./params[0] );   //substitute with 1/kstar
	else tempparams = params;  //otherwise use the value from the random walk
	likely_i = likelihood(tempparams, p22);  //compute the likelihood
      }
      else if(i != 0 && weight == 1.) likely_i = likely_i1;

      p22->set_kstar( params1[0] );  //set the value of kstar
      if( params[5] == -1 && sigma[5] == 0 && upper[5] == -1 && lower[5] == -1 )  //if sigma_v = -1
	tempparams1 = substitute_value( params1, n_params, 5, 1./params1[0] );   //substitute with 1/kstar
      else tempparams1 = params1;  //otherwise use the random value
      likely_i1 = likelihood(tempparams1, p22);

      if( likely_i1 <= likely_i || exp(-0.5*(likely_i1 - likely_i)) >
	  gsl_rng_uniform(r) ){
	out << weight << "\t" << likely_i;   //print the weight and likelihood
	for(j=0; j<n_params; ++j){   //print the parameters and copy the new parameters into the old
	  out << "\t" << params[j];
	  params[j] = params1[j];
	}
	out << std::endl;
	weight = 1.;
      }
      else if(i == max_step - 1){   //save the last step
	out << weight << "\t" << likely_i;   //print the weight and likelihood
	for(j=0; j<n_params; ++j){   //print the parameters and copy the new parameters into the old
	  out << "\t" << params[j];
	  params[j] = params1[j];
	}
      }
      else ++weight;

      if((i+1)%100000 == 0 and verbose.getValue() ) 
	std::cout << "loop number: " << i+1 << std::endl;
    }

    out.close();
    /*==========================================================================*/
    /* deallocations                                                            */
    /*==========================================================================*/
    /* deallocate the parameters */
    delete[] params;
    delete[] params1;
    delete[] upper; 
    delete[] lower; 
    delete[] sigma; 
    delete p22;
    /* deallocate the random */
    gsl_rng_free(r);
    gsl_rng_free(rg);
    /* deallocate common variable */
    dealloc_comm();
  } 
  catch (TCLAP::ArgException &e){  // catch any exceptions
    std::cerr << "Error: " << e.error() << " for arg " << e.argId() << std::endl; 
    exit(10);
  }
  catch (inputexception &e){
    std::cerr << "Error reading input file '" << e.what() << "' (Probably it does not exists)." << std::endl;
    exit(5);
  }
  catch (outputexception &e){
    std::cerr << "The output file '" << e.what() << "' already exists. ";
    std::cerr << "Delete it, move it, rename it or change the output root." << std::endl;
    exit(6);
  }

  return 0;
}

