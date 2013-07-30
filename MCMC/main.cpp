/*==========================================================================
 * Version: 5.00           date:xx/xx/xx                                    
 *                                                                          
 * Author: Francesco Montesano, MPE, Garching                               
 *                                                                          
 * Purpose: Markov Chain Monte Carlo
 *==========================================================================*/

#include "main.h"

int main(int argc, char* argv[])
{
  // read the inifile
  if(argc < 2){
    std::cerr << "No inifile provided." << std::endl;
    std::cerr << "Usage: " << std::endl;
    std::cerr << "  " << argv[0] << " inifile.ini" << std::endl;
    exit(1);
  }

  ParseIni ini_file(argv[1]);

  //verbose: set to false by default
  if(ini_file.get_param("verbose", &common::verbose)!=0) 
    common::verbose = false;

  //initialise the theory reading the linear and 1loop power spectrum
  Theory theory(ini_file);

  int n_datasets; //number of datasets
  if(ini_file.get_param("n_datasets", &n_datasets) != 0)
    ini_file.ini_error("n_datasets", 2);
  if(common::verbose) 
    std::cout << "Reading " << n_datasets << " datasets" << std::endl;

  //vector of likelihoods
  std::vector<Likelihood> likelihoods; 
  //create a likelihood for each dataset
  //the constructor wants the likelihood number (starting from 1)
  //the inifile and the theory object (it saves a local copy) 
  for(int i=0; i<n_datasets; ++i){
    likelihoods.push_back(Likelihood(i+1, ini_file, theory));
  }
  if(common::verbose){
    std::cout << "Datasets read and likelihoods initialized" << std::endl;
    std::cout << std::endl << "Start the mcmc chain" << std::endl;
  }

  //initialise mcmc engine creating the parameter names and reading the
  //parameters limits
  MCMC mcmc(likelihoods, ini_file);

  mcmc.run();

  exit(0);
}

