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
  if(argc < 1){
    std::cerr << "No inifile provided." << std::endl;
    std::cerr << "Usage: " << std::endl;
    std::cerr << "  " << argv[0] << " inifile.ini" << std::endl;
    exit(1);
  }

  ParseIni ini_file(argv[1]);

  //initialise the theory reading the linear and 1loop power spectrum
  Theory theory(ini_file);

  int n_datasets; //number of datasets
  if(ini_file.get_parami<int>("n_datasets", &n_datasets) != 0){
    std::cerr << "'n_datasets' is not in the inifile" << std::endl;
    exit(2);
  }

  //vector of likelihoods
  std::vector<Likelihood> likelihoods; 
  //create a likelihood for each dataset
  //the constructor wants the likelihood number (starting from 1)
  //the inifile and the theory object (it saves a local copy) 
  for(int i=0; i<n_dataset; ++i){
    likelihoods.push_back(Likelihood(i+1, ini_file, theory));
  }

  //initialise mcmc engine creating the parameter names and reading the
  //parameters limits
  MCMC mcmc(likelihoods, ini_file);

  mcmc.run();

  exit(0);
}

