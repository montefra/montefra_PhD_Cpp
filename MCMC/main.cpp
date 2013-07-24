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

  int n_datasets; //number of datasets
  if(ini_file.get_parami<int>("n_dataset", &n_dataset) != 0){
    std::cerr << "'n_dataset' is not in the inifile" << std::endl;
    exit(2);
  }

  //vectors of likelihoods.
  //Upcasting to a base class
  std::vector<Read_pk_files> pk_data; 

  for(int i=0; i<n_dataset; ++i){
    std::string pk_dataset = std::string("pk_dataset")+to_string(i+1);  //dataset file name
    pk_data.pop_back(

  }


  //data objects
  //likelihood objec that uses the data
  //mcmc that uses the likelihood

  exit();
}

