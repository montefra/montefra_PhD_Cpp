/*==========================================================================
 * Version: 5.00           date:xx/xx/xx                                    
 *                                                                          
 * Author: Francesco Montesano, MPE, Garching                               
 *                                                                          
 * Purpose: likelihood class. Basically contains a data object and a theory
 * object and take care of connecting those two to the mcmc engine
 *==========================================================================*/

#include "likelihood.h"

/*==========================================================================
 * Constructor
 * Parameters
 * ----------
 *  dataset: int
 *    dataset number
 *  ini: ParseIni object
 *    container of the inifile
 *  theory: Theory object
 *    used to get the model
 *==========================================================================*/
Likelihood::Likelihood(int dataset, ParseIni &ini, Theory &theory)
{
  this->dataset=dataset;
  this->ini=&ini;
  this->theory=&theory;
  //read the dataset
  std::string pk_dataset = std::string("pk_dataset")+
    common::to_string(this->dataset);
  if(common::verbose)
    std::cout << "Reading dataset :" << dataset << std::endl;
  data = new Dataset(pk_dataset);

  //get the paramnames from the theory and create append the dataset number to
  //the short name
  retrieve_paramnames();

  //get the values of k where to evaluate the power spectrum and allocate it
  k = data->get_k();
  n_k = k->size;
  //allocate pktheory of the same size of the k
  pktheory = gsl_vector_alloc(n_k); 
}

/*==========================================================================
 * get the paramnames from the theory object, append the dataset number 
 * to the paramnames and save them
 *==========================================================================*/
void Likelihood::retrieve_paramnames(){
  //get the paramnames from the theory
  paramnames = theory->get_paramnames();
  long_names = theory->get_long_names();

  for(size_t i=0; i<paramnames.size(); ++i){
    std::string temp = paramnames[i]+common::to_string(dataset); //add the dataset number
    mcmc2theory[paramnames[i]] = temp;  //add the element to the map
    paramnames[i] = temp;  //save the new paramname into the vector
  }
}

/*==========================================================================
 * return the likelihood with the given parameters
 * Parameters
 * ----------
 *  params: map
 *    map with parameter name as key and its value as value
 *==========================================================================*/
double Likelihood::get_like(std::map<std::string, double> &params){
  std::map<std::string, double> to_theory;  //parameters to pass to the theory
  std::map<std::string, std::string>::iterator it; //interator for the mcmc2theory map
  //save in to_theory the parameters to pass to theory.get_model to compute the
  //model
  for(it=mcmc2theory.begin(); it!=mcmc2theory.end(); ++it)
    to_theory[it->first] = params[it->second];

  //get the model power spectrum give the parameters
  for(size_t i=0; i<n_k; ++i) 
    gsl_vector_set(pktheory, i, theory->get_model(gsl_vector_get(k, i), to_theory));

  //convolve the model power spectrum
  pkconvolved = data->convolve(pktheory);

  //compute the chi^2
  double chi2 = data->get_chi2(pkconvolved);

  return(exp(-0.5*chi2));
}
