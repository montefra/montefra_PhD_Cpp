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
  std::string dataset_file;
  if(ini.get_param(pk_dataset, &dataset_file) != 0)
    ini.ini_error(pk_dataset, 20);
  if(common::verbose)
    std::cout << "Reading dataset #" << dataset << " from " << dataset_file << std::endl;
  data = new Dataset(dataset_file);

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

  //if same_alpha==true all likelihoods are computed using the same alpha
  //if same_params==true all likelihoods are computed using the same parameters 
  bool same_alpha, same_params;
  if(ini->get_param("same_alpha", &same_alpha) != 0)
    same_alpha = false;
  if(ini->get_param("same_params", &same_params) != 0)
    same_params = false;

  int tdataset=dataset;  //local dataset number
  if(same_params) tdataset = 1;  //use 1 for every likelihood

  std::vector<std::string>::iterator it;
  std::string temp;

  if(same_alpha){  //look for the value called alpha and use dataset=1 only for that
    for(it=paramnames.begin(); it!=paramnames.end(); ++it){
      //add the dataset number
      if((*it).compare("alpha")==0)
        temp = *it+"1";
      else
        temp = *it+common::to_string(tdataset); 
      mcmc2theory[*it] = temp;  //add the element to the map
      *it = temp;  //save the new paramname into the vector
    }
  }
  else{
    for(it=paramnames.begin(); it!=paramnames.end(); ++it){
      temp = *it+common::to_string(tdataset); 
      mcmc2theory[*it] = temp;  //add the element to the map
      *it = temp;  //save the new paramname into the vector
    }
  }
}

/*==========================================================================
 * return the likelihood with the given parameters
 * Parameters
 * ----------
 *  params: map
 *    map with parameter name as key and its value as value
 *==========================================================================*/
double Likelihood::get_chi2(std::map<std::string, double> params){
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

  return(chi2);
}
