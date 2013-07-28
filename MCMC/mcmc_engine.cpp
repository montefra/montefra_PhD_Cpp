/*==========================================================================
 * Version: 5.00           date:xx/xx/xx                                    
 *                                                                          
 * Author: Francesco Montesano, MPE, Garching                               
 *                                                                          
 * Purpose: initialise the parameter and run the mcmc 
 *==========================================================================*/

#include "mcmc_engine.h"

/*==========================================================================
 * Constructor
 * get the number of datasets used and the inifile objects containing the
 * parameter limits
 * Parameters
 * ----------
 *  n_datasets: int
 *    number of datasets
 *  ini: ParseIni object
 *    ini object containing the the parametes limits and step
 *==========================================================================*/
MCMC::MCMC(std::vector<Likelihood> likes, ParseIni ini):
  this->likes = likes;
{
  n_likes=likes.size();

  //read the file root
  if(ini.get_param<std::string>("file_root", &file_root)!=0) exit(70);
  //read the number of steps of the chain
  if(ini.get_param<size_t>("n_steps", &n_steps)!=0) exit(71);

  init_random();
  //get the paramnames from the likelihoods
  //save them to file and read the parameters from the ini
  initialise_paramnames(ini);
}


/*==========================================================================
 * retrieve the parameter names from the likelihoods
 *==========================================================================*/
void MCMC::get_paramnames(){
  //get all the paramnames from the first likelihood
  paramnames = likes[0].get_paramnames();
  long_names = likes[0].get_long_names();
  //get the paramnames from the other likelihoods
  //if any already exists skip it
  for(size_t i=1; i<likes.size(); ++i){
    std::vector<std::string> temp_pn(likes[i].get_paramnames());
    std::vector<std::string> temp_ln(likes[i].get_long_names());
    for(size_t j=0; j<temp_pn.size(); ++j){
      if(std::find(paramnames.begin(), paramnames.end(), temp_pn[j]) !=
          paramnames.end()) continue;
      else{
        paramnames.push_back(temp_pn[j]);
        long_names.push_back(temp_ln[j]);
      }
    }
  }
  n_parameters = paramnames.size();  //save the number of parameters
}

/*==========================================================================
 * create the short and long paramnames, save them in the paramnames file
 *==========================================================================*/
void MCMC::save_paramnames(){
  std::string ofile = file_root+".paramnames"
  ofstream out(ofile.c_str());
  for(size_t i=0; i<n_parameters; ++i) 
    out << paramnames[i] << "\t" << long_names[i] << std::endl;
  out.close();
}

/*==========================================================================
 * read the paramnames limits and steps
 * Parameters
 * ----------
 *  ini: ParseIni object
 *    inifile containing the parameters 
 *==========================================================================*/
void MCMC::read_params(ParseIni ini){
  double tstart, tmin, tmax, tstep; //temporary doubles
  //read the parameter limits one by one and save into the maps
  std::vector<std::string>::iterator it;
  for(it = paramnames.begin(); it != paramnames.end(); ++it){
    if(get_mcmc_params(*it, &tstart, &tmin, &tmax, &tstep)!=0){
      std::cerr << "Parameter " << *it;
      std::cerr << " not found or incomplete" << std::endl;
      exit(75);
    }
    //if the parameters is to be changed (step size >0), start at a random point into the
    //allowed parameter range
    if(tstep > 0.) tstart = tmin + (tmax-tmin)*gsl_rng_uniform(r); 
    start[*it] = tstart;
    min[*it] = tmin;
    max[*it] = tmax;
    step[*it] = tstep;
  }
}

/*==========================================================================
 * run the mcmc chain
 *==========================================================================*/
void MCMC::run(){
  //open the output file
  std::string ofile = file_root+".txt"
  ofstream out(ofile.c_str());
  out.setf(ios_base::scientific);
  out.precision(6);
  out.width(9);

  //likelihoods
  double like1=0, like2=0;
  //compute the likelihood at the starting point
  std::map<std::string, double> params1(start);
  for(size_t i=0; i<n_likes; ++i) like1 += likes[i].get_like(params1);



  close(out)
}
