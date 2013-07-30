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
MCMC::MCMC(std::vector<Likelihood> &likes, ParseIni &ini){
  this->likes = likes;
  n_likes=likes.size();

  //read the file root
  if(ini.get_param("file_root", &file_root)!=0)
    ini.ini_error("file_root", 60);
  //read the number of steps of the chain
  if(ini.get_param("n_steps", &n_steps)!=0){
    size_t def_n_step = 300000;
    std::cerr << "n_steps not found of found empty in the inifile ";
    std::cerr << ini.get_fname() << ". Defaulting to " << def_n_step << std::endl;
    n_steps = def_n_step;
  }

  init_random();
  //get the paramnames from the likelihoods
  //save them to file and read the parameters from the ini
  initialise_paramnames(ini);
}

/*==========================================================================
 * initialiase the random generators
 *==========================================================================*/
void MCMC::init_random(){
  r = gsl_rng_alloc(gsl_rng_taus2);
  rg = gsl_rng_alloc(gsl_rng_taus2);
  gsl_rng_set(r, time(NULL));
  gsl_rng_set(rg, time(NULL)+1);
}

/*=========================================================================a
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
  std::string ofile = file_root+".paramnames";
  if(common::fileexists(ofile)){
    std::cerr << "File " << ofile << " already exists";
    std::cerr << "Delete, move or rename it" << std::endl;
    exit(61);
  }

  std::ofstream out(ofile.c_str());
  if(!out.is_open()){
    std::cerr << "Failed to open output file " << ofile << std::endl;
    exit(62);
  }
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
void MCMC::read_params(ParseIni &ini){
  double tstart, tmin, tmax, tstep; //temporary doubles
  //read the parameter limits one by one and save into the maps
  std::vector<std::string>::iterator it;
  for(it = paramnames.begin(); it != paramnames.end(); ++it){
    if(ini.get_mcmc_params(*it, &tstart, &tmin, &tmax, &tstep)!=0)
      ini.ini_error(*it, 63);
    //if the parameters is to be changed (step size >0), start at a randomly
    //displaced point within the range
    if(tstep > 0.) tstart = tstart + gsl_ran_gaussian(rg, tstep); 
    if(tstart<tmin) tstart=tmin;
    else if(tstart>tmax) tstart=tmax;
    start[*it] = tstart;
    min[*it] = tmin;
    max[*it] = tmax;
    step[*it] = tstep;
  }
}

/*==========================================================================
 * initialise the parameter names (mostly call the three above functions)
 *==========================================================================*/
void MCMC::initialise_paramnames(ParseIni &ini){
  get_paramnames();
  save_paramnames();
  long_names.clear();  //clear the long names
  read_params(ini);  //get the paramnames limits and steps from the inifile
}

/*==========================================================================
 * run the mcmc chain
 *==========================================================================*/
void MCMC::run(){
  //open the output file
  std::string ofile = file_root+".txt";
  if(common::fileexists(ofile)){
    std::cerr << "File " << ofile << " already exists";
    std::cerr << "Delete, move or rename it" << std::endl;
    exit(64);
  }

  std::ofstream out(ofile.c_str());
  if(!out.is_open()){
    std::cerr << "Failed to open output file " << ofile << std::endl;
    exit(65);
  }
  out.setf(std::ios_base::scientific);
  out.precision(6);
  out.width(9);

  //likelihoods
  double like1=1., like2=1.;
  size_t weights=0; //number of times the chains stays in one point
  //compute the likelihood at the starting point
  std::map<std::string, double> params1(start);
  for(size_t i=0; i<n_likes; ++i) like1 *= likes[i].get_like(params1);
  weights=1;

  for(size_t j=1; j<n_steps; ++j){ //loop over the mcmc steps
    std::map<std::string, double> params2 = new_parameters(params1);
    for(size_t i=0; i<n_likes; ++i) like2 *= likes[i].get_like(params2);

    //criteria to accept the new likelihood: eitheir is larger or the 
    //new-old ratio is larger than a random number. In this case print
    //the old weight, likelihood and parameters, set the weights to 1 
    //and save params2 and like2 into params1 and like1
    if(like1 < like2 || like1/like2 > gsl_rng_uniform(r)){
      out << weights << "\t" << like1;
      //loop over the paramnames vector to make sure that the parameters are
      //printed always in the same order
      for(std::vector<std::string>::iterator it=paramnames.begin();
          it!=paramnames.end(); ++it){  
        out << "\t" << params1[*it];
        params1[*it] = params2[*it];
      }
      out << std::endl;
      weights = 1;
      like1 = like2;
    }
    else ++weights;

    if((j+1)%100000 == 0 and common::verbose) 
      std::cout << "loop number: " << j+1 << std::endl;
  }

  //if weights larger than one, means that the last set of parameter of the
  //chain has not been saved to file
  if(weights>1){
    out << weights << "\t" << like1;
    //loop over the paramnames vector to make sure that the parameters are
    //printed always in the same order
    for(std::vector<std::string>::iterator it=paramnames.begin();
        it!=paramnames.end(); ++it){  
      out << "\t" << params1[*it];
    }
    out << std::endl;
  }
  out.close();
}

/*==========================================================================
 * randomly generate a new set of parameters
 * Parameters
 * ----------
 *  params: map<string, double>
 *    map containing the parameters to substitute
 * output
 * ------
 *  newparams: map<string, double>   
 *    map with the new parameters
 *==========================================================================*/
std::map<std::string, double> MCMC::new_parameters(std::map<std::string,
    double> &params){
  std::map<std::string, double> newparams;
  std::map<std::string, double>::iterator it; 
  for(it=params.begin(); it!=params.end(); ++it){
    for(;;){
      newparams[it->first] = it->second + gsl_ran_gaussian(rg, step[it->first]);
      if(newparams[it->first] >= min[it->first] && 
          newparams[it->first] <= max[it->first]) break;
    }
  }
  return newparams;
}

