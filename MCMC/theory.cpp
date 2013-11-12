/*==========================================================================
 * Version: 5.00           date:xx/xx/xx                                    
 *                                                                          
 * Author: Francesco Montesano, MPE, Garching                               
 *                                                                          
 * Purpose: read the linear and 1loop power spectrum and save them into 
 * gsl spline objects. 
 * Provide method(s) to that return the model
 *==========================================================================*/

#include "theory.h"

/*==========================================================================
 * Constructor
 *
 * Parameters
 * ----------
 *  ini: ParseIni object
 *    ini parser objects containing. This code need that ini contains 'plin'
 *    and 'p1loop' files
 *==========================================================================*/
Theory::Theory(ParseIni ini){
  std::string file_names; //names of the model power spectra

  //linear power spectrum
  if(ini.get_param("plin", &file_names)!=0)
    ini.ini_error("plin", 10);
  spllin = gslf::read_2_spline(file_names);
  //1loop power spectrum
  if(ini.get_param("p1loop", &file_names)!=0)
    ini.ini_error("p1loop", 11);
  spl1l = gslf::read_2_spline(file_names);

  //allocate the accelerators
  acclin = gsl_interp_accel_alloc();
  acc1l = gsl_interp_accel_alloc();

  //fill the paramnames
  fill_paramnames();
}
/*==========================================================================
 * Destructor
 * deallocate stuff
 *==========================================================================*/
Theory::~Theory(){
  //free splines
  gsl_spline_free(spllin);
  gsl_spline_free(spl1l);
  //free the interpolation accelerators
  gsl_interp_accel_free(acclin);
  gsl_interp_accel_free(acc1l);
}

/*==========================================================================
 * fill the paramnames with the names of the parameters used here
 * the order in the paramnames vector is the same as the one accepted by
 * 'get_model'
 *==========================================================================*/
void Theory::fill_paramnames(){
    paramnames.push_back("k_star");
    long_names.push_back("k_\\star [h/Mpc]");
    paramnames.push_back("amc");
    long_names.push_back("A_{\\mathrm{MC}}");
    paramnames.push_back("alpha");
    long_names.push_back("\\alpha");
    paramnames.push_back("bias"); 
    long_names.push_back("b");
    paramnames.push_back("noise");
    long_names.push_back("n [(Mpc/h)^3]");
}

/*==========================================================================
 * compute the model in wavenumber k given the model parameter in double
 * array. The parameter names are defined in 'fill_paramnames()'
 * Parameters
 * ----------
 *  k: double
 *    wavenumber where to evaluate the model
 *  params: map
 *    map containing the parameter name as key and its value as value
 * output
 * ------
 *  pk: double
 *    model power spectrum evaluated in k
 *==========================================================================*/
double Theory::get_model(double k, std::map<std::string, double> params){
  double ak = params["alpha"]*k;
  double pk = gsl_spline_eval(spllin, ak, acclin);   // P_lin(alpha*k)
  pk *= exp(-pow(k/params["k_star"], 2));  //*exp(-(k/kstar)^2)
  //+ A_MC P_1loop(alpha*k)
  pk += params["amc"] * gsl_spline_eval(spl1l, ak, acc1l);
  pk *= pow(params["bias"], 2);  // *bias^2
  pk += params["noise"];   //add noise
  pk /= pow(params["alpha"], 3);  // /alpha^3: so not absorbed into bias and noise
  return(pk);                     
}
