/*==========================================================================
 * Version: 5.00           date:xx/xx/xx                                    
 *                                                                          
 * Author: Francesco Montesano, MPE, Garching                               
 *                                                                          
 * Purpose: read the pk data to use in the mcmc chain
 *==========================================================================*/

#include "read_file.h"

/*==========================================================================
 * Constructor of the object
 * Parameter
 * ---------
 *  pk_file: string with the file name containing the name of the file to
 *  read in and the bins to consider
 *==========================================================================*/
Dataset::Dataset(std::string &pk_dataset){
  
  //import the dataset file into a inifile object
  dataset = new ParseIni(pk_dataset);
  //get the number of bins
  get_bin_numbers(*dataset);
  //read the measured power spectrum
  read_pk(*dataset);
  //read and invert the covariance matrix
  read_invert_cov(*dataset);
  //read the window file files
  read_window(*dataset);

  //allocata the gsl vectors and matrices needed for the window matrix
  //convolution and chi^2
  alloc_gsl();
}

/*==========================================================================
 * Destructor of the object
 *==========================================================================*/
Dataset::~Dataset(){
  gsl_vector_free(data);
  gsl_matrix_free(invcov);
  gsl_matrix_free(Wij);
  gsl_vector_free(kj);
  gsl_vector_free(W0j);
  gsl_vector_free(G2i);
  gsl_vector_free(convolved_model);
  gsl_vector_free(tempv);
}

/*==========================================================================
 * read the number of bins to read from the files and to use for the
 * likelihood
 * Parameter
 * ---------
 *  ini: inifile
 *==========================================================================*/
void Dataset::get_bin_numbers(ParseIni ini){
  if(ini.get_param("kitot", &n_bins.kitot) != 0)
    ini.ini_error("kitot", 30);
  if(ini.get_param("kimin", &n_bins.kimin) != 0)
    ini.ini_error("kimin", 31);
  if(ini.get_param("kimax", &n_bins.kimax) != 0)
    ini.ini_error("kimax", 32);
  if(ini.get_param("kjtot", &n_bins.kjtot) != 0)
    ini.ini_error("kjtot", 33);
  if(ini.get_param("kjmin", &n_bins.kjmin) != 0)
    ini.ini_error("kjmin", 34);
  if(ini.get_param("kjmax", &n_bins.kjmax) != 0)
    ini.ini_error("kjmax", 35);

  n_bins.kiuse = n_bins.kimax-n_bins.kimin;
  n_bins.kjuse = n_bins.kjmax-n_bins.kjmin;
}

/*==========================================================================
 * read measured power
 * Parameter
 * ---------
 *  ini: inifile
 *==========================================================================*/
void Dataset::read_pk(ParseIni ini){
  std::string pk_file;
  if(ini.get_param("pk_file", &pk_file)!=0)
    ini.ini_error("pk_file", 36);
  if(common::verbose)
    std::cout << "Reading power spectrum file " << pk_file << std::endl;

  double dummy, pk;
  //number of skipped lines and lines saved into the gsl_vector
  int skip_counter=0, saved_lines=0;

  //allocate the gls vector that will contain the data
  data = gsl_vector_alloc(n_bins.kiuse);   

  std::ifstream in(pk_file.c_str(), std::ifstream::in);
  std::string line;   //string containing the line
  while(getline(in, line)){   //read a line of the input file
    if(line.find("#") == 0)  //if in the line there is a '#' skip the line (it's likely a comment)
      continue;
    if(skip_counter++<n_bins.kimin) continue;  //skip first kimin bins
    else if(saved_lines < n_bins.kiuse){  //save all the kiuse after kimin
      std::stringstream(line) >> dummy >> pk >> dummy >> dummy;
      gsl_vector_set(data, saved_lines++, pk);
    }
    else break;  //after reading all of them break out of the while loop
  }
  in.close();
}

/*==========================================================================
 * read and invert the covariance matrix
 * Parameter
 * ---------
 *  ini: inifile
 *==========================================================================*/
void Dataset::read_invert_cov(ParseIni ini){
  std::string covmat_file; //file name of the covariance
  if(ini.get_param("covmat_file", &covmat_file)!=0)
    ini.ini_error("covmat_file", 37);
  if(common::verbose){
    std::cout << "Reading and inverting the convariance from file ";
    std::cout << covmat_file << std::endl;
  }

  //read and cut the covariance matrix
  gsl_matrix * temp_cov = gslf::read_cut_gsl_matrix(covmat_file, n_bins.kitot,
      n_bins.kitot, n_bins.kiuse, n_bins.kiuse, n_bins.kimin, n_bins.kimin);
  //allocata the space for the inverse
  invcov = gsl_matrix_alloc(n_bins.kiuse,n_bins.kiuse);  

  //invert the covariance using LU decomposiation
  // permutation for LU decomposition
  gsl_permutation *p = gsl_permutation_alloc(n_bins.kiuse);
  int signum;    // required by LU decomposition
  if(gsl_linalg_LU_decomp(temp_cov, p, &signum)!= 0){
    std::cerr << "LU decomposition failed" << std::endl;
    exit(38);
  }
  if(gsl_linalg_LU_invert(temp_cov, p, invcov)!= 0){    // inversion
    std::cerr << "Matrix inversion failed" << std::endl;
    exit(39);
  }
  gsl_permutation_free(p);
  gsl_matrix_free(temp_cov);   //free the memory of the temporary matrix

  // unbias inverse covariance matrix
  unbias_inverse_cov(ini); 
}
/*==========================================================================
 * Unbias the inverse covariance
 * Parameter
 * ---------
 *  ini: inifile
 *==========================================================================*/
void Dataset::unbias_inverse_cov(ParseIni ini){
  // unbias the inverse covariance matrix as in Hartlap et al. 2007 
  // multiplying it my 1-D = 1 - (n_bins+1)/(n_mocks-1)
  // If a correlation 'r_mocks' is given and positive, D is multiplied by
  // (1+r_mocks^2)/2 according to Percival et al 2013

  int n_mocks;  //number of mocks used to estimate che covariance matrix
  double r_mocks;   //correlation between mocks
  int err_n_mocks, err_r_mocks;  //errors where reading the inifile

  // read the number of mocks and correlation
  err_n_mocks = ini.get_param("n_mocks", &n_mocks);
  err_r_mocks = ini.get_param("r_mocks", &r_mocks);

  if(err_r_mocks==0 && n_mocks>n_bins.kiuse+2){
    double D = (n_bins.kiuse + 1.)/(n_mocks-1.);
    if(err_r_mocks==0 && r_mocks>0) // if the correlation is given and positive
      D *= (1.+pow(r_mocks, 2))/2.;
    if(gsl_matrix_scale(invcov, 1.-D) != 0){
      std::cerr << "Error multiplying the inverse covariance matrix to unbias it";
      std::cerr << std::endl;
      exit(41);
    }
  }
  else{
    std::cerr << "The number of mocks has not been provided or is not invertible";
    std::cerr << "The inverse covariance matrix is biased if computed from simulations." << std::endl;
  }
}

/*==========================================================================
 * read the files for the window matrix
 * Parameter
 * ---------
 *  ini: inifile
 *==========================================================================*/
void Dataset::read_window(ParseIni ini){
  std::string file_names; //window matrix file names
 
  if(common::verbose)
    std::cout << "Reading the window matrix files" << std::endl;
  //read the window matrix Wij
  if(ini.get_param("wij_file", &file_names)!=0)
    ini.ini_error("wij_file", 42);
  Wij = gslf::read_cut_gsl_matrix(file_names, n_bins.kitot,
      n_bins.kjtot, n_bins.kiuse, n_bins.kjuse, n_bins.kimin, n_bins.kjmin);

  //read kj
  if(ini.get_param("kj_file", &file_names)!=0)
    ini.ini_error("kj_file", 43);
  kj = gslf::read_cut_gsl_vector(file_names, n_bins.kjtot, n_bins.kjuse, n_bins.kjmin);

  //read W0j
  if(ini.get_param("w0j_file", &file_names)!=0)
    ini.ini_error("w0j_file", 44);
  W0j = gslf::read_cut_gsl_vector(file_names, n_bins.kjtot, n_bins.kjuse, n_bins.kjmin);

  //read whole G20i file
  if(ini.get_param("G2i_file", &file_names)!=0)
    ini.ini_error("G2i_file", 45);
  gsl_vector *temp_G = gslf::read_gsl_vector(file_names, n_bins.kitot+1);
  G20 = gsl_vector_get(temp_G, 0);  //get the first element
  G2i = gsl_vector_alloc(n_bins.kiuse);
  //copy the interesting part of the vector
  gsl_vector_view vview = gsl_vector_subvector(temp_G, n_bins.kimin+1,
      n_bins.kiuse);
  gsl_vector_memcpy(G2i, &(vview.vector));
  gsl_vector_free(temp_G);
}

/*==========================================================================
 * Allocate gls vectors and matrices used when convolving and computing the
 * chi^2
 *==========================================================================*/
void Dataset::alloc_gsl(){
  convolved_model = gsl_vector_alloc(n_bins.kiuse);
  tempv = gsl_vector_alloc(n_bins.kiuse);
}

/*==========================================================================
 * return the input gsl vector convolved with the window matrix
 * Parameters
 * ----------
 *  model: gsl vector
 *    power spectrum to be convolved
 * output
 * ------
 *  convolved_model: gsl vector
 *==========================================================================*/
gsl_vector *Dataset::convolve(gsl_vector *model){
  if(gsl_blas_dgemv(CblasNoTrans, 1., Wij, model, 0., convolved_model) !=0) 
    exit(46);
  // integral constraint: add a rescaled version of the window function to
  // convolved_model
  double c;  //constant used to rescale the window function
  if(gsl_blas_ddot(W0j, model, &c)) exit(47);  // c = sum_j W0j P(kj)
  // convolved_model = convolved_model - c*G^2(ki)/G^2(0)
  if(gsl_blas_daxpy(-c/G20, G2i, convolved_model)) exit(48);  

  return(convolved_model);
}

/*==========================================================================
 * compute the chi^2 give the model
 * Parameters
 * ----------
 *  model: gsl vector
 *    model power spectrum
 * output
 * ------
 *  chi2: double
 *    chi^2 from the data stored here and a model
 *==========================================================================*/
double Dataset::get_chi2(gsl_vector *model){
  double chi2;
  if(gsl_vector_sub(model, data) != 0) exit(49);  //model = model - data
  if(gsl_blas_dgemv(CblasNoTrans, 1., invcov, model, 0., tempv)!= 0)
    exit(50);  // tempv=invcov (model-data)  
  // chi^2=(model-data)^{T}tempv
  if(gsl_blas_ddot(model, tempv, &chi2)!= 0 ) exit(51);

  return(chi2);
}
