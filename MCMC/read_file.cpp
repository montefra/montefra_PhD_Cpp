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
Read_pk_files::Read_pk_files(std::string pk_dataset){
  
  //import the dataset file into a inifile object
  dataset = new ParseIni(pk_dataset);
  //get the number of bins
  get_bin_numbers(dataset);
  //read the measured power spectrum
  read_pk(dataset);
  //read and invert the covariance matrix
  read_invert_cov(dataset);
  //read the window file files
  read_window(dataset);

  //allocata the gsl vectors and matrices needed for the window matrix
  //convolution and chi^2
  alloc_gsl();
}


/*==========================================================================
 * read the number of bins to read from the files and to use for the
 * likelihood
 * Parameter
 * ---------
 *  ini: inifile
 *==========================================================================*/
void Read_pk_files::get_bin_numbers(ParseIni ini){
    if(ini.get_param<int>("kitot", &n_bins.kitot) != 0) exit(10);
    if(ini.get_param<int>("kimin", &n_bins.kimin) != 0) exit(11);
    if(ini.get_param<int>("kimax", &n_bins.kimax) != 0) exit(12);
    if(ini.get_param<int>("kjtot", &n_bins.kjtot) != 0) exit(13);
    if(ini.get_param<int>("kjmin", &n_bins.kjmin) != 0) exit(14);
    if(ini.get_param<int>("kjmax", &n_bins.kjmax) != 0) exit(15);

    n_bins.kiuse = n_bins.kimax-n_bins.kimin;
    n_bins.kjuse = n_bins.kjmax-n_bins.kjmin;
}

/*==========================================================================
 * read measured power
 * Parameter
 * ---------
 *  ini: inifile
 *==========================================================================*/
void Read_pk_files::read_pk(ParseIni ini){
  std::string pk_file;
  if(ini.get_param<std::string>("pk_file", &pk_file)!=0) exit(20);

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
    else if(saved_counter < n_bins.kiuse){  //save all the kiuse after kimin
      std::stringstream(line) >> dummy >> pk >> dummy >> dummy;
      gsl_vector_set(data, saved_counter++, pk);
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
void Read_pk_files::read_invert_cov(ParseIni ini){
  std::string covmat_file; //file name of the covariance
  if(ini.get_param<std::string>("covmat_file", &covmat_file)!=0) exit(21);
  //read and cut the covariance matrix
  gsl_matrix * temp_cov = read_cut_matrix(covmat_file, n_bins.kitot,
      n_bins.kitot, n_bins.kiuse, n_bins.kiuse, n_bins.kimin, n_bins.kimin);
  //allocata the space for the inverse
  invcov = gsl_matrix_alloc(n_bins.kiuse,n_bins.kiuse);  

  //invert the covariance using LU decomposiation
  // permutation for LU decomposition
  gsl_permutation *p = gsl_permutation_alloc(n_bins.kiuse);
  int signum;    // required by LU decomposition
  if(gsl_linalg_LU_decomp(temp_cov, p, &signum)!= 0){
    std::cerr << "LU decomposition failed" << std::endl;
    exit(22);
  }
  if(gsl_linalg_LU_invert(temp_cov, p, invcov)!= 0){    // inversion
    std::cerr << "Inversion failed" << std::endl;
    exit(23);
  }
  gsl_permutation_free(p);
  gsl_matrix_free(temp_cov);   //free the memory of the temporary matrix

  // unbias the inverse covariance matrix as in Hartlap et al. 2007 
  // multiplying it my (n_mocks-n_bins-2)/(n_mock-1)
  int n_mocks;  //number of mocks used to estimate che covariance matrix
  if(ini.get_param<int>("n_mocks", &n_mocks) == 0){
    // if n-p-2<0 it's not possible to unbias the inverse
    if(n_mocks < n_bins.kiuse+2){
      std::cerr << "The covariance matrix is singular and should not be ";
      std::cerr << "possible to invert it" << std::endl;
      exit(24);
    }
    else{
      if(gsl_matrix_scale(invcov, (n_mocks-n_bins.kiuse-2.)/(n_mocks-1.)) != 0){
        std::cerr << "Error multiplying the inverse covariance matrix to unbias it";
        std::cerr << std::endl;
        exit(25);
      }
    }
  }
  else{
    std::cerr << "The number of mocks has not been provided. ";
    std::cerr << "The inverse covariance matrix is biased." << std::endl;
  }
}

/*==========================================================================
 * read the files for the window matrix
 * Parameter
 * ---------
 *  ini: inifile
 *==========================================================================*/
void Read_pk_files::read_window(ParseIni ini){
  std::string file_names; //window matrix file names
 
  //read the window matrix Wij
  if(ini.get_param<std::string>("wij_file", &file_names)!=0) exit(30);
  Wij = read_cut_matrix(wij_file, n_bins.kitot,
      n_bins.kjtot, n_bins.kiuse, n_bins.kjuse, n_bins.kimin, n_bins.kjmin);

  //read kj
  if(ini.get_param<std::string>("kj_file", &file_names)!=0) exit(31);
  kj = read_cut_gsl_vector(file_names, n_bins.kjtot, n_bins.kjuse, n_bins.kjmin);

  //read W0j
  if(ini.get_param<std::string>("w0j_file", &file_names)!=0) exit(32);
  W0j = read_cut_gsl_vector(file_names, n_bins.kjtot, n_bins.kjuse, n_bins.kjmin);

  //read whole G20i file
  if(ini.get_param<std::string>("G2i_file", &file_names)!=0) exit(33);
  gsl_vector *temp_G = read_cut_gsl_vector(file_names, n_bins.kitot+1);
  G20 = gsl_vector_get(gsl_vector, 0);  //get the first element
  G2i = gsl_vector_alloc(n_bins.kiuse);
  //copy the interesting part of the vector
  gsl_vector_memcpy(G2i, &(gsl_vector_subvector(temp_G, n_bins.kimin+1,
          n_bins.kiuse).vector));
  gsl_vector_free(temp_G);
}
