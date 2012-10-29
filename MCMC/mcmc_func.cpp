/*==========================================================================*/
/* Version: 3.00           date:17/01/12                                    */
/*                                                                          */
/* Author: Francesco Montesano, MPE, Garching                               */
/*                                                                          */
/* Purpose : function definitions for the main in 'mcmc.cpp'                */
/*==========================================================================*/

#include "mcmc.h"


/*==========================================================================*/
/*Function: read data from file:                                            */
/* input :name file and int 'which' to know what kind of file is load:      */
/* 0) linear P(k): k, P_{lin}(k) and first ordert PT contribution:          */
/* 1) input data: k, P(k), var                                              */
/* 2) covariance matrix: only the 5h column of interest                     */
/* 3) kj (only one column): saved in two columns                            */
/* k, P_{1loop}(k)                                                          */
/*==========================================================================*/
int input_data(std::string input_file, std::vector<in_data> &indata, int which)
{
  in_data read;
  double dummy;

  ifstream in(input_file.c_str());
  for(;;){
    if(which==0) in >> read.x >> read.y;   // read the theoretical expectation values
    else if(which==1) in >> read.x >> read.y >> read.z;  // it reads the third column in the data file
    else if(which==2) in >> dummy >> dummy >> dummy >> dummy >> read.x;  // it reads the covariance matrix
    else if(which==3){  //read the file with the kj
      in >> read.x; 
      read.y=read.x;
    }
    if(in.eof() == true) break;
    else indata.push_back(read);
  }
  in.close();

  return(indata.size());
}

/*==========================================================================*/
/* Reads the file containing a window matrix of size sizei*sizej and gives  */
/* back a nki*nkj gsl_matrix starting from the kimin and kjmin bins         */
/*==========================================================================*/

gsl_matrix * winmat(std::string winfile, int kimin, int kjmin, int nki, int nkj, int sizei, int sizej)
{
  gsl_matrix *temp;  //temporary matrix to read the input matrix
  FILE *f;   //file unit
  gsl_matrix_view vwin;   //gsl matrix view of a subsample of the full matrix
  gsl_matrix *win; //gsl matrix with the window function

  temp = gsl_matrix_alloc(sizei,sizej);  //allocate the temporary matrix
  win = gsl_matrix_alloc(nki,nkj);  //allocate the window matrix

  f = fopen(winfile.c_str(), "r");
  if(gsl_matrix_fscanf(f, temp) !=0){
    std::cerr << "The file " << winfile << "containing the matrix of the window function can't be read properly." << std::endl;
    exit(40);
  }
  fclose(f);
  vwin = gsl_matrix_submatrix(temp, kimin, kjmin, nki, nkj);   //take a view of the ki*kj submatrix
  gsl_matrix_memcpy(win, &(vwin.matrix));   //copy the part of the matrix into the window matrix

  gsl_matrix_free(temp);   //free the memory of the temporary matrix

  return win;
}

/*==========================================================================*/
/* Gets a vector in_data and save the second column it into a gsl vector    */
/* where the first columns kmin < k < kmax. The gsl vector is returned.     */
/* x_or_y = 0 indata.x copied, x_or_y = 1 indata.y copied                   */
/* if mod_kmin=true, the number of discarted bins for k < kmin is           */
/* substitute to kmin                                                       */
/*==========================================================================*/
gsl_vector * vec2gslvec(std::vector<in_data> &indata, double *kmin, double kmax, int x_or_y, bool mod_kmin)
{
  int i;   //loop integers
  int min=0, n=0;   //mimimun index and number of elements in the output vector
  gsl_vector * v;   //gsl vector 

  if(kmax == -1. && *kmin ==-1.) n = indata.size();   //if kmin==-1 and kmax==-1, all the indata kept
  else for(i=0; i<int(indata.size()); ++i){
      if(indata[i].x < *kmin) ++min;
      else ++n;
      if(indata[i].x > kmax) break; //otherwise vector check for the number of k<kmax
  }
  v = gsl_vector_alloc(n);   //allocate the vector gslvec with the right size of elements
  if(x_or_y == 0) for(i=0; i<n; ++i) gsl_vector_set(v, i, indata[min+i].x);  //set the vector from indata.x
  else for(i=0; i<n; ++i) gsl_vector_set(v, i, indata[min+i].y);  //set the vector from indata.y

  if(mod_kmin == true) *kmin = min;   //subsitute kmin with the number of bins with k<kmin if required
  return v;   //return the number of elements in the gsl_vector
}

/*==========================================================================*/
/* Function: take the covariance matrix in a vector, square it, discard the */
/* first mininv bins, cut such that it is diminv*diminv. Then invert it:    */
/* if diag=true as diagonal; if diag=false using the LU decomposition       */
/*==========================================================================*/
gsl_matrix * invert(std::vector<in_data> &indata, int diminv, int mininv, bool diag)
{
  int i, j;    // loop integers
  int dvec;   // size of the in_data vector
  gsl_permutation *p; // permutation for LU decomposition
  gsl_matrix *inv, *temp;   // gsl matrix required by the inversion algorithms
  int signum;    // required by LU decomposition

  dvec = int(sqrt(indata.size()));   //save the size of indata
  temp = gsl_matrix_alloc(diminv, diminv);   //allocate a temporary matrix
  inv = gsl_matrix_alloc(diminv, diminv);   //allocate the output inverse convariance matrix
  for(i=0; i<diminv; ++i) for(j=0; j<diminv; ++j) gsl_matrix_set(temp, i, j, indata[dvec * (mininv+i) + (mininv+j)].x); //square the matrix considering only the first diminv bins in each dimension

  if(diag == true){ //if the matrix is to be considered as diagonal
    gsl_matrix_set_zero(inv);   //set to 0 the inverse matrix
    for(i=0; i<diminv; ++i) gsl_matrix_set(inv, i, i, 1./gsl_matrix_get(temp, i, i));   //invert the diagonal of temp into inv
  }
  else{
    p = gsl_permutation_alloc(diminv);    // allocate permutation
    if(gsl_linalg_LU_decomp(temp, p, &signum)!= 0){   // LU decomposition
      std::cerr << "LU decomposition failed" << std::endl;
      exit(11);
    }
    if(gsl_linalg_LU_invert(temp, p, inv)!= 0){    // inversion
      std::cerr << "Inversion failed" << std::endl;
      exit(12);
    }
    gsl_permutation_free(p);
  }
  gsl_matrix_free(temp);

  return inv;
}

/*==========================================================================*/
/* read file 'fname', and save the first two columns into a gsl spline      */
/* object that is returned                                                  */
/*==========================================================================*/
gsl_spline * th2spl(std::string fname)
{
  std::vector<double> x,y;  //temporary vector to read the file
  gsl_spline * spl;   //output spline
  
  x.push_back(0);   //add point P(k=0)=0 to the input file: avoid errors when the spline is evaluated
  y.push_back(0);   //and shouldn't influence the output 
  
  if( readfirst2columns(fname, 2, x, y) == 1 ) 
    throw inputexception(fname);   //read the input file and save in two vectors 

  x.push_back(x[x.size()-1] *100. );   //add point P(k=very large)=0 to the input file: avoid errors when the spline is evaluated
  y.push_back(0);   //and shouldn't influence the output 

//  vector<double>::iterator iter = x.begin();  //get the interator of the 
//  xit = x.insert( xit, x[0]*0.99 );

  spl = gsl_spline_alloc(gsl_interp_linear, x.size());
  if(gsl_spline_init(spl, &x[0], &y[0], x.size()) !=0){
    std::cerr << "Initialization of the spline " << gsl_spline_name(spl);
    std::cerr << " for the data in file " << fname << " failed" << std::endl;
    exit(41);
  }

  x.clear();
  y.clear();
  return spl;  //return the spline
}

/*==========================================================================*/
/* Function: calculate the likelihood given the parameters stored in the    */
/* struc 'parameters'                                                       */
/*==========================================================================*/

double likelihood(double *par, pmc *p22)
{
  double tCt, like; // t^{t}inv t product

  /* the array that contain the theory power given the parameters 'par'*/
  for(size_t i=0; i<kj->size; ++i) gsl_vector_set( theory, i,
      p_obj(gsl_vector_get(kj, i), par, p22) );

#ifdef WINMAT    
  /* if window matrix is used */
  /* convolve the theory power spectrum with the window function */
  if(gsl_blas_dgemv(CblasNoTrans, 1., window, theory, 0., wtheory) !=0) exit(20);
  /* integral constraint: add a rescaled version of the window function to wtheory */
  double c;  //constant used to rescale the window function
  if(gsl_blas_ddot(W0j, theory, &c)) exit(21);  // c = sum_j W0j P(kj)
  if(gsl_blas_daxpy(-c/G20, G2, wtheory)) exit(22);  // wtheory = wtheory - c*G^2(ki)/G^2(0)
#else
  /* if window matrix is *NOT* used */
  if(gsl_vector_memcpy(wtheory, theory) != 0) exit(20);
#endif

  // set to 0 the the temporary array
  gsl_vector_set_zero(tempv);
  if( par[4] < 0 ){   //if the value of the bias is negative
    /* marginalization on the amplitude: Lewis 2002, F2 */
    gsl_matrix_set_zero(tempm);
    /* t^{T}invcov t */
    if(gsl_blas_dgemv(CblasNoTrans, 1., invcov, wtheory, 0., tempv)!= 0) exit(23);  // tempv=invcov t
    if(gsl_blas_ddot(wtheory, tempv, &tCt)!= 0) exit(24);  // tCt=t^{T}tempv
    /* invcov - (invcov tt^{T} invcov)/tCt = invcov - (invcov t (invcov^{T}
     * t)^{T})/tCt */
  // set to 0 the the temporary array
    gsl_vector_set_zero(tempv);
    if(gsl_blas_dgemv(CblasTrans, 1., invcov, wtheory, 0., tempv)!= 0)
      exit(25);  // tempv=invcov^{T} t
    if(gsl_blas_dger(1., wtheory, tempv, tempm)!= 0) exit(26);  // tempm=t tempv^{T}
    if(gsl_matrix_memcpy(tempm2, invcov)!= 0) exit(27);   // copy invcov to tempm2
    if(gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, -1./tCt, invcov, tempm, 1.,
	  tempm2)!= 0) exit(28);  // tempm2=invcov - (invcov tempm)/tCt
    /* d^{T}tempm2 d */
    if(gsl_blas_dgemv(CblasNoTrans, 1., tempm2, data, 0., tempv)!= 0) exit(29);  // tempv=tempm2 d
    if(gsl_blas_ddot(data, tempv, &like)!= 0) exit(30);  // like=d^{T}tempv
    like += log(tCt);   //get the total likelihood
  }
  else{  //otherwise the bias is included in the likelihood determination
    if( gsl_vector_sub( wtheory, data ) != 0 ) exit(32);  //theory = theory - data
    if( gsl_blas_dgemv(CblasNoTrans, 1., invcov, wtheory, 0., tempv)!= 0)
      exit(23);  // tempv=invcov (theory-data)
    if( gsl_blas_ddot(wtheory, tempv, &like)!= 0 ) exit(24);  // like=(theory-data)^{T}tempv
  }

  return(like);
}

/*==========================================================================*/
/* Given the index i, such that k=data[i].x and the set of parameters       */
/* computs the P_{obj}(k)                                                   */
/*==========================================================================*/

double p_obj(double k, double *par, pmc *p22){
  double theo;   // contains the varius pieces of the calculation (too long for a single line)

  theo = gsl_spline_eval(spllin, par[2]*k, acclin); // * exp(-k*k/(par[0]*par[0]));  //exp(-(k/kstar)^2) P_lin(alpha*k)
  theo += par[1] * p22->pMC(par[2]*k);                      //+ A_MC P_1loop(alpha*k)
  if( par[4] > 0 ) theo *= par[4] * par[4];  //if the bias is a parameter
  theo *= rsd( k, par, cosmo);     //rds 
  theo += par[3];   //add noise
  return( theo );                     
}

/*==========================================================================
 * angular average redshift space distortions: analytical expression
 * from Ariel.
 * it includes the Kaiser boost factor and finger of god.
 * Parameters
 * ----------
 * k: double
 *   wave number 
 * par: array of doubles
 *   parameters for the step of the mcmc
 * cosmo: cosmology struct
 *   cosmology related parameters. In this case f = ln(D)/ln(a)   
 * output
 * ------
 * rsd: double
 *   value of the correction for the given set of parameters and k
 *==========================================================================*/
double rsd( double k, double *par, cosmology c){
  double rsd;   //value to return
  if( c.f <= 0 ) rsd = 1;  //if no redshift space distortion wanted
  else{  //if required
    double beta = c.f/par[4];
    double x = k * c.f * par[5]; 
    rsd = atan(x)/x * ( pow(beta+x*x, 2) - 4.*beta*beta) + 2.*beta*beta;
    rsd += pow(beta-x*x, 2) / (x*x+1.); 
    rsd /= 2.*pow( x, 4 );
  }
  return( rsd );
}


/*==========================================================================*/
/* prints the parameter names and returns the short name to be used in the  */
/* *.ini file                                                               */
/* WARNING: update it when changing the number or the name of the parameters*/
/*==========================================================================*/
int paramnames(std::string ofile, std::vector<std::string> &sname){
  int i;   //loop integers
  int n_params = 6;
  std::vector<std::string> lname;   //vectors of strings that contain long the name of the parameters in latex format

  sname.push_back("k_star");
  lname.push_back("k_\\star [h/Mpc]");
  sname.push_back("amc");
  lname.push_back("A_{\\mathrm{MC}}");
  sname.push_back("alpha");
  lname.push_back("\\alpha");
  sname.push_back("noise");
  lname.push_back("n [(Mpc/h)^3]");
  sname.push_back("bias"); 
  lname.push_back("b");
  sname.push_back("sigmav");
  lname.push_back("sigma_{\\mathrm{v}} [Mpc/h]");

  ofstream out(ofile.c_str());
  for(i=0; i<n_params; ++i) out << sname[i] << "\t" << lname[i] << std::endl;
  out.close();

  lname.clear();  //clear the long names, not used anymore

  return n_params;
}

/*==========================================================================*/
/* read the values of the starting point, of upper and lower limit          */
/* and of the step size of the parameters from 'inifile'. The parameter     */
/* names are stored in 'sname'                                              */
/*==========================================================================*/
void readini(std::string inifile, std::vector<std::string> sname, 
    double *params, double *upper, double *lower, double *sigma)
{
  size_t i,j;  //loop integers
  std::vector<std::string> inis;   //vector containing the parameters initial guess, boundaries and sigmas
  std::string par;   //which parameter is wanted and the numbers associated to it

  inis = readlines(inifile);  //get a vector of strings with only the important lines of the ini file 

  for(i=0; i<sname.size(); ++i){   //loop through the hard coded parameters
    par = "param["+sname[i]+"]";  //create the path to search
    for(j=0; j<inis.size(); ++j){  //loop through the inifile lines
      if(inis[j].compare(0,par.size(), par) == 0){  //if a parameter found
	stringstream iniparse;  //parse the four numerical values
	iniparse << inis[j].substr(inis[j].find("=")+1);   //get everything after the '='
	iniparse >> params[i];   //copy the initial guess
	iniparse >> lower[i];    //copy the lower bound
	iniparse >> upper[i];    //copy the upper bound
	iniparse >> sigma[i];    //copy the step size
	iniparse.str("");  //clear iniparse
	break;  //if the parameter found, go to the next one
      }
    }
  }
}

/*==========================================================================
 * Substitute the a value in the input array with the desired one
 * Parameters
 * ----------
 * array: pointer to double 
 *   array where to modify
 * n: size_t
 *   size of array
 * i: size_t
 *   position to substitute
 * val: double
 *   array[i] = val
 * output
 * ------
 * oarray: pointer to double
 *   pointer to the array with the substituted value
 *==========================================================================*/
double * substitute_value( double *array, size_t n, size_t i, double val ){
  double * oarray = new double[n];  //allocate the output array
  for( size_t j=0; j<n; ++j) oarray[j] = array[j];   //make a local copy of the array
  oarray[i] = val;   // substitute the value 
  return( oarray ); // return the pointer to the array
}

/*==========================================================================*/
/* deallocate common variables                                              */
/*==========================================================================*/
void dealloc_comm()
{
  /* deallocate the matrices and vectors for the likelihood */
  gsl_vector_free(theory);
  gsl_vector_free(wtheory);
  gsl_vector_free(tempv);
  gsl_matrix_free(tempm);
  gsl_matrix_free(tempm2);
  /* deallocate the matrices and vectors storing the power spectra, window and covariance matrices */
  gsl_vector_free(data);   //P(k)
  gsl_matrix_free(invcov); //inversed of the covariance matrix
  gsl_vector_free(kj);     //kj
#ifdef WINMAT    
  gsl_vector_free(G2);     //G2(ki)
  gsl_vector_free(W0j);    //W0j
  gsl_matrix_free(window); //window function
#endif
  gsl_spline_free(spllin);   //gsl spline of the linear and the 1loop ps
  gsl_interp_accel_free(acclin);   //gsl accellation of the linear and the 1loop ps
}

