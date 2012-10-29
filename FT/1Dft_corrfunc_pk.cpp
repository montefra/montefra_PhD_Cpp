/*==========================================================================*/
/* version: 0.01, 13/07/2011                                                */       
/* author: Francesco Montesano (montefra@mpe.mpg.de)                        */
/*                                                                          */
/* Purpose: compute the 1D fourier transform of the power spectrum into the */
/* correlation function or viceversa using gsl integration routines         */
/* P(k) = 4\pi \int_0_\infty dr r^2 xi(r) j_0(rk)                           */
/* xi(r) = 1/(2\pi^2) \int_0_\infty dk k^2 P(k) j_0(rk)                     */
/*==========================================================================*/

#include "1Dft_corrfunc_pk.h"

/*==========================================================================*/
/* Fourier transform the power spectrum into the correlation function       */
/* 'k' and 'pk' are the input wavenumber and P(k) with dimention 'n' or 'nk'*/
/* 'xi' contains the output correlation function.                           */
/* all pointers have to point to already allocated vectors                  */
/*==========================================================================*/
double* pk2corfunc(double* k, double* pk, double* xi, size_t n)
//xi(r) must have 'n' elements and is evaluated in r=\pi/k
//the pointer to the values of r returned
{
  size_t i;   //loop variable
  double* r;  //the array of x where to evaluate the correlation function

  try {r = new double[n];}    //allocate r
  catch(bad_alloc& e) {cerr << "Error allocating the array for 'r' in " << __FILE__ << " at line " << __LINE__-1 << endl;}
  for(i=0; i<n; ++i) r[i] = M_PI / k[n-i-1];   //fill the input x array

  ftQawf(k, pk, n, r, xi, n);

  for(i=0; i<n; ++i) xi[i] /= 2.*M_PI*M_PI;   //set the normalisation
  return r;   //return the pointer to r
}
void pk2corfunc(double* k, double* pk, size_t nk, double* r, double* xi, size_t nr)
//'xi' is evaluated at the given 'x' scales; 'xi' and 'x' must have 'nx' elements
{
  size_t i;   //loop variable

  ftQawf(k, pk, nk, r, xi, nr);

  for(i=0; i<nr; ++i) xi[i] /= 2.*M_PI*M_PI;   //set the normalisation
}

/*==========================================================================*/
/* Fourier transform the correlation function into the power spectrum       */
/* 'r' and 'xi' are the input scale and xi(r) with dimension 'n' or 'nr'    */
/* 'pk' contains the output correlation function.                           */
/* all pointers have to point to already allocated vectors                  */
/*==========================================================================*/
double* corfunc2pk(double* r, double* xi, double* pk, size_t n)
//'pk' must have 'n' elements and is evaluated in k=\pi/x
//the pointer to the values of k returned
{
  size_t i;   //loop variable
  double* k = new double[n];  //allocate the array of k where to evaluate the power spectrum
  
  for(i=0; i<n; ++i) k[i] = M_PI / r[n-i-1];   //fill the input k array

  ftQawf(r, xi, n, k, pk, n);

  for(i=0; i<n; ++i) pk[i] *= 4.*M_PI;   //set the normalisation
  return k;   //return the pointer to k
}
void corfunc2pk(double* r, double* xi, size_t nr, double* k, double* pk, size_t nk)
//'Pk' is evaluated at the given 'k' wavenumbers; 'Pk' and 'k' must have 'nk' elements
{
  size_t i;   //loop variable

  ftQawf(r, xi, nr, k, pk, nk);

  for(i=0; i<nk; ++i) pk[i] *= 4.*M_PI;   //set the normalisation
}

/*==========================================================================*/
/* compute f(x) = int_0^\infty dk k^2 fk(k) j_0(k*x)                        */
/* j_0(k*x) is the spherical bessel function of order 0                     */
/* it uses the gsl function 'gsl_integration_qawo' to compute the integral  */
/* f(x) = int_k[0]^k[nk-1] dk k/x fk(k) sin(k*x)                            */
/* 'nx' is the size of 'x' and 'fx' and 'nk' is the size of 'k' and 'fk'    */
/* 'fx' needs to be already allocated                                       */
/*==========================================================================*/
void ftQawf(double* k, double* fk, size_t nk, double* x, double* fx, size_t nx)
{
  size_t i;  //loop integers
  size_t sizew = 500, sizetable=10;      //size of the gsl integration workspace 'w' and number of levels in the qawo table
  int integral_status=0;       //output error message from the qawo integration
  int error_gsl=0;   //error status of gsl
  double l = k[nk-1]-k[0];   //range of the integration
  gslfparams fparams;    //struct with the parameters to be pass to the gsl function
  double* kfk;    //array containing k*fk. Needed for the interpolation
  double epsabs=1.e-5, epsrel=1.e-4;   //set the default precision
  double integrationerror;    //error returned from the integration routine

  gsl_function f;    // gsl function declaration
  gsl_error_handler_t* default_handler;    // declaration of the default gsl error handler
  gsl_integration_workspace* w;     //gsl integration workspace
  gsl_integration_qawo_table* qawo_table;  //table for the qawo routine

  /* deactivate standard error handler */
  default_handler = gsl_set_error_handler_off();   //turn off the error handler as save it in 'default_handler'

  /* allocate and initialise spline and allocate accelerator */
  fparams.spl = gsl_spline_alloc(gsl_interp_akima, nk);   //alloc a gsl_spline that will go into 'function'
  if(fparams.spl == 0){
    cerr << __FILE__ << ": " << __LINE__-2 << ". ERROR: the allocation of the spline object failed" << endl;
    exit(10);
  }
  try {kfk = new double[nk];}    //allocate kfk
  catch(bad_alloc& e) {cerr << "Error allocating the array for k*f(k) in " << __FILE__ << " at line " << __LINE__-1 << endl;}
  for(i=0; i<nk; ++i) kfk[i] = k[i]*fk[i];    //compute k*f(k) and fill 'kfk'
  error_gsl = gsl_spline_init(fparams.spl, k, kfk, nk);
  if( error_gsl != 0 ){
      cerr << __FILE__ << ": " << __LINE__-1 << ". ";
      if( error_gsl == 4 ) cerr << "ERROR: x values must be monotonically increasing" << endl;
      else cerr << "ERROR: initialisation of the gsl spline failed" << endl;
      exit(20);
  }
  delete[] kfk;     //deallocate kfk
  fparams.acc = gsl_interp_accel_alloc();   //allocate the accelerator for the interpolation
  if(fparams.spl == 0){
    cerr << __FILE__ << ": " << __LINE__-2 << ". ERROR: the allocation of the interpolation accelerator failed" << endl;
    exit(11);
  }
  /* pass the lower and upper limit of k to the function parameters */
  fparams.a = k[0];
  fparams.b = k[nk-1];

  /* assign the function and the parameters to the gsl function */
  f.function = &function;   //assign the function 'function' to integrate to the gsl function 'f'
  f.params = &fparams;       //assigne the parameter struct to the gsl function 'f'

  w = alloc_integration_workspace(sizew, __FILE__, __LINE__);

  qawo_table = alloc_qawo_table(x[0], l, GSL_INTEG_SINE, sizetable, __FILE__, __LINE__);  //allocate the qawo table and check for errors

  /* integration */
  for(i=0; i<nx; ++i){   //loop through the value of x and compute f(x)

    if( gsl_integration_qawo_table_set(qawo_table, x[i], l, GSL_INTEG_SINE) != 0 ){    //changes the parameter x in sin(x*k)
      cerr << __FILE__ << ": " << __LINE__-1 << ". ";
      cerr << "The change of the frequency in the qawo table failed for i=" << i << " and x=" << x[i] << endl;
      exit(21);
    }

    for(;;){
      //integral_status = gsl_integration_qag(&f, k[0], k[nk-1], epsabs, epsrel, sizew, GSL_INTEG_GAUSS61, w, &fx[i], &integrationerror);
      integral_status = gsl_integration_qawo(&f, k[0], epsabs, epsrel, sizew, w, qawo_table, &fx[i], &integrationerror);

      if(integral_status == 0) break;   //if everything is fine get out of the infinit loop
      else if(integral_status == GSL_ETABLE){  //if the number of levels is insufficient for the requested accuracy
	gsl_integration_qawo_table_free(qawo_table); //frees all the memory associated with the qawo table 'qawo_table'
	sizetable += 2;   //increase the size of the qawo table
	qawo_table = alloc_qawo_table(x[i], l, GSL_INTEG_SINE, sizetable, __FILE__, __LINE__);  //reallocate the qawo table 
	continue;   //reevaluate the integral
      }
      else if(integral_status == GSL_EMAXITER){    //the maximum number of subdivisions was exceeded.
	gsl_integration_workspace_free(w);         //free the memory from the gsl integration workspace 'w'
	sizew *= 2;   //double the size of the workspace
	w = alloc_integration_workspace(sizew, __FILE__, __LINE__);  //reallocate 'w'
	continue;   //reevaluate the integral
      }
      else if(integral_status == GSL_EROUND){  // cannot reach tolerance because of roundoff error, or roundoff error was detected in the extrapolation table.
	epsrel *= 10.;  //decrease the required relative precision
	continue;   //reevaluate the integral
      }
      else{
	cerr << "Error " << integral_status << "happenden while integrating. I don't know if can be dynamically solved" << endl;
	exit(66);
      }
    }
    epsrel=1.e-4;   //reset the required relative precision
    fx[i] /= x[i];   //divide by x the fourier transform
  }
  gsl_set_error_handler(default_handler);    //restore default error handler

  gsl_integration_qawo_table_free(qawo_table); //frees all the memory associated with the qawo table 'qawo_table'
  gsl_integration_workspace_free(w);         //free the memory from the gsl integration workspace 'w'
  gsl_interp_accel_free(fparams.acc);    //free the gsl spline accelerator
  gsl_spline_free(fparams.spl);   //free the gsl spline object
}

/*==========================================================================*/
/* function to be passed to the gsl integration routine                     */
/*==========================================================================*/
double function(double x, void* params)
{
  gslfparams fparams = *(gslfparams *) params;    //assign the parameters to a local variable
  double y=0.;  //output from the spline evaluation routine

  if(x>=fparams.a && x<=fparams.b){
    if( gsl_spline_eval_e(fparams.spl, x, fparams.acc, &y) != 0 ){
      cerr << __FILE__ << ": " << __LINE__-1 << ". ";
      cerr << "Error in interpolating f(x) in x=" << x << ". Make sure that it isn't outside the boundaries." << endl;
      exit(30);
    }
  }
  return y;
}

/*==========================================================================*/
/* misc functions                                                           */
/*==========================================================================*/

/* allocate a qawo table and check for errors                               */
/* for the defination of the first four quantities see gsl manual           */
/* the last two are name of the file and the number of the line where the   */
/* function is called: it's use for the error message                       */
/* can be given with the __FILE__ and __LINE__ preprocessor definitions     */
gsl_integration_qawo_table* alloc_qawo_table(double omega, double L, enum gsl_integration_qawo_enum sine,
                                             size_t n, string filename, int linenumber)
{
  gsl_integration_qawo_table* qawo_table;   //declare qawo table
  qawo_table = gsl_integration_qawo_table_alloc(omega, L, sine, n);   //allocate the qawo table as sin(x*k)
  if(qawo_table == 0){    //check allocation errors
    cerr << filename << ": " << linenumber << ". The allocation of the qawo table failed" << endl;
    exit(13);
  }
  return qawo_table;
}

/* allocate a integration workspace and check for errors                    */
/* for the defination of the first quantity see gsl manual                  */
/* the last two are name of the file and the number of the line where the   */
/* function is called: it's use for the error message                       */
/* can be given with the __FILE__ and __LINE__ preprocessor definitions     */
gsl_integration_workspace* alloc_integration_workspace(size_t n, string filename, int linenumber)
{
  gsl_integration_workspace* w;   //declare the integration workspace
  w = gsl_integration_workspace_alloc(n);    //allocate the integration workspace
  if(w == 0){
    cerr << filename << ": " << linenumber << ". The allocation of the integration workspace failed" << endl;
    exit(12);
  }
  return w;
}
