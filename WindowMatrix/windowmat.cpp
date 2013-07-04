/*==========================================================================*/
/* Version: 1.20           date:07/10/10                                    */
/*                                                                          */
/* Author: Francesco Montesano, MPE, Garching                               */
/*                                                                          */
/* Purpose: given a sferical averaged window function in k-space |G(k)|^2   */
/* computes the matrix version of it W_{ij}: P_m(ki) = sum_{ij} P(kj) W_{ij}*/
/*                                                                          */
/*==========================================================================*/

#include "/data01/montefra/Common_code/com_included_pack.h"
#include "gsl/gsl_errno.h"
#include "gsl/gsl_integration.h"
#include "gsl/gsl_spline.h"
#include "gsl/gsl_multifit.h"

#define EPS 3.0e-11

int n_outfiles = 4;   //number of output file names

struct in_data{
  double x, y;    // contain the input
};

struct par{  // structure containing the parameters for the integration
  gsl_spline *spl;   // declaration of the interpolation scheme
  gsl_interp_accel *acc;   // declaration of the accelarator
  double ki, kj;   // ki and kj 
  double kmin, kmax;   //miminum and maximum k integration (and window function) limits
  double gf[4];  //gauss fit enabled: if gf[0] negative not enabled, if positive enabled. gf[1]=A, gf[2]=mu, gf[3]=sigma
};

/* function declaration */
void help(char name[], int bins, int gf);        // help function
void ofilenames(string *ofiles);   //build the output file names
int check_file(string file_name);
int input_data(string input_file, vector<in_data> &indata, double kmin, double kmax);  // read the input file and give back a vector
void gauleg(double x1, double x2, double x[], double w[], int n);  // Numerical recipes: Gauss-Legendre n-point quadrature formula
void gaussfit(double *x, double *y, double *gf);  //gaussian fit
double func(double t, void *params);    //function to integrate


/*==========================================================================*/
/*                                 MAIN                                     */
/*==========================================================================*/

int main(int argc, char* argv[])
{
  int i, j;   //loop integers
  string input[2];    //file with the window function and with the value of k in wich it will computed
  string *output;   //file of the window function and of the value of k where to perform the internal sum
  vector<in_data> ks;   //vector containing (temporarly) the spherical averaged window function and the values of k where to evaluate the l.h.s. of the power spectrum
  in_data temp;  //temporary structure
  int dimwin, dimks;   // dimension of the above vectors
  double *k, *G;     //save winndow to two arrays for the win for the spline
  double *x, *y;     //temporary arrays for the gaussian fit
  int N_kj = 250;  //number of kj where to evaluate W_{ij}
  double exponent = 2.;   //kj^exponent in the integration.
  int gf = 10.;  //if positive enable gaussian fit
  double kjmin=-1, kjmax=-1;  // minimum and maximum values for the values of kj: if negative the minimum and maximum values from the window function used
  double kimin=-1, kimax=-1;  // minimum and maximum values for the values of ki from the "input[1]" file: if negative whole range
  double kwmin=-1, kwmax=-1;  // minimum and maximum values for the values of k for the window function and integration limits: if negative whole range considered. If kwmin negative gaussian fit done to the first bins of the window function to get G^2(0)
  double *kj, *w;   //values of kj and corresponding weights as given from gauled
  double **Wij;   //matrix representation of the window function |G(k)|^2
  double N=0.;  //normalization of each line of the window matrix
  par pars;    //parameters for the integration
  gsl_function f;  // gsl function 
  gsl_integration_workspace *ws;  // gsl integration workspace
  int wssize = 10000;    // size of gsl integration workspace
  double result, abserr;   // result and error from the numerical integration
  ofstream out;   //output output stream to file
  
  gsl_error_handler_t * gsl_set_error_handler_off();

  if(argc == 2 && (strcmp(argv[1], "-h") == 0 || strcmp(argv[1], "--help") == 0) ) help(argv[0], N_kj, gf);

  if(argc < 4){
    cout << "At least the input, output file prefix needed" << endl;
    cout << "For more information type '" << argv[0] << " --help'" << endl;
    exit(1);
  }
  output = new string[n_outfiles];

  for(i=0; i<2; ++i) input[i] = argv[i+1];  // get input file names
  output[0] = argv[3];  // get output file prefix
  ofilenames(output);   //construct output file names
  
  for(i=0; i<2; ++i) if(check_file(input[i]) == 0){
      cout << "Input file " << input[i] << " doesn't exist."  << endl;
      exit(2);
    }
  for(i=0; i<4; ++i) if(check_file(output[i]) == 1){
      cout << "Output file " << output[i] << " already exists."  << endl;
      exit(3);
    }

  if(argc>=5) N_kj = atoi(argv[4]);     // get the number of bins for the kj
  if(argc>=6) gf = atof(argv[5]);      // get if gaussian fit must be done
  if(argc>=7) kjmin = atof(argv[6]);   // get the minimum kj
  if(argc>=8) kjmax = atof(argv[7]);   // get the maximum kj
  if(argc>=9) kimin = atof(argv[8]);   // get the minimum ki 
  if(argc>=10) kimax = atof(argv[9]);  // get the maximum ki 
  if(argc>=11) kwmin = atof(argv[10]); // get the minimum kw 
  if(argc>=12) kwmax = atof(argv[11]); // get the maximum kw 

  /* ----------------------- */

  dimwin = input_data(input[0], ks, kwmin, kwmax);  //read the window function temporarly to ks, save into arrays
  k = new double[dimwin];   
  G = new double[dimwin];   
  for(i=0; i<int(dimwin); ++i){
    k[i] = ks[i].x;
    G[i] = ks[i].y;
  }
  ks.clear();   //clear the window function vector
  temp.x=0.;   //set the element in the first position of dimks containting ki=0
  ks.push_back(temp); 
  dimks = input_data(input[1], ks, kimin, kimax);  //read the file with the output ki

  /* allocate the k_j and weight arrays and the Wij matrix */
  kj = new double[abs(N_kj)];   
  w = new double[abs(N_kj)];   
  Wij = new double*[dimks];
  for(i=0; i<int(dimks); ++i) Wij[i] = new double[abs(N_kj)];
  for(i=0; i<abs(N_kj); ++i){
    kj[i]=0.;
    w[i]=0.;
  }

  /* gauleg */
  if(kjmin<0.) kjmin = k[0];        //if negative set the minumum and maximum value of kj as the 
  if(kjmax<0.) kjmax = k[dimwin-1]; //extrema of k of the window function
  if(N_kj < 0){   //if negative do gauleg in logarithmically spaced bins
    kjmin = log(kjmin);
    kjmax = log(kjmax);
    exponent = 3.;  //if logarithmic integration extra kj needed
  }
  gauleg(kjmin, kjmax, kj, w, abs(N_kj)); //get the gauled kj and w
  if(N_kj < 0) for(i=0; i<abs(N_kj); ++i) kj[i] = exp(kj[i]);   //if negative reconvert kj to linear

  /* set integration */
  pars.kmin = k[0];        //if negative set the minumum and maximum value of kj as the 
  pars.kmax = k[dimwin-1]; //extrema of k of the window function
  if(gf < 0.) pars.gf[0] = k[1];           //if no gaussian fit to be done
  else{
    x = new double[gf];
    y = new double[gf];
    for(i=0; i<gf; ++i){
      x[i] = k[i];
      y[i] = G[i];
    }
    pars.gf[0] = gf;   //if gaussian enabled integrate from k=0
    gaussfit(x, y, pars.gf);
    pars.gf[0] = 0.;   //if gaussian enabled integrate from k=0

    delete [] x;   //deallocate
    delete [] y;
  }

  /* create a spline of the |G(k)|^2 */
  pars.spl = gsl_spline_alloc(gsl_interp_cspline, dimwin);
  if(gsl_spline_init(pars.spl, k, G, dimwin)){
    cerr << "Error in the initialization of the interpolation routine with the '" << gsl_spline_name(pars.spl) << "' interpolation scheme" << endl;
    exit(11);
  }
  pars.acc = gsl_interp_accel_alloc();

  ws = gsl_integration_workspace_alloc(wssize);  // allocate the gsl workspace

  double integr;
  /* compute Wi_{j} */
  for(i=0; i<dimks; ++i){
    pars.ki = ks[i].x;   //save ki for the integration
    for(j=0; j<abs(N_kj); ++j){
      pars.kj = kj[j];   //save kj for the integration

      f.function = &func;   // associate the function
      f.params = &pars;  //associate the parameters

      //integr = gsl_integration_qag(&f, -1., 1, 1.e-4, 1.e-4, wssize, 4, ws, &result, &abserr);

      if(gsl_integration_qag(&f, -1., 1, 0., 8.e-3, wssize, 3, ws, &result, &abserr)){
      //if(integr != 0){
        cerr << "error with the integration at the loop " << i << ", " << j << ", corresponding to ki=" << ks[i].x << " and kj=" << kj[j] << endl;
	exit(12);
      }

      Wij[i][j] = w[j] * pow(kj[j],exponent) * result / (4*M_PI*M_PI);  //W_{ij}
      N += Wij[i][j]; 
    }
    if(N>0.) for(j=0; j<abs(N_kj); ++j) Wij[i][j] /= N;
    N = 0.; 
  }
  gsl_integration_workspace_free(ws);   // free the gsl workspace

  /* print the output files  */
  out.open(output[0].c_str());  //file with the window function Wij i!=0
  out.setf(ios_base::scientific);
  out.precision(6);
  out.width(9);
  for(i=1; i<dimks; ++i){
    for(j=0; j<abs(N_kj); ++j) out << Wij[i][j] << "\t";
    out << endl;
  }
  out.close();

  out.open(output[1].c_str());  //file with the kj
  out.setf(ios_base::scientific);
  out.precision(6);
  out.width(9);
  for(i=0; i<abs(N_kj); ++i) out << kj[i] << endl;
  out.close();

  out.open(output[2].c_str());  //file with the W0j
  out.setf(ios_base::scientific);
  out.precision(6);
  out.width(9);
  for(i=0; i<abs(N_kj); ++i) out << Wij[0][i] << endl;
  out.close();
  
  out.open(output[3].c_str());  //file with the G^2([0,ki])
  out.setf(ios_base::scientific);
  out.precision(6);
  out.width(9);
  pars.kj = 0.;   //save kj for the integration
  for(i=0; i<dimks; ++i){
    pars.ki = ks[i].x;   //save ki for the integration
    out << func(0., &pars) << endl;
  }
  out.close();

  delete [] k;
  delete [] G;
  delete [] output;

  return 0;
}


/*==========================================================================*/
/*                               FUNCTIONS                                  */
/*==========================================================================*/

void ofilenames(string *ofiles)
/* Given the output file name prefix in ofiles[0] returns the output file names. 
The array of strings must be already allocated with the right number of elements.*/
{
  int i;  //loop integers
  ostringstream s1;    // string stream to build the input linear and 1loop PS file name
  string ends[] = {".Wij.dat", ".kj.dat", ".W0j.dat", ".G20i.dat"};

  for(i=n_outfiles-1; i>=0; --i){
    s1 << ofiles[0] << ends[i];   //create the file name
    ofiles[i] = s1.str();        //save the output file name
    s1.str("");   //empty s1
  }
}

/*==========================================================================*/
int check_file(string file_name)
/* Input: name of the file to be checked
   Output: 0 if the file doesn't exists; 1 if it exists*/
{
  ifstream check(file_name.c_str());
  if (check.is_open())
  {
    check.close();
    return(1);
  }
  else return(0);
}

/*==========================================================================*/
int input_data(string input_file, vector<in_data> &indata, double kmin, double kmax)
/* get the file name input_file, return a vector of structure in_data
only the first two out of three columns are read*/
{
  in_data read;
  double dummy;

  ifstream in(input_file.c_str());
  for(;;){
    in >> read.x >> read.y >> dummy;
    if (in.eof() == true) break;
    else if(kmin>=0. && read.x<kmin) continue;  //if k < kmin skip 
    else if(kmax>=0. && read.x > kmax) break;   //if k > kmax stop reading
    else indata.push_back(read);
  }
  in.close();

  return(indata.size());
}

/*==========================================================================*/
void gauleg(double x1, double x2, double x[], double w[], int n)
/*Given the lower and upper limits of integration x1 and x2, and given n, this routine returns
arrays x[1..n] and w[1..n] of length n, containing the abscissas and weights of the Gauss-
Legendre n-point quadrature formula.*/
{
  int m,j,i;
  double z1,z,xm,xl,pp,p3,p2,p1;

  m=(n+1)/2;
  xm=0.5*(x2+x1);
  xl=0.5*(x2-x1);
  for(i=1;i<=m;i++) {
    z=cos(3.141592654*(i-0.25)/(n+0.5));
    do{
      p1=1.0;
      p2=0.0;
      for(j=1;j<=n;j++) {
	p3=p2;
	p2=p1;
	p1=((2.0*j-1.0)*z*p2-(j-1.0)*p3)/j;
      }
      pp=n*(z*p1-p2)/(z*z-1.0);
      z1=z;
      z=z1-p1/pp;
    } while(fabs(z-z1) > EPS);
    x[i-1]=xm-xl*z;
    x[n-i]=xm+xl*z;
    w[i-1]=2.0*xl/((1.0-z*z)*pp*pp);
    w[n-i]=w[i-1];
  }
}

/*==========================================================================*/
void gaussfit(double *x, double *y, double *gf)
/* Given arrays x and y=f(x) of gf[0] elements, it performs a gaussian fit and returns
the amplitud A=gf[1], the mean mu=gf[2] and the stddev sigma=gf[3]*/
{
  int i;    //loop integers
  gsl_multifit_linear_workspace *w;
  gsl_vector *ly, *c;   //input log(y), output coefficients
  gsl_matrix *X, *cov;  //input matrix of x values, output covariance
  double chisq;    // chi^2 output
  
  w = gsl_multifit_linear_alloc(gf[0], 3);   //allocate the workspace for the multifit
  ly = gsl_vector_alloc(gf[0]);  //allocate the temporary matrix
  c = gsl_vector_alloc(3);  //allocate the temporary matrix
  X = gsl_matrix_alloc(gf[0],3);  //allocate the temporary matrix
  cov = gsl_matrix_alloc(3,3);  //allocate the temporary matrix

  for(i=0; i<gf[0]; ++i) gsl_vector_set(ly, i, log(y[i]));  //save log(y) gsl vector
  for(i=0; i<gf[0]; ++i){   //save x^2, x and 1 in the matrix
    gsl_matrix_set(X, i, 0, x[i]*x[i]);
    gsl_matrix_set(X, i, 1, x[i]);
    gsl_matrix_set(X, i, 2, 1.);
  }

  if(gsl_multifit_linear(X, ly, c, cov, &chisq, w)){
    cerr << "Polinomial fit failed" << endl;
    exit(30);
  }

  gf[3] = sqrt(-1./(2.*gsl_vector_get(c,0)));    //sigma
  gf[2] = gsl_vector_get(c,1) * gf[3]*gf[3];     //mu
  gf[1] = exp(gsl_vector_get(c,2) + gf[2]*gf[2]/(gf[3]*gf[3]));   //amplitude

  gsl_vector_free(ly);   //deallocate gsl vectors and matrixes
  gsl_vector_free(c);
  gsl_matrix_free(X);
  gsl_matrix_free(cov);
  gsl_multifit_linear_free(w);   //deallocate the workspace
}

/*==========================================================================*/
double func(double t, void *params)
/* function to be integrated */
{
  par pars = *(par *) params;   //copy params to struct par
  double y=0.;    //output 
  double kk;     //|ki-kj|

  kk = sqrt(pars.ki*pars.ki + pars.kj*pars.kj - 2*t*pars.ki*pars.kj);  //compute |ki-kj|

  if(kk>=pars.kmin && kk<=pars.kmax){
    if(gsl_spline_eval_e(pars.spl, kk, pars.acc, &y)){
      cerr << "Error to spline the function |G(k)|^2 in k="<< kk << endl;
      exit(20);
    }
  }
  else if(kk>=pars.gf[0] && kk <=pars.kmin) y = pars.gf[1]*exp(- pow(kk-pars.gf[2],2.) / (2.*pars.gf[3]*pars.gf[3]) );
  return y;  // return the value of the function
}

/*==========================================================================*/
/* Help function: instrusctions to run the code: just type ./name --help    */
/* and an help file will be printed on screen                               */
/*==========================================================================*/

void help(char name[], int bins, int gf)
{
  cout << endl;
  cout << "/*==========================================================================*/" << endl;
  cout << "/* Version: 1.20           date:07/10/10                                    */" << endl;
  cout << "/*                                                                          */" << endl;
  cout << "/* Author: Francesco Montesano, MPE, Garching                               */" << endl;
  cout << "/*                                                                          */" << endl;
  cout << "/* Purpose: given a sferical averaged window function in k-space |G(k)|^2   */" << endl;
  cout << "/* computes the matrix version of it W_{ij}:                                */" << endl;
  cout << "/* P_m(ki) = sum_{ij} P(kj) W_{ij} + c |G(ki)|^2                            */" << endl;
  cout << "/*==========================================================================*/" << endl;
  cout << endl << endl;

  cout << "This program reads the spherical averaged window function and a list of ki values" << endl;
  cout << "and prints the window matrix, the values of kj, of W(0,kj) and of |G([0,ki])|^2" << endl;
  cout << "Usage" << endl;
  cout << name << " win_infile ki_infile outfile_prefix (N_kj gf kjmin kjmax kimin kimax kwmin kwmax)" << endl << endl;

  cout << "Parameters" << endl << "----------" << endl;
  cout << "win_infile: string" << endl << "  File name of the input window function G^2(k) (three columns: k, G^2(k), unused)." << endl;
  cout << "ki_infile: string" << endl << "  File name containing the value of ki (three columns: ki, unused, unused)." << endl;
  cout << "outfile_prefix: string" << endl << "  Prefix of the output file names. They will they are:" << endl;
  cout << "  outfile_prefix.Wij.dat, outfile_prefix.kj.dat, outfile_prefix.W0j.dat and outfile_prefix.G20i.dat." << endl;
  cout << "N_kj: integer" << endl << "  Number kj values used. If positive linear integration, if negative logarithmic. Default: " << bins << "." <<endl;
  cout << "gf: integer" << endl << "  if positive enable gaussian fit of the first 'gf' bins of the window function" << endl;
  cout << "  and the integration from k=0. Default " << gf << "." << endl;
  cout << "kjmin, kjmax: double" << endl << "  minimum and maximum values of kj. If negative, the extrema of the input window function used. Default kjmin=-1, kjmax=-1" << endl;
  cout << "kimin, kimax: double" << endl << "  minimum and maximum values of ki. If negative, all the ki used. Default kimin=-1, kimax=-1" << endl;
  cout << "kwmin, kwmax: double" << endl << "  minimum and maximum values for the integration. If negative, a gaussian fit to the first bins of the window function will be performed and the maximum set to the maximum of the input window function. Default kwmin=-1, kwmax=-1" << endl;

  exit(0);
}
